function [fftCfg, noiseCfg] = estNoiseFloor_pipe(fftCfg)
%ESTNOISEFLOOR_PIPE 逐距离门估计噪底功率, 标记强干扰距离门与弱chirp
%   输出: fftCfg.rxPowSum -- 各距离门跨chirp的通道平均功率 [512距离 × 512多普勒]
%         noiseCfg.floorPow -- 每个距离门的噪底功率(线性值)
%         noiseCfg.heavyFlag -- 强干扰距离门标记(分布拖尾过宽的门)
%         noiseCfg.weakChirpFlag -- 噪底突变距离门标记
%         noiseCfg.floorDbAvg -- 远端距离门平均噪底(dB), 全局参考电平
    % 取前512个距离bin, 对64个通道求平均功率谱: 距离-多普勒二维功率图
    fftCfg.rxPowSum = squeeze(mean(abs(fftCfg.data_sva(1 : 512, :, :)) .^ 2, 2));
    dopBinNum = size(fftCfg.rxPowSum, 2);
    noiseCfg.heavyFlag = zeros(dopBinNum, 1);
    % ---------- 逐距离门统计多普勒维功率分布, 用直方图中位数估计噪底 ----------
    for rngIdx = 1 : size(fftCfg.rxPowSum, 1)
        % 一个距离门内512个多普勒bin绝大多数是纯噪声, 目标只占少数bin,
        % 因此用分布的中位数(而非均值)作噪底, 可不受强目标影响
        logPow = log10(fftCfg.rxPowSum(rngIdx, :));
        % 溢出/异常尖峰限幅, 防止撑爆直方图动态范围
        spikeIdx = logPow > 1e3;
        if ~isempty(spikeIdx)
            logPow(spikeIdx) = 100;
            fftCfg.rxPowSum(rngIdx, spikeIdx) = 1e100;

        end

        pLo = min(logPow);
        pHi = max(logPow);

        % 在[最小,最大]对数功率范围内等宽划分256个直方图bin:
        % binScale为每bin宽度对应的对数功率跨度, binEdge为各bin边界
        binScale = (dopBinNum / 2 - 1) / (pHi - pLo + 1e-9);
        binEdge = pLo + (0 : dopBinNum / 2) / binScale;
        histCnt = zeros(1, dopBinNum/2);
        % 逐多普勒bin落入对应直方图bin计数
        for k = 1 : dopBinNum
            binPos = round((logPow(k) - pLo) * binScale) + 1;
            histCnt(binPos) = histCnt(binPos) + 1;
        end

        % 由累计分布CDF定位 中位数/Q1/Q3/最大值 所在bin
        histCdf = cumsum(histCnt);
        [medBin] = find(histCdf >= dopBinNum/2, 1);
        [q1Bin] = find(histCdf >= dopBinNum/4, 1);
        [maxBin] = find(histCdf >= dopBinNum, 1);
        [q3Bin] = find(histCdf >= dopBinNum*0.75, 1);
        if rngIdx < 0
            % 恒假分支(保留代码): 若启用则以Q3(上四分位)作噪底, 偏保守
            noiseCfg.floorPow(rngIdx, 1) = 10.^((binEdge(q3Bin) + binEdge(q3Bin + 1))/2);
        else
            % 实际生效: 取中位数bin边界中点, 换回线性功率
            noiseCfg.floorPow(rngIdx, 1) = 10.^((binEdge(medBin) + binEdge(medBin + 1))/2);
        end
        % 噪底上界: 最大值bin处功率(含目标/干扰的贡献)
        noiseCfg.floorPowFull(rngIdx, 1) = 10.^((binEdge(maxBin) + binEdge(maxBin + 1))/2);
        % 分布宽度: 中位数与Q1在对数功率轴上的间距, 反映门内干扰/杂波强度
        noiseCfg.spreadVal(rngIdx) = (medBin - q1Bin) / binScale;
        if((medBin - q1Bin) / binScale > 2.3/10)
            noiseCfg.heavyFlag(rngIdx) = 1;     % 拖尾过宽 -> 判为强干扰距离门
        end
    end

    % ---------- 标记噪底突变的距离门 ----------
    noiseCfg.weakChirpFlag = zeros(dopBinNum,1);

    % 相邻距离门噪底跳变(此处db()/2等效半功率轴)超过2.5的门,
    % 跳变起始点的后一个门被标记, 供CFAR环节区分处理
    floorJump = diff(db(noiseCfg.floorPow)/2);
    jumpPos = find(floorJump > 2.5);
    noiseCfg.weakChirpFlag(jumpPos + 1) = 1;

    % 取约3/4量程处(376~396)远端距离门的平均噪底(dB), 该处通常无目标, 作全局参考
    noiseCfg.floorDbAvg = 10*log10(mean(noiseCfg.floorPow(512 * 3 /4 - 10 : 512 * 3 /4 + 10)));

end
