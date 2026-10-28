function [fftCfg, noiseCfg] = estNoiseFloor(fftCfg)
%ESTNOISEFLOOR 噪底估计: rxPowSum(原功能, 含全部判断) + 逐 RX 扩展(仅噪声统计)
%   输入:
%       fftCfg.data_sva : 加窗处理后的 3D 频域数据, 复数 [rangeBin, RX, dopplerBin]
%   输出:
%   【原有功能: 对 rxPowSum(RX 平均功率谱)的噪声统计 + 全部判断, 与原代码一致】
%       fftCfg.rxPowSum        对 RX 平均的每 [距离,多普勒] 功率(线性) [512×dop]
%       noiseCfg.floorPow      rxPowSum 的 medBin 噪底(中位数法) [512×1]
%       noiseCfg.floorPowFull  rxPowSum 的 maxBin 噪底(最大尾部法) [512×1]
%       noiseCfg.spreadVal     分布宽度(强干扰判据) [512×1]
%       noiseCfg.heavyFlag     强干扰距离门标记(1/0) [512×1]
%       noiseCfg.weakChirpFlag 弱chirp/噪底突跳位置标记(1/0) [512×1]
%       noiseCfg.floorDbAvg    rxPowSum 噪底(dB)后段距离平均(标量)
%   【新增扩展: 逐 RX 通道噪声统计(只做噪声, 不做判断)】
%       noiseCfg.floorPowCh    各通道 medBin 噪底 [512×nRX]
%       noiseCfg.floorPowFullCh 各通道 maxBin 噪底 [512×nRX]
%       noiseCfg.figNoise      8 通道噪底图句柄(本函数内绘制, 每通道 4 条曲线)

    %------------------------------------------------------------
    % 1) rxPowSum: 对 RX 平均的每 [距离,多普勒] 功率(线性)
    %------------------------------------------------------------
    nRngUse = size(fftCfg.data_sva, 1) - 1;   % 有效距离 bin 数(= nV-1; 1024 实采样 → 512, 读数据不硬编码)
    fftCfg.rxPowSum = squeeze(mean(abs(fftCfg.data_sva(1 : nRngUse, :, :)) .^ 2, 2));

    dopBinNum = size(fftCfg.rxPowSum, 2);

    noiseCfg.heavyFlag = zeros(dopBinNum, 1);

    %------------------------------------------------------------
    % 2) 原有功能: 对 rxPowSum 逐距离门噪声统计 + 全部判断
    %------------------------------------------------------------
    for rngIdx = 1 : size(fftCfg.rxPowSum, 1)
        logPow = log10(fftCfg.rxPowSum(rngIdx, :));

        spikeIdx = logPow > 1e3;
        if ~isempty(spikeIdx)
            logPow(spikeIdx) = 100;
            fftCfg.rxPowSum(rngIdx, spikeIdx) = 1e100;
        end

        pLo = min(logPow);
        pHi = max(logPow);

        binScale = (dopBinNum / 2 - 1) / (pHi - pLo + 1e-9);
        binEdge  = pLo + (0 : dopBinNum / 2) / binScale;

        histCnt = zeros(1, dopBinNum/2);
        for k = 1 : dopBinNum
            binPos = round((logPow(k) - pLo) * binScale) + 1;
            binPos = min(max(binPos, 1), dopBinNum/2);   % 越界保护(不改分布)
            histCnt(binPos) = histCnt(binPos) + 1;
        end

        histCdf = cumsum(histCnt);
        [medBin] = find(histCdf >= dopBinNum/2, 1);
        [q1Bin]  = find(histCdf >= dopBinNum/4, 1);
        [maxBin] = find(histCdf >= dopBinNum,   1);
        [q3Bin]  = find(histCdf >= dopBinNum*0.75, 1);

        noiseCfg.floorPow(rngIdx, 1) = 10.^((binEdge(medBin) + binEdge(medBin + 1))/2);
        noiseCfg.floorPowFull(rngIdx, 1) = 10.^((binEdge(maxBin) + binEdge(maxBin + 1))/2);

        noiseCfg.spreadVal(rngIdx) = (medBin - q1Bin) / binScale;
        if((medBin - q1Bin) / binScale > 2.3/10)
            noiseCfg.heavyFlag(rngIdx) = 1;
        end
    end

    noiseCfg.weakChirpFlag = zeros(dopBinNum,1);
    floorJump = diff(db(noiseCfg.floorPow)/2);
    jumpPos = find(floorJump > 2.5);
    noiseCfg.weakChirpFlag(jumpPos + 1) = 1;

    noiseCfg.floorDbAvg = 10*log10(mean( ...
        noiseCfg.floorPow(round(nRngUse * 3/4) - 10 : round(nRngUse * 3/4) + 10) ...
    ));

    %------------------------------------------------------------
    % 3) 新增扩展: 逐 RX 通道噪声统计(只统计噪声, 不做任何判断)
    %    算法与上面完全相同, 数据换成"单通道功率谱"
    %------------------------------------------------------------
    nRx = size(fftCfg.data_sva, 2);
    for rxIdx = 1 : nRx
        rxPow = squeeze(abs(fftCfg.data_sva(1 : nRngUse, rxIdx, :)) .^ 2);  % 本通道 [nRngUse×dop]
        for rngIdx = 1 : size(rxPow, 1)
            logPow = log10(rxPow(rngIdx, :));

            spikeIdx = logPow > 1e3;
            if ~isempty(spikeIdx)
                logPow(spikeIdx) = 100;
                rxPow(rngIdx, spikeIdx) = 1e100;
            end

            pLo = min(logPow);
            pHi = max(logPow);

            binScale = (dopBinNum / 2 - 1) / (pHi - pLo + 1e-9);
            binEdge  = pLo + (0 : dopBinNum / 2) / binScale;

            histCnt = zeros(1, dopBinNum/2);
            for k = 1 : dopBinNum
                binPos = round((logPow(k) - pLo) * binScale) + 1;
                binPos = min(max(binPos, 1), dopBinNum/2);
                histCnt(binPos) = histCnt(binPos) + 1;
            end

            histCdf = cumsum(histCnt);
            [medBin] = find(histCdf >= dopBinNum/2, 1);
            [maxBin] = find(histCdf >= dopBinNum,   1);

            noiseCfg.floorPowCh(rngIdx, rxIdx) = 10.^((binEdge(medBin) + binEdge(medBin + 1))/2);
            noiseCfg.floorPowFullCh(rngIdx, rxIdx) = 10.^((binEdge(maxBin) + binEdge(maxBin + 1))/2);
        end
    end

    %------------------------------------------------------------
    % 4) 8 通道噪底绘图: 每通道 4 条曲线(本通道 medBin/maxBin + rxPowSum medBin/maxBin)
    %    各 RX 分格纵轴统一(全局范围), 便于一眼对比哪个 RX 不好
    %------------------------------------------------------------
    screenSize = get(0, 'ScreenSize');
    figW = screenSize(3) * 0.98;  figH = screenSize(4) * 0.80;
    figPos = [(screenSize(3)-figW)/2, (screenSize(4)-figH)/2, figW, figH];

    fChDb      = 10*log10(noiseCfg.floorPowCh);       % 本通道 medBin (dB)
    fFullChDb  = 10*log10(noiseCfg.floorPowFullCh);   % 本通道 maxBin (dB)
    fAvgDb     = 10*log10(noiseCfg.floorPow);         % rxPowSum medBin (dB)
    fFullAvgDb = 10*log10(noiseCfg.floorPowFull);     % rxPowSum maxBin (dB)
    all4 = [fChDb(:); fFullChDb(:); fAvgDb(:); fFullAvgDb(:)];
    yl4  = [min(all4), max(all4)];  yl4 = yl4 + 0.06*diff(yl4)*[-1, 1];  % 图6 统一纵轴(留 6% 边距)
    yl5  = [min(fChDb(:)), max(fChDb(:))];  yl5 = yl5 + 0.06*diff(yl5)*[-1, 1];  % 图7 同

    figN = figure('Position', figPos, 'Name', 'fig16_NoiseFloor');
    tiledlayout(2, ceil(nRx/2), 'TileSpacing', 'compact', 'Padding', 'compact');
    binAxis0 = 0 : size(fChDb, 1) - 1;              % 0 基横轴(与 RD 图/gate_R 约定一致)
    for rxIdx = 1 : nRx
        nexttile;
        plot(binAxis0, fChDb(:, rxIdx), 'b-'); hold on;
        plot(binAxis0, fFullChDb(:, rxIdx), 'r-');
        plot(binAxis0, fAvgDb, 'k--');
        plot(binAxis0, fFullAvgDb, 'm--');
        ylim(yl4);                          % 各 RX 统一显示范围
        title(sprintf('RX %d', rxIdx));
        xlabel('距离 bin');
        ylabel('噪底 (dB)');
        legend({'floorPowCh 本通道 medBin', 'floorPowFullCh 本通道 maxBin', ...
                'floorPow rxPowSum medBin', 'floorPowFull rxPowSum maxBin'}, ...
               'Location', 'best', 'FontSize', 6);
        grid on;
    end
    sgtitle('8 通道逐距离门噪底（蓝/红=本通道 medBin/maxBin, 黑/品红虚=rxPowSum medBin/maxBin）');
    noiseCfg.figNoise = figN;

    %------------------------------------------------------------
    % 5) 展示用图: 仅各 RX 的 medBin 噪底(每个 RX 单独一格, 纵轴统一)
    %------------------------------------------------------------
    figN2 = figure('Position', figPos, 'Name', 'fig17_NoiseFloor_med');
    tiledlayout(2, ceil(nRx/2), 'TileSpacing', 'compact', 'Padding', 'compact');
    for rxIdx = 1 : nRx
        nexttile;
        plot(binAxis0, fChDb(:, rxIdx), 'b-', 'LineWidth', 1);
        ylim(yl5);                          % 各 RX 统一显示范围
        title(sprintf('RX %d', rxIdx));
        xlabel('距离 bin');
        ylabel('噪底 (dB)');
        grid on;
    end
    sgtitle('各 RX 通道逐距离门噪底（medBin）');
    noiseCfg.figNoise2 = figN2;
end
