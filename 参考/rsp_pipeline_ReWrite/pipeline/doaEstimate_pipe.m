function [tgtOut, detCfg] = doaEstimate_pipe(detCfg, tgtOut)
%DOAESTIMATE_PIPE 方位DBF + 子阵俯仰DBF测角, 结果写回目标列表
%   输入: detCfg.doaChanIn -- 各目标64通道×8chirp复数快照(每目标一行)
%         tgtOut          -- 检测输出的目标列表
%   输出: tgtOut 中补充 azimuth/elevation/doaMag/azimuthList/sideLobe
%   依赖: AntCaliData.aziSteerVec/eleSteerVec -- 方位/俯仰导向矢量(天线标定)
%         calValidChannel_pipe -- 子阵->通道映射(外部函数)
    load('./AdcAndParam.mat','AntCaliData');
    aziStartIdx = 1;                              % 方位导向矢量起始索引(对应-75°)
    subArrTotal = 12;                             % 子阵总数
    eleGrpSteerVec = zeros(subArrTotal, 51);      % 各子阵的俯仰组导向矢量 [12 × 51俯仰格]
    chanTotal = 64;
    aziGrid = -75:75; eleGrid = -25:25;           % 方位/俯仰搜索网格(1°步进)
    subArrUsed = 12;                               % 俯仰测角实际使用的子阵数
    % ---------- 预计算各子阵的俯仰组导向矢量 ----------
    % 取该子阵所有通道的俯仰导向矢量求和后归一化, 作为子阵级俯仰导向矢量
    for subIdx = 1 : subArrTotal
        chanSel = calValidChannel_pipe(subIdx, chanTotal);
        grpSum = sum(AntCaliData.eleSteerVec(chanSel,:), 1);
        eleGrpSteerVec(subIdx, :) = grpSum ./ abs(grpSum);
    end
    %% ---------- 方位DBF: 对全部目标一次性做波束形成 ----------
    % 取目标快照的有效通道(第一子阵通道选择), 与方位导向矩阵相乘,
    % 得到每个目标在[-75°,75°]网格上的DBF谱(dB)
    chanSel = calValidChannel_pipe(1, 64);
    aziDbfOut = detCfg.doaChanIn(:, chanSel) * 1 / length(chanSel) * ...
                    AntCaliData.aziSteerVec(chanSel, aziStartIdx : aziStartIdx + length(aziGrid) - 1);
    doaRes.aziDbfSpec = db(abs(aziDbfOut));       % [目标数 × 151方位格]
    doaRes.eleDbfSpec = zeros(size(tgtOut, 1), length(eleGrid));
    doaRes.item = struct('aziPos', [], 'elePos', [], 'aziSubOff', [], 'eleSubOff', [], 'aziPow', [], 'sideLobeLv', []);
    doaRes.item = repmat(doaRes.item, size(tgtOut, 1), 1);
    % ---------- 逐目标在方位谱上找峰(前4个), 并记录亚格偏移与副瓣电平 ----------
    for tgtIdx = 1 : size(tgtOut, 1)
        [pkPow, pkPos, pkSubOff] = searchPeaksPwr_pipe(doaRes.aziDbfSpec(tgtIdx, :).', 4);
        doaRes.item(tgtIdx).aziPos = single(pkPos);         % 前4个方位峰的格索引
        doaRes.item(tgtIdx).aziSubOff = single(pkSubOff);   % 各峰的亚格偏移
        doaRes.item(tgtIdx).aziPow = single(pkPow);         % 各峰功率
        doaRes.item(tgtIdx).sideLobeLv = single(max(pkPow) - pkPow);   % 相对主峰的副瓣电平
    end
    % ---------- 逐目标俯仰测角: 先子阵合成, 再子阵级俯仰DBF ----------
    for tgtIdx = 1 : size(tgtOut, 1)

        subArrOut = zeros(1, subArrUsed);
        % 每个子阵: 用目标主方位峰处的导向矢量对本子阵通道做匹配滤波(子阵合成),
        % 把相位参考统一到目标方位上
        for subIdx = 1 : subArrUsed
            steerPos = doaRes.item(tgtIdx).aziPos(1) + aziStartIdx - 1;    % 主峰对应的导向矢量索引
            chanSel = calValidChannel_pipe(subIdx, 64);
            steerVec = AntCaliData.aziSteerVec(chanSel, steerPos).';
            subArrOut(subIdx) = steerVec * detCfg.doaChanIn(tgtIdx, chanSel).';
        end

        % 子阵输出 × 子阵俯仰导向矢量 -> [-25°,25°]俯仰谱, 取第1个峰(幅度谱findpeaks)
        eleDbfBuf = subArrOut * eleGrpSteerVec(1:subArrUsed,:);
        pkPow = zeros(4, 1);
        pkPos = zeros(4, 1);
        pkSubOff = zeros(size(eleDbfBuf, 1), 1);

        elePow = db(abs(eleDbfBuf));
        [pkPow, pkPos, pkSubOff] = searchPeaksPwr_pipe(elePow.', 1);

        doaRes.item(tgtIdx).elePos = single(pkPos);         % 俯仰峰格索引
        doaRes.item(tgtIdx).eleSubOff = single(pkSubOff);   % 俯仰亚格偏移
    end
    % 方位谱逐目标归一化(减去各自最大值), 便于对比
    doaRes.aziDbfSpec = single(doaRes.aziDbfSpec - repmat(max(doaRes.aziDbfSpec, [], 2), 1, size(doaRes.aziDbfSpec, 2)));

    % ---------- 结果写回目标列表 ----------
    for tgtIdx = 1 : size(tgtOut, 1)
        % 方位/俯仰角 = 网格角度 + 亚格偏移(一阶精化)
        tgtOut(tgtIdx).azimuth = aziGrid(doaRes.item(tgtIdx).aziPos(1)) + doaRes.item(tgtIdx).aziSubOff(1);
        tgtOut(tgtIdx).elevation = eleGrid(doaRes.item(tgtIdx).elePos(1)) + doaRes.item(tgtIdx).eleSubOff(1);

        tgtOut(tgtIdx).doaMag = doaRes.item(tgtIdx).aziPow(1);              % 主方位峰功率
        tgtOut(tgtIdx).azimuthList = aziGrid(doaRes.item(tgtIdx).aziPos).'; % 前4个候选方位列表
        tgtOut(tgtIdx).sideLobe = doaRes.item(tgtIdx).sideLobeLv.';         % 对应副瓣电平
    end

end
