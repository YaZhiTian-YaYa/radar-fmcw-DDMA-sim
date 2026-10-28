function [AntCaliData, rowLut] = genAntCaliData(A_manifold, FFT_d0, lambda, cCal, aziRange, eleRange)
%GENANTCALIDATA 生成虚拟天线标定数据(理论导向矢量 × 0°校准), 结构体名与参考一致
%   [AntCaliData, rowLut] = genAntCaliData(A_manifold, FFT_d0, lambda, cCal, aziRange, eleRange)
%   替代参考 AdcAndParam.mat 中的实测 AntCaliData:
%     .aziSteerVec      [64 × nAzi]  方位导向矢量(每通道 × 方位格)
%     .eleSteerVec      [64 × nEle]  俯仰导向矢量
%     .aziSteerVecRange [−75 75 1]   方位格范围(与参考同名同格式)
%     .eleSteerVecRange [−25 25 1]   俯仰格范围
%   理论导向矢量: 通道 c 位于 (x_c, z_c), 对角度 θ 的响应相位 = −k·(坐标)·sinθ;
%   再乘 0°校准表 cCal(共轭已含在 cCal 内) —— 等效把通道相位失配/码常量校准进表内,
%   与参考"实测导向矢量天然含标定"的用法一致。
%   rowLut: [行数 × 64] 通道按 z 坐标(俯仰行)分组查找表, 供 DOA 子阵划分
%   (参考阵列: 方位最长行 30 元/俯仰 12 行; 我们的阵列: 最长行 36 元/俯仰 9 行,
%    全部由 A_manifold 读数据得出, 不人为输入)
    xC = A_manifold(1, :);                       % 通道 x (m) —— 导向矢量必须用物理米坐标
    zC = A_manifold(3, :);                       % 通道 z (m) → 俯仰行划分依据
    k  = 2*pi/lambda;

    % ---- 俯仰行划分(读数据): 按 z 分组 ----
    [zU, ~, iz] = unique(zC);
    nRow = numel(zU);
    rowLut = false(nRow, numel(xC));
    for rr = 1:nRow
        rowLut(rr, :) = (iz == rr).';
    end

    % ---- 角度格(与参考同格式) ----
    aziGrid = aziRange(1) : aziRange(3) : aziRange(2);
    eleGrid = eleRange(1) : eleRange(3) : eleRange(2);

    % ---- 理论导向矢量(方位含 0°校准; 俯仰纯理论) ----
    % 符号约定: 本仿真 beat 相位 = +2π·f_c·τ → 通道相位含 −k·x·sinθ,
    % 导向矢量取共轭匹配 +k·x·sinθ (DBF 峰在真实角度处)
    % 校准施加次数: 方位 DBF = 快照×aziSteerVec → cCal 施加 1 次(校准);
    %   俯仰 DBF = 子阵输出(已含 aziSteerVec 的 1 次)×eleGrpSteerVec —— 若
    %   eleSteerVec 也含 cCal 则施加 2 次(cCal²) → 俯仰系统性偏 ~+0.5°,
    %   故 eleSteerVec 用纯理论(俯仰维的校准由方位匹配滤波已带的那次完成)
    cCalCol = cCal(:);
    aziSteerVec = exp(+1i * k * (xC(:) * sind(aziGrid)))  .* cCalCol;   % [64 × nAzi] 含校准
    eleSteerVec = exp(+1i * k * (zC(:) * sind(eleGrid)));               % [64 × nEle] 纯理论

    AntCaliData = struct('aziSteerVec', aziSteerVec, 'eleSteerVec', eleSteerVec, ...
                         'aziSteerVecRange', aziRange, 'eleSteerVecRange', eleRange);
end
