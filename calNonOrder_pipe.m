function [nonCohImage, nonCohOrder] = calNonOrder_pipe(imageData, ddmList, ddmNum)
%CALNONORDER_PIPE DDM-MIMO速度解模糊: 逐阶对齐跨组平均, 选积累功率最大的阶
%   输入: imageData -- [rngUse距离 × subBins(=16)子带 × grpTotal(=16)组] 功率图
%         ddmList   -- 8个TX的DDM多普勒频移(单位: 组)
%         ddmNum    -- 候选阶数(=组数16, 即速度模糊周期)
%   输出: nonCohImage -- 各单元16个候选阶中的最大积累功率
%         nonCohOrder -- 取得最大值的阶数(对应目标无模糊速度所在分组)

    nonCohImage = zeros(size(imageData));
    % 遍历16种候选模糊阶 idx0:
    % 若目标真实速度对应第idx0阶, 则其8个TX的回波分别落在
    % mod(ddmList + idx0 - 1, 16) + 1 所指示的8个组内;
    % 按该假设把8个组的数据在子带维对齐后平均, 即得该阶下的非相参积累。
    % 【子带边界劈分吸收(2026-11-04)】谱线分数偏移跨 16bin 组边界时能量劈进相邻
    % 两个 bin → 劈分泄漏落邻子带, 低 ncSNR(≈4dB, 300m 级)目标邻带泄漏得分反超
    % 真带 → 选错阶 → 检出位置偏 141bin(实测 12dBsm@300m: 真位 174.6/检出 35)。
    % 修: 每组先取子带 ±1 三列最大(组界回绕)再跨组平均 —— 真带收回劈分能量, 各
    % 候选阶噪声同步抬高(+1.2dB 量级), 对比度不劣化; 高 ncSNR 场景结论不变
    for idx0 = 1 : ddmNum
        validIdx0 = mod(ddmList + idx0 - 1, ddmNum) + 1;   % 该阶下8个TX各自所在的组号
        blk = imageData(:, :, validIdx0);                   % [rng × subBins × 8TX]
        blk3 = max(blk, circshift(blk, [0, 1, 0]));         % 子带 +1 邻域(回绕)
        blk3 = max(blk3, circshift(blk, [0, -1, 0]));       % 子带 −1 邻域
        nonCohImage(:, :, idx0) = mean(blk3, 3);
    end
    % 跨16阶取最大: 目标真实阶对齐后能量叠加最强, 最大值对应的阶即解模糊结果
    [nonCohImage, nonCohOrder] = max(nonCohImage, [], 3);
    sumImage = mean(imageData,3);          % 不对齐直接平均(参考量)

    nonCohRatio = nonCohImage./sumImage;   % 对齐积累相对普通平均的提升比(当前未输出)
end
