function [accumCfg] = nonCohAccumulate_fun_pipe(fftCfg)
%NONCOHACCUMULATE_FUN_PIPE 非相参积累: DDM分组重排 + 逐阶对齐平均并选阶
%   输出: accumCfg.nonCohMap -- 非相参积累功率图 [512距离 × 32多普勒子带]
%         accumCfg.ddmOrdMap -- 每个距离-多普勒单元选中的DDM阶数(速度解模糊用)

    % 512个多普勒bin按DDM分组结构重排为 [512距离 × 32子带 × 16组]
    % (多普勒bin序号 = 子带序号 + (组号-1)*32)
    imageData = reshape(fftCfg.rxPowSum, [512, ...
                                             32, ...
                                             16]);
    % ddmList = [0,14,12,9,8,5,4,1]: 8个发射通道各自的DDM多普勒频移(单位: 组)
    % 对16种候选模糊阶做"按频移对齐 + 跨组平均", 取积累功率最大的阶
    [nonCohImage, nonCohOrder] = calNonOrder_pipe(imageData, ...
                            [0,14,12,9,8,5,4,1], ...
                            16);
    accumCfg.nonCohMap = single(squeeze(nonCohImage));   % 选阶后的非相参积累图
    accumCfg.ddmOrdMap = single(squeeze(nonCohOrder));   % 逐单元最优DDM阶数

end
