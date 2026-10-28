function [accumCfg] = nonCohAccumulate_fun_pipe(fftCfg, ddmList, dims)
%NONCOHACCUMULATE_FUN_PIPE 非相参积累: DDM分组重排 + 逐阶对齐平均并选阶
%   (移植自 参考\rsp_pipeline_ReWrite\pipeline\nonCohAccumulate_fun_pipe.m,
%    ddmList/dims 改为入参 —— 用我们自己的 DDMA 编码与主函数维度)
%   输入: fftCfg.rxPowSum -- RX平均功率图 [距离 × 多普勒]
%         ddmList        -- 8个TX的DDM多普勒频移(单位: 组, = find(encodeSeq)-1)
%         dims           -- 公共维度结构体(主函数唯一定义)
%   输出: accumCfg.nonCohMap -- 非相参积累功率图 [距离 × 32多普勒子带]
%         accumCfg.ddmOrdMap -- 每个距离-多普勒单元选中的DDM阶数(速度解模糊用)

    % 多普勒bin按DDM分组结构重排为 [距离 × 子带 × 组] (bin序号 = 子带 + (组-1)×子带bin数)
    % 注意第一维是【距离数】(读数据取 size, 勿写 dopBinNum —— 原参考流水线两维恰好
    % 都是 512 才碰巧能跑, nChirps 改成 256 等其它长度时会 reshape 报错)
    imageData = reshape(fftCfg.rxPowSum, [size(fftCfg.rxPowSum, 1), ...
                                             dims.subBins, ...
                                             dims.grpTotal]);
    % 对 grpTotal 种候选模糊阶做"按频移对齐 + 跨组平均", 取积累功率最大的阶
    [nonCohImage, nonCohOrder] = calNonOrder_pipe(imageData, ...
                            ddmList, ...
                            dims.grpTotal);
    accumCfg.nonCohMap = single(squeeze(nonCohImage));   % 选阶后的非相参积累图
    accumCfg.ddmOrdMap = single(squeeze(nonCohOrder));   % 逐单元最优DDM阶数

end
