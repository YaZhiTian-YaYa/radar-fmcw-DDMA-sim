function alpha = interpolation_fun_pipe(input)
%INTERPOLATION_FUN_PIPE 亚bin偏移估计(三点加权差分)
%   输入: input -- 5点功率序列, 顺序为 [前2bin, 前1bin, 当前bin, 后1bin, 后2bin]
%   输出: alpha -- 当前bin的亚bin偏移, 限幅到[-0.75, 0.75]
%   注意: 公式仅使用前3个点(前2bin/前1bin/当前bin), 偏移恒>=当前bin一侧
    if(1)
        alpha = (-input(1) + input(3)) / (input(1) + input(2) + input(3));
    end

    % 偏移限幅, 防止插值越过一个bin
    if alpha < -0.75
        alpha = -0.75;
    elseif alpha > 0.75
        alpha = 0.75;
    end
end
