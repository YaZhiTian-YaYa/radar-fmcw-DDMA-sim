function data = sva_fun_pipe(data, range_svaEn, doppler_svaEn)
%SVA_FUN_PIP 空间变异变迹(SVA)副瓣抑制
%   对指定维相邻三个频点 X(k-1), X(k), X(k+1) 做自适应加权:
%       Y(k) = X(k) - alpha/2 * (X(k-1) + X(k+1)),   alpha ∈ [0,1]
%   该式等效于逐bin自适应选择窗系数的余弦类加窗: alpha=0 即不加窗,
%   alpha=1 相当于Hann窗; 系数由信号自身按输出模值最小(副瓣最小)原则确定,
%   因此既抑制副瓣又不展宽主瓣
%   range_svaEn / doppler_svaEn: 距离维 / 多普勒维处理使能
    if(range_svaEn == 0 && doppler_svaEn == 0), return; end

    if(range_svaEn == 1)
        % ---- 距离维: 非周期信号, 上下边界各补一行零 ----
        data = cat(1, zeros(size(data(end, :, :))), data, zeros(size(data(1, :, :))));
        data1 = data(1 : end - 2, :, :);        % 左邻bin X(k-1)
        data2 = data(2 : end - 1, :, :);        % 中心bin X(k)
        data3 = data(3 : end - 0, :, :);        % 右邻bin X(k+1)
        % 最优加权系数: alpha = 2*Re{X(k)*conj(X(k-1)+X(k+1))} / |X(k-1)+X(k+1)|^2
        alpha = 2 * real(data2 .* conj(data1 + data3)) ./ (abs(data1 + data3) .^ 2 + 1e-9);
        alpha(alpha < 0) = 0;  alpha(alpha > 1) = 1; alpha(1, :, :) = 1.0; alpha(end, :, :) = 1.0;  % 限幅到[0,1]; 补零的边界行取满窗
        data = data2 - alpha / 2 .* (data1 + data3);   % 自适应加权输出
    end

    if(doppler_svaEn == 1)
        % ---- 多普勒维: 频谱周期循环, 边界做循环延拓 ----
        data = cat(3, data(:, :, end), data, data(:, :, 1));
        data1 = data(:, :, 1 : end - 2);        % 左邻bin X(k-1)
        data2 = data(:, :, 2 : end - 1);        % 中心bin X(k)
        data3 = data(:, :, 3 : end - 0);        % 右邻bin X(k+1)
        % 同距离维: 求最优加权系数并限幅到[0,1]
        alpha = 2 * real(data2 .* conj(data1 + data3)) ./ (abs(data1 + data3) .^ 2 + 1e-9);
        alpha(alpha < 0) = 0;  alpha(alpha > 1) = 1;
        data = data2 - alpha / 2 .* (data1 + data3);   % 自适应加权输出
    end
end
