function data = sva_fun_pipe(data, range_svaEn, doppler_svaEn)
%SVA_FUN_PIPE 自适应旁瓣抑制/泄漏抵消（SVA）
%   输入:
%     data          : 3D 频域数据（通常维度可理解为 [range(距离) , channel/virtual array , doppler(多普勒)])
%     range_svaEn   : 距离维 SVA 开关（1=启用，0=不启用）
%     doppler_svaEn : 多普勒维 SVA 开关（1=启用，0=不启用）
%   输出:
%     data          : 经过 SVA 处理后的数据（与输入同尺寸）

    % 若距离维和多普勒维都不启用，直接返回原数据
    if(range_svaEn == 0 && doppler_svaEn == 0), return; end

    % ---------------------------
    % 1) 距离维方向的 SVA
    % ---------------------------
    if(range_svaEn == 1)
        % 为了让“边界 bin”也能用到三点邻域（prev / curr / next），
        % 这里在距离维（第1维）做首尾补零：
        %   data(end,:,:) 作为“最后一个距离门”
        %   data(1,:,:)     作为“第一个距离门”
        % cat(1, ...) 表示沿第1维拼接
        data = cat(1, zeros(size(data(end, :, :))), data, zeros(size(data(1, :, :))));

        % 取三点窗口：
        %   data1 对应“左邻居”（prev）
        %   data2 对应“当前 bin”（curr）
        %   data3 对应“右邻居”（next）
        data1 = data(1 : end - 2, :, :);
        data2 = data(2 : end - 1, :, :);
        data3 = data(3 : end - 0, :, :);

        % 计算自适应系数 alpha：
        %   形式：alpha = 2 * real( data2 * conj(data1 + data3) ) / |data1 + data3|^2
        %
        % 直观理解：
        %   - 用“左+右邻居的复值和”去预测/解释当前 bin 的成分
        %   - alpha 越大，说明当前 bin 与左右邻居的相关性越强
        %   - 再根据相对强度把“可由泄漏解释的部分”估出来
        %
        % 其中 +1e-9 是为了防止分母为 0 导致数值爆炸
        alpha = 2 * real(data2 .* conj(data1 + data3)) ./ (abs(data1 + data3) .^ 2 + 1e-9);

        % 将 alpha 限幅到 [0, 1]
        %   - 小于0：置0
        %   - 大于1：置1
        % 另外把边界（第一行和最后一行）的 alpha 强制为 1.0
        % 这是为了避免补零/边界导致的不合理系数
        alpha(alpha < 0) = 0;
        alpha(alpha > 1) = 1;
        alpha(1, :, :) = 1.0;
        alpha(end, :, :) = 1.0;

        % 最终 SVA 更新（“泄漏抵消”核心公式）：
        %   data = data2 - (alpha/2) * (data1 + data3)
        %
        % 含义：
        %   - 把当前 bin(data2)中“能由左右邻居解释/对应的泄漏成分”
        %     估计出来并减掉
        data = data2 - alpha / 2 .* (data1 + data3);
    end

    % ---------------------------
    % 2) 多普勒维方向的 SVA
    % ---------------------------
    if(doppler_svaEn == 1)
        % 为了处理多普勒维（第3维）边界，同样做首尾循环拼接：
        %   data(:,:,end) 作为“尾端邻居”补到前面
        %   data(:,:,1)     作为“头端邻居”补到后面
        % cat(3, ...) 表示沿第3维拼接
        data = cat(3, data(:, :, end), data, data(:, :, 1));

        % 取三点窗口（沿第3维）：
        data1 = data(:, :, 1 : end - 2);  % prev
        data2 = data(:, :, 2 : end - 1);  % curr
        data3 = data(:, :, 3 : end - 0);  % next

        % 同样计算自适应系数 alpha（公式与距离维相同，只是沿第3维取窗口）
        alpha = 2 * real(data2 .* conj(data1 + data3)) ./ (abs(data1 + data3) .^ 2 + 1e-9);

        % alpha 限幅到 [0, 1]
        alpha(alpha < 0) = 0;
        alpha(alpha > 1) = 1;

        % SVA 更新：从当前 bin 中减去估计的泄漏成分
        data = data2 - alpha / 2 .* (data1 + data3);
    end
end