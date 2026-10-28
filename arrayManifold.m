function A = arrayManifold(antenna_t, antenna_r)
%ARRAYMANIFOLD 由 TX/RX 天线相位中心求虚拟阵列流型（虚拟阵元坐标矩阵）
%   A = arrayManifold(antenna_t, antenna_r)
%   输入:
%     antenna_t  [3 × nTX] TX 天线相位中心坐标 [x;y;z] (m)
%     antenna_r  [3 × nRX] RX 天线相位中心坐标 [x;y;z] (m)
%   输出:
%     A  [3 × nTX*nRX] 虚拟阵元坐标(TX+RX 相位中心矢量叠加)。
%        列顺序【TX 内循环 / RX 外循环】:
%        TX1RX1, TX2RX1, ..., TX8RX1, TX1RX2, TX2RX2, ..., TX8RX8
%        (与 rf_sim_core 内部虚拟阵列、后续 DOA 模块的通道顺序一致)
    nTX = size(antenna_t, 2);
    nRX = size(antenna_r, 2);
    A = zeros(3, nTX*nRX);
    idx = 1;
    for j = 1:nRX            % RX 外循环
        for i = 1:nTX        % TX 内循环
            A(:, idx) = antenna_t(:, i) + antenna_r(:, j);
            idx = idx + 1;
        end
    end
end
