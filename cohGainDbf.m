function gDb = cohGainDbf(zRow, lambda, eleDeg)
%COHGAINDBF 通道相参积累增益(阵列 DBF 口径, 通道幅度归一化)——供 snrAnaCalc 的 cohGainDb 实参
%   gDb = 10lg( |Σ_c e^{jk·z_c·sinθel}|² / N )
%   【归一化口径(用户 2026-10-29: 算法方案里 DBF 权重即归一化)】通道幅度取单位
%   (不含 TX_A/RX_A 锥削) —— 逐天线损耗差的均值损耗已由 chTerm(rfInfo.chLossDb)
%   显式列项, 勿在相参增益里重复计入; 归一化后 θel=0 恒等于公式值 10lg(N)。
%   物理口径: 与仿真一维 x 方位 FFT 同构的"强行 N 元相参"—— 方位角任意(x 相位由
%   角度 FFT 取峰自行对齐, 不进本式), 俯仰角≠0 时 z 行相位(k·z·sinθel)对不齐 →
%   增益自动低于 10lg(N)(此时收益变小属预期, 用户口径)。
%   输入: zRow = A_manifold(3,:).'(虚拟阵 z 坐标列向量, m), lambda 波长(m),
%         eleDeg 目标俯仰角(deg)
    phZ = 2*pi/lambda * zRow(:) * sind(eleDeg);         % 残余 z 相位 (rad)
    gDb = 10*log10(abs(sum(exp(1j*phZ))).^2 / numel(phZ));
end
