function cohStat = cohStatCalc(nChan, nFft, nRXavg)
%COHSTATCALC 相参口径统计项 cohStat 的蒙特卡洛理论值(角度维取峰 → 噪底参考抬升)
%   cohStat = 10lg(nChan·F_nc / F_coh):
%     F_coh = 相参噪底中位: nChan 通道复高斯 → nFft 点方位 FFT(补零, 与 jointAccum
%             的角度 FFT 同构)逐次取 max|X|², 其全体中位(每角度 bin 噪声均值=nChan·σ²)
%             —— bin 间相关使有效独立数介于 nChan 与 nFft 之间(64/128 时 ≈85,
%             max 中位 ≈ 5.06×单bin均值)
%     F_nc  = 非相参噪底中位: nRXavg 路 |·|² 功率平均的中位(Gamma 型, 与逐距离门
%             噪底 floorPow 同构, ≈0.96σ²)
%   依赖雷达结构(nChan/nFft/nRXavg)—— 换雷达(如 16T16R: nChan=256)自动重算,
%   勿写死常数(2026-10-29 用户要求; 95 当前 64/128/8 → ≈−7.2 dB)。
%   只含 max-of-bin 统计的理论值; 残差(理论略优于仿真, ~2 dB 口径小量)由审计
%   如实显示, 勿凑差为零。固定种子(rng 7), 结果可复现。
    if nargin < 2 || isempty(nFft),   nFft   = 128; end
    if nargin < 3 || isempty(nRXavg), nRXavg = 8;   end
    Ntry = 2e5;  chunk = 2e4;  rng(7);             % 固定种子(纯统计量, 可复现)
    mxAll = zeros(Ntry, 1);
    for q = 1:chunk:Ntry
        n  = min(chunk, Ntry - q + 1);
        Z  = (randn(nChan, n) + 1i*randn(nChan, n)) / sqrt(2);   % nChan 通道, E|z|²=1
        X  = fft([Z; zeros(nFft - nChan, n)]);                   % nFft 点方位 FFT(补零)
        mxAll(q:q+n-1) = max(abs(X).^2, [], 1).';
    end
    F_coh = median(mxAll);                                       % 相参噪底中位(σ²=1 口径)
    nc = mean(abs((randn(nRXavg, Ntry) + 1i*randn(nRXavg, Ntry))/sqrt(2)).^2, 1).';
    F_nc = median(nc);                                           % 非相参噪底中位
    cohStat = 10*log10(nChan * F_nc / F_coh);
    fprintf('cohStatCalc      : 理论 cohStat = %+.2f dB(蒙特卡洛 %d 次: nChan=%d/nFft=%d/nRXavg=%d, max中位 %.2f×单bin均值, 非相参中位 %.3f)\n', ...
            cohStat, Ntry, nChan, nFft, nRXavg, F_coh/nChan, F_nc);
end
