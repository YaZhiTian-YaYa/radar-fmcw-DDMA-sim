% 逐TX窗相位假设检验: 目标多普勒分数偏移 δ 使码线逐TX相位≠公共相位 → 12% 残差
% 修法候选: 原子逐TX乘 hann 核相位 W(chirpPos_k - δ)。检验 J 修正后是否塌到噪声底。
diary('C:\Users\NJ\AppData\Local\Temp\probe95\diagWph.log'); clc; diary on;
addpath('F:\Work3\CH_NEW_95_射频_DDMA');
saveFigDir = 'C:\Users\NJ\AppData\Local\Temp\probe95\cohfigs';
if ~exist(saveFigDir, 'dir'), mkdir(saveFigDir); end
T = @(gR,gD,az,el,rcs,ph) struct('gate_R',gR,'gate_D',gD,'azi',az,'ele',el, ...
    'rcsDb',rcs,'rcsPhDeg',ph,'vTan',0,'isCorner',0);
C = struct('gate_R',5,'gate_D',0,'azi',0,'ele',0,'rcsDb',-5,'rcsPhDeg',0,'vTan',0,'isCorner',1);
targetsOverride = [C, T(285,0,-0.5,0,10,0), T(285,0,0.5,0,2,90), T(257,40,-8,0,5,0)];  % P7
seedOverride = 4;
batchState.idx = 1;  batchState.okRun = true;
run('F:\Work3\CH_NEW_95_射频_DDMA\main_sim.m');
aziGrid = AntCaliData.aziSteerVecRange(1):AntCaliData.aziSteerVecRange(3):AntCaliData.aziSteerVecRange(2);
[~, aziRowBest] = max(sum(rowLut, 2));
chanAzi = find(rowLut(aziRowBest, :));
A_dml = conj(AntCaliData.aziSteerVec(chanAzi, :));
fitFun = @(th) steerFitAt(A_dml, aziGrid, th);
dAz = abs([tgtOut.azimuth] - (-0.5));
nRow = size(detCfg.doaChanIn, 1);
cand = find(dAz < 1.5 & (1:numel(tgtOut)) <= nRow);
[~, im] = max([tgtOut(cand).doaMag]);  t = cand(im);
x = double(detCfg.doaChanIn(t, chanAzi)).';
C0 = real(x' * x);
% ---- 通道的 TX 归属与该 cell 的码线 bin ----
txOf = ones(64, 1);  for c = 1:64, txOf(c) = floor((c - 1) / 8) + 1; end   % TX内循环/RX外循环
txCh = txOf(chanAzi);                                       % 36 通道各自的 TX
rIdx0 = double(tgtOut(t).rangIdx);  sb0 = mod(double(tgtOut(t).dopplerIdx) - 1, dims.subBins) + 1;
binsK = squeeze(double(accumCfg.chirpPosMap(rIdx0, sb0, :)));  % 8 根码线被抓 bin [8×1]
% 目标多普勒分数位置: 检测的 dopplerIdx 是解模糊后 bin; 各码线的目标真值位置 = binsK - alpha
alD = double(tgtOut(t).doppler_alpha);                      % 亚bin偏移(约定: alpha>0 峰偏前侧低bin)
fprintf('WPH 点#%d: dopIdx %d alpha %.3f, 码线 bins %s\n', t, tgtOut(t).dopplerIdx, alD, mat2str(binsK.'));
% ---- 检验 A: 原子按 TX 乘 hann 核(幅度×相位) @目标线位置 binsK-al ----
% hann 核(DFT of 0.5-0.5cos): 对频率 f(bin 单位, 相对 bin0 中心), N 点 DFT 在 bin k 的值
Nfft = dims.dopBinNum;
hannKer = @(f, k) 0.5*dirichlet(f-k, Nfft) - 0.25*dirichlet(f-k-1, Nfft) - 0.25*dirichlet(f-k+1, Nfft);
% dirichlet 近似(sinc): D(v)=sin(pi v)/(Nfft*sin(pi v/Nfft))*e^{j pi v (1-1/Nfft)}
dirichlet = @(v, N) (sin(pi*v) ./ (N * sin(pi*v/N + eps)) .* exp(1i*pi*v*(1 - 1/N)));
hannKer = @(f, k) 0.5*dirichlet(f-k, Nfft) - 0.25*dirichlet(f-k-1, Nfft) - 0.25*dirichlet(f-k+1, Nfft);
fTgt = -alD;                                                % 目标频率相对主 bin 的偏移(bin)
wTx = hannKer(fTgt, 0) * ones(1, 8);                        % 各码线相对其 bin 中心同偏移 → 同值?
% 注意: 各码线 chirpPos 不同但目标线位置 = chirpPos + fTgt → 每根线的核值都 = hannKer(fTgt, 0) → 公共因子!
% → 若两目标同 cell 同 δ, 逐TX窗相位是公共因子, 会被 LS 幅度吸收 —— 假设不成立?
fprintf('WPH 各码线核值(应公共): %s\n', mat2str(abs([hannKer(fTgt,0) hannKer(fTgt,0)]), 3));
% ---- 检验 B: 残差到底长什么样 —— 对 x 做逐TX相位剥离后再看能量占比 ----
% 对每 TX 的通道组求 x 相对 a(-0.5) 拟合的相位残差(公共幅度扣除后)
At = [fitFun(-0.5), fitFun(0.5)];
aT = inv2x2(At'*At) * (At'*x);
r0 = x - At * aT;
fprintf('WPH 双角(±0.5)残差占比 %.3f; 单角(-0.5)残差占比 %.3f\n', ...
    real(r0'*r0)/C0, 1 - abs(fitFun(-0.5)'*x)^2/real(fitFun(-0.5)'*fitFun(-0.5))/C0);
% 逐 TX 看残差能量分布
eTx = zeros(1, 8);
for k = 1:8
    m = (txCh == k);
    eTx(k) = real(r0(m)' * r0(m));
end
fprintf('WPH 残差逐TX能量: %s (均匀=噪声/结构, 集中=特定TX)\n', mat2str(round(eTx/sum(eTx)*100)).');
% 逐 TX 的单原子(-0.5)投影相位(相对全通道LS相位): 若逐TX相位不同 → 窗相位结构实锤
a1v = fitFun(-0.5);  al1 = (a1v' * x) / real(a1v' * a1v);
phTx = zeros(1, 8);
for k = 1:8
    m = (txCh == k);
    pk = angle((a1v(m)' * x(m)) / (al1 * real(a1v(m)' * a1v(m)) * ones(sum(m),1)));
    phTx(k) = rad2deg(mean(angle(x(m) .* conj(a1v(m)) ) - angle(al1)));
end
fprintf('WPH 逐TX相对相位(deg): %s\n', mat2str(round(phTx,1)).');
% ---- 决定性核对: 真角度幅度恢复 + 邻空cell噪声底 ----
At2 = [fitFun(-0.5), fitFun(0.5)];
aLs = inv2x2(At2'*At2) * (At2' * x);
fprintf('WPH 真角度LS幅度: [%6.2f %6.2f] dB (真值差8dB), 残差占比 %.3f\n', ...
    20*log10(abs(aLs(1))+eps), 20*log10(abs(aLs(2))+eps), ...
    real((x-At2*aLs)'*(x-At2*aLs))/C0);
sbE = sb0 + 2;  if sbE > dims.subBins, sbE = sb0 - 2; end
xE = double(accumCfg.chanSnapshot(chanAzi, rIdx0, sbE));
fprintf('WPH 邻空cell能量占比(噪声参考) %.3f\n', real(xE'*xE)/C0);
diary off;
