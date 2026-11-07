% DDMA 泄露取证 v2: 抬升能量落点按"到最近码线距离"分桶 + 码线自身峰变化
diary('C:\Users\NJ\AppData\Local\Temp\probe95\diagLeak2.log'); clc; diary on;
addpath('F:\Work3\CH_NEW_95_射频_DDMA');
T = @(gR,gD,az,el,rcs,ph) struct('gate_R',gR,'gate_D',gD,'azi',az,'ele',el, ...
    'rcsDb',rcs,'rcsPhDeg',ph,'vTan',0,'isCorner',0);
C = struct('gate_R',5,'gate_D',0,'azi',0,'ele',0,'rcsDb',-5,'rcsPhDeg',0,'vTan',0,'isCorner',1);
targetsOverride = [C, T(285,0,0,0,10,0)];
seedOverride = 4;
saveFigDir = 'C:\Users\NJ\AppData\Local\Temp\probe95\cohfigs';
if ~exist(saveFigDir,'dir'), mkdir(saveFigDir); end
batchState.idx = 1;  batchState.okRun = true;
run('F:\Work3\CH_NEW_95_射频_DDMA\main_sim.m');

rT = squeeze(10*log10(fftCfg.rxPowSum(286, :) + eps));   % 目标行
rE = squeeze(10*log10(fftCfg.rxPowSum(350, :) + eps));   % 空行
linesK = squeeze(double(accumCfg.chirpPosMap(286, 1, :)));  % [8×1]
fprintf('LEAK2 空行中位 %7.2f | 目标行全行中位 %7.2f (抬升 %+.2f)\n', ...
    median(rE), median(rT), median(rT) - median(rE));
% 到最近码线的距离
dMin = zeros(1, 256);
for b = 1:256
    dMin(b) = min(abs(mod(b - linesK + 128, 256) - 128));
end
edges = [0 1 2 4 8 16 32 128];
for q = 1:numel(edges) - 1
    mm = dMin >= edges(q) & dMin < edges(q + 1);
    if any(mm)
        fprintf('LEAK2 距线 %3d-%3d: n=%3d 中位 %7.2f 抬升 %+.2f dB\n', ...
            edges(q), edges(q + 1) - 1, sum(mm), median(rT(mm)), median(rT(mm)) - median(rE));
    end
end
% 码线峰高(线损检验: OU 与白噪的 E|e^{jφ}| 相同 → 峰应基本不变)
for k = 1:8
    seg = rT(mod(linesK(k) - 2 : linesK(k) + 2, 256) + 1);
    fprintf('LEAK2 线%d bin%3d 峰 %7.2f\n', k, linesK(k), max(seg));
end
diary off;
