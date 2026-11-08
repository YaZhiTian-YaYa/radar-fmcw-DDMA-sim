% ASCII ring-metric: per-distance-ring medians (Hann decay check) + farCZ
diary('C:\Users\NJ\AppData\Local\Temp\probe95\ringChk.log'); clc; diary on;
addpath('F:\Work3\CH_NEW_95_射频_DDMA');
saveFigDir = 'C:\Users\NJ\AppData\Local\Temp\probe95\cohfigs';
batchState.idx = 1;  batchState.okRun = true;
run('F:\Work3\CH_NEW_95_射频_DDMA\main_sim.m');
rxI = 7;
rowT = squeeze(20*log10(abs(fftCfg.data_sva(144, rxI, :)) + eps)); rowT = rowT(:).';
rowE = squeeze(20*log10(abs(fftCfg.data_sva(300, rxI, :)) + eps)); rowE = rowE(:).';
fMed = median(rowE);
[~, pk] = maxk(rowT, 8);
dL = zeros(1,256);
for b = 1:256, dL(b) = min(abs(mod(b-1 - (pk-1) + 128, 256) - 128)); end
fprintf('RG fMed=%.2f | ', fMed);
for q = [3 5; 6 9; 10 20; 21 45]
    fprintf('d%d-%d:%+.1f ', q(1), q(2), median(rowT(dL>=q(1) & dL<=q(2))) - fMed);
end
fprintf('| far(d>=8) %+.1f\n', median(rowT(dL>=8)) - fMed);
diary off;
