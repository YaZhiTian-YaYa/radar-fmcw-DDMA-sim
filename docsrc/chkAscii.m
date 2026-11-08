% ASCII-only state check (bypass CJK console corruption)
diary('C:\Users\NJ\AppData\Local\Temp\probe95\chkAscii.log'); clc; diary on;
addpath('F:\Work3\CH_NEW_95_射频_DDMA');
saveFigDir = 'C:\Users\NJ\AppData\Local\Temp\probe95\cohfigs';
batchState.idx = 1;  batchState.okRun = true;
run('F:\Work3\CH_NEW_95_射频_DDMA\main_sim.m');
rxI = 7;
rowT = squeeze(20*log10(abs(fftCfg.data_sva(144, rxI, :)) + eps));
rowE = squeeze(20*log10(abs(fftCfg.data_sva(300, rxI, :)) + eps));
rowT = rowT(:).';  rowE = rowE(:).';
fMed = median(rowE);
[~, pk] = maxk(rowT, 8);
mL = false(1,256);
for k = 1:8, mL(mod((pk(k)-1-2):(pk(k)-1+2),256)+1) = true; end
g16 = unique(mod(pk-1,16));
mE = false(1,256);
for s = (0:15)*16 + g16(1)
    if ~ismember(mod(s,256), pk-1), mE(mod(mod(s,256)-2:mod(s,256)+2,256)+1) = true; end
end
cz = ~mL & ~mE;
v = rowT(cz);
fprintf('ASC fMed=%.2f nCZ=%d min=%.2f q25=%.2f med=%.2f q75=%.2f max=%.2f nOver6=%d\n', ...
    fMed, numel(v), min(v), prctile(v,25), median(v), prctile(v,75), max(v), sum(v > fMed+6));
fprintf('ASC cr=%g lines=%s\n', prm.cr, mat2str(sort(pk-1).'));
diary off;
