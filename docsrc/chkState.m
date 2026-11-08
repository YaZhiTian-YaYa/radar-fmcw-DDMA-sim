% 状态核验: 热噪参考行的滤波增益 + 目标行切面 + 码线位置 (干净打印)
diary('C:\Users\NJ\AppData\Local\Temp\probe95\chkState.log'); clc; diary on;
addpath('F:\Work3\CH_NEW_95_射频_DDMA');
saveFigDir = 'C:\Users\NJ\AppData\Local\Temp\probe95\cohfigs';
batchState.idx = 1;  batchState.okRun = true;
run('F:\Work3\CH_NEW_95_射频_DDMA\main_sim.m');
rxI = 7;
fl = 10*log10(mean(noiseCfg.floorPowCh, 2) + eps);  fl = fl(:);
fprintf('ST cr=%g | 逐行噪底: 行144 %.2f 行200 %.2f 行300 %.2f 行350 %.2f 行400 %.2f 行500 %.2f\n', ...
    prm.cr, fl(144), fl(200), fl(300), fl(350), fl(400), fl(500));
rowT = squeeze(20*log10(abs(fftCfg.data_sva(144, rxI, :)) + eps));
rowE = squeeze(20*log10(abs(fftCfg.data_sva(300, rxI, :)) + eps));
fMed = median(rowE);
[~, pk] = maxk(rowT, 8);
fprintf('ST 参考行300中位 %.2f | 码线0基:\n', fMed);
disp(sort(pk-1).');
mL = false(1,256);
for k = 1:8, mL(mod((pk(k)-1-2):(pk(k)-1+2),256)+1) = true; end
g16 = unique(mod(pk-1,16));
mE = false(1,256);
for s = (0:15)*16 + g16(1)
    if ~ismember(mod(s,256), pk-1), mE(mod(mod(s,256)-2:mod(s,256)+2,256)+1) = true; end
end
cz = ~mL & ~mE;
fprintf('ST 干净区 n=%d 中位 %.2f (+%.2f vs 行300) | 超+6dB bin=%d\n', ...
    sum(cz), median(rowT(cz)), median(rowT(cz))-fMed, sum(rowT>fMed+6 & cz));
diary off;
