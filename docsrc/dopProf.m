% 目标行多普勒剖面分桶(距码线远近) —— 与用户实测"两侧都低, 远处才抬肩"对照
addpath('F:\Work3\CH_NEW_95_射频_DDMA');
saveFigDir='C:\Users\NJ\AppData\Local\Temp\probe95\cohfigs'; batchState.idx = 1; batchState.okRun = true;
run('F:\Work3\CH_NEW_95_射频_DDMA\main_sim.m');
r0 = 429; rE = 350; rxI = 7;
rowT = squeeze(20*log10(abs(fftCfg.data_sva(r0, rxI, :)) + eps));
rowE = squeeze(20*log10(abs(fftCfg.data_sva(rE, rxI, :)) + eps));
[~, pkI] = maxk(rowT, 8);
dMin = zeros(1,256);
for b=1:256, dMin(b)=min(abs(mod(b-pkI+128,256)-128)); end
edges = [1 2 4 8 16 32 64 128];
for q=1:numel(edges)-1
    mm = dMin>=edges(q) & dMin<edges(q+1);
    if any(mm), fprintf('DOP 距线 %3d-%3d: 中位 %7.2f 差 %+.2f\n', edges(q), edges(q+1)-1, median(rowT(mm)), median(rowT(mm))-median(rowE)); end
end
