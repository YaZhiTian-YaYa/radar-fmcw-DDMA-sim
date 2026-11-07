addpath('F:\Work3\CH_NEW_95_射频_DDMA');
saveFigDir='C:\Users\NJ\AppData\Local\Temp\probe95\cohfigs'; batchState.idx = 1; batchState.okRun = true;
run('F:\Work3\CH_NEW_95_射频_DDMA\main_sim.m');
r0 = 429; rE = 350; rxI = 7;
rowT = squeeze(20*log10(abs(fftCfg.data_sva(r0, rxI, :)) + eps));
rowE = squeeze(20*log10(abs(fftCfg.data_sva(rE, rxI, :)) + eps));
[~, pkI] = maxk(rowT, 8);
dMin = zeros(1,256);
for b=1:256, dMin(b)=min(abs(mod(b-pkI+128,256)-128)); end
fprintf('W2 目标行峰 %7.2f | 空行 %7.2f\n', max(rowT), median(rowE));
for q = [1 2 3 4 8 16 32 128]
    mm = dMin>=q & dMin<q*2;
    if any(mm), fprintf('W2 距线 %3d-%3d: n=%3d 中位 %7.2f 差 %+.2f\n', q, q*2-1, sum(mm), median(rowT(mm)), median(rowT(mm))-median(rowE)); end
end
mAll = true(1,256); mAll(mod((pkI-4:pkI+3),256)+1)=false;
fprintf('W2 排线±4后中位 %7.2f 差 %+.2f\n', median(rowT(mAll)), median(rowT(mAll))-median(rowE));
