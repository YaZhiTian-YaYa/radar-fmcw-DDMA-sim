% 泄漏行(PN裙边)探针: 行1(bin0直流泄漏)+行3 近距行 多普勒剖面 vs 空行
diary('C:\Users\NJ\AppData\Local\Temp\probe95\leakRow.log'); clc; diary on;
addpath('F:\Work3\CH_NEW_95_射频_DDMA');
saveFigDir = 'C:\Users\NJ\AppData\Local\Temp\probe95\cohfigs';
if ~exist(saveFigDir,'dir'), mkdir(saveFigDir); end
batchState.idx = 1;  batchState.okRun = true;
run('F:\Work3\CH_NEW_95_射频_DDMA\main_sim.m');
rE = 350; rxI = 7;
rowE = squeeze(20*log10(abs(fftCfg.data_sva(rE, rxI, :)) + eps));
for r0 = [1 3 5]
    rowT = squeeze(20*log10(abs(fftCfg.data_sva(r0, rxI, :)) + eps));
    [~, pk] = maxk(rowT, 3);
    m = true(1,256);
    for k=1:3, m(mod((pk(k)-4:pk(k)+3),256)+1) = false; end
    fprintf('LK 行%d: 行峰 %7.2f | 排峰±4中位 %7.2f | 空行 %7.2f | 差 %+.2f\n', ...
        r0, max(rowT), median(rowT(m)), median(rowE), median(rowT(m))-median(rowE));
end
diary off;
