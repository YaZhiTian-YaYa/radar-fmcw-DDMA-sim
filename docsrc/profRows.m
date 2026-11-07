addpath('F:\Work3\CH_NEW_95_射频_DDMA');
saveFigDir='C:\Users\NJ\AppData\Local\Temp\probe95\cohfigs'; batchState.idx = 1; batchState.okRun = true;
run('F:\Work3\CH_NEW_95_射频_DDMA\main_sim.m');
rows = (340:510);
med = zeros(numel(rows),1);
for q = 1:numel(rows)
    r = squeeze(20*log10(abs(fftCfg.data_sva(rows(q), 7, :)) + eps));
    [~, pk] = maxk(r, 16);
    m = true(1,256); m(mod((pk-3:pk+3),256)+1) = false;
    med(q) = median(r(m));
end
sel = @(a,b) round(med(rows>=a & rows<=b).',1);
fprintf('PROF 340-360: %s\n', mat2str(sel(340,360)));
fprintf('PROF 420-425: %s\n', mat2str(sel(420,425)));
fprintf('PROF 426-433: %s\n', mat2str(sel(426,433)));
fprintf('PROF 434-440: %s\n', mat2str(sel(434,440)));
fprintf('PROF 441-450: %s\n', mat2str(sel(441,450)));
fprintf('PROF 490-510: %s\n', mat2str(sel(490,510)));
