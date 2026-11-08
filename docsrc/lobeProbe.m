% 主峰裂瓣量化: 8 根码线各自 ±12bin 内高于 门限 的旁瓣清单(位置/相对线高)
diary('C:\Users\NJ\AppData\Local\Temp\probe95\lobeProbe.log'); clc; diary on;
addpath('F:\Work3\CH_NEW_95_射频_DDMA');
saveFigDir = 'C:\Users\NJ\AppData\Local\Temp\probe95\cohfigs';
if ~exist(saveFigDir,'dir'), mkdir(saveFigDir); end
batchState.idx = 1;  batchState.okRun = true;
run('F:\Work3\CH_NEW_95_射频_DDMA\main_sim.m');
rxI = 7;
rowT = squeeze(20*log10(abs(fftCfg.data_sva(144, rxI, :)) + eps));
rE = 350; rowE = squeeze(20*log10(abs(fftCfg.data_sva(rE, rxI, :)) + eps));
fMed = median(rowE);
[~, pk] = maxk(rowT, 8);
fprintf('LB 噪底中位 %7.2f | 8线: %s\n', fMed, mat2str(sort(pk-1)));
nLobe = 0;
for k = 1:8
    lo = pk(k)-12; hi = pk(k)+12;
    idx = mod(lo:hi, 256)+1; seg = rowT(idx); segOff = mod(lo:hi, 256) - (pk(k)-1);
    seg(abs(segOff) <= 1) = -Inf;                     % 主瓣±1 除外
    thr = fMed + 6;
    cand = find(seg > thr);
    for q = 1:numel(cand)
        b0 = mod(lo-1+cand(q), 256);
        nLobe = nLobe + 1;
        fprintf('LB 线%d@%d 的瓣 %+3dbin: %7.2f (线%+6.1f dB, 高于门限%+5.1f)\n', ...
            k, pk(k)-1, segOff(cand(q)), seg(cand(q)), seg(cand(q))-rowT(pk(k)), seg(cand(q))-thr);
    end
end
fprintf('LB 合计裂瓣 %d 个(门限=噪底+6dB, 主瓣±1除外)\n', nLobe);
diary off;
