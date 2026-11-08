% 逐距离门噪底剖面(实测口径 fig24/fig5 同款): 目标±25bin 形状抽取
% 形状判据: 目标±2 = Hann 窗零点凹陷(纯热噪), ±8-12 = 弥散台座肩, 与 99 实测(144/152 凹, 136 肩)对比
diary('C:\Users\NJ\AppData\Local\Temp\probe95\profFloor.log'); clc; diary on;
addpath('F:\Work3\CH_NEW_95_射频_DDMA');
saveFigDir = 'C:\Users\NJ\AppData\Local\Temp\probe95\cohfigs';
if ~exist(saveFigDir,'dir'), mkdir(saveFigDir); end
batchState.idx = 1;  batchState.okRun = true;
run('F:\Work3\CH_NEW_95_射频_DDMA\main_sim.m');
fl = 10*log10(mean(noiseCfg.floorPowCh, 2) + eps);   % 逐距离门噪底(通道平均, dB)
fl = fl(:);
prow = mean(abs(fftCfg.data_sva(:, 7, :)).^2, 3);  prow = prow(:);
prow(1:12) = 0;
[~, r0] = max(prow);                                 % 目标行(自动, 同 obsShoulder)
base = median(fl(r0-25:r0-15));
fprintf('PF 目标行 %d (%.1f m), 基线(−25~−15bin 中位) %7.2f\n', r0, (r0-1)*dims.DeltaR, base);
for b = r0-14 : r0+14
    tag = '';
    if abs(b-r0)==2, tag = ' <- 凹(窗零点)?'; end
    if abs(b-r0)==8 || abs(b-r0)==10 || abs(b-r0)==12, tag = ' <- 肩?'; end
    fprintf('PF %+3d: %7.2f (%+.2f)%s\n', b-r0, fl(b), fl(b)-base, tag);
end
vL = fl(r0-2); vR = fl(r0+2);
fprintf('PF 汇总: 谷L(−2) %+.2f 谷R(+2) %+.2f | 肩L(−10) %+.2f 肩R(+10) %+.2f | 肩−谷 %+.2f / %+.2f\n', ...
    vL-base, vR-base, fl(r0-10)-base, fl(r0+10)-base, fl(r0-10)-vL, fl(r0+10)-vR);
diary off;
