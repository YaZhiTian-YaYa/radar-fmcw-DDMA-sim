% 噪底平坦度检查: 全程逐距离门噪底的波动统计 + 分段定位(泄漏区/目标区/远场)
diary('C:\Users\NJ\AppData\Local\Temp\probe95\flatChk.log'); clc; diary on;
addpath('F:\Work3\CH_NEW_95_射频_DDMA');
saveFigDir = 'C:\Users\NJ\AppData\Local\Temp\probe95\cohfigs';
if ~exist(saveFigDir,'dir'), mkdir(saveFigDir); end
batchState.idx = 1;  batchState.okRun = true;
run('F:\Work3\CH_NEW_95_射频_DDMA\main_sim.m');
fl = 10*log10(mean(noiseCfg.floorPowCh, 2) + eps);  fl = fl(:);
nR = numel(fl);
seg = @(a,b) fl(a:b);
medAll = median(fl(30:end));                     % 远场参考(排泄漏区)
% 分段: 泄漏区 1-30, 中场 30-130, 目标区 130-160, 远场 160-end
fprintf('FL VCO支路: %s | 远场参考中位 %7.2f\n', ...
    char(exist('enableMeasShoulder','var')*'(见打印)'+' '), medAll); %#ok<NASGU>
segDef = [30 130; 130 160; 160 min(300,nR); 300 nR];
for s = 1:size(segDef,1)
    a = segDef(s,1); b = min(segDef(s,2), nR);
    if b <= a, continue; end
    v = seg(a,b);
    fprintf('FL 段 %3d-%3d: 中位 %7.2f | 波动max %7.2f | min %7.2f | 峰峰 %5.2f dB\n', ...
        a, b, median(v), max(v), min(v), max(v)-min(v));
end
% 逐 bin 打印 120-170(目标区)与 300-360(纯远场) 的细节
fprintf('FL 目标区细节:\n');
for b = 130:10:170, fprintf('FL bin %3d: %7.2f (%+.2f)\n', b, fl(b), fl(b)-medAll); end
fprintf('FL 远场抽样(每10bin):\n');
for b = 300:10:400, fprintf('FL bin %3d: %7.2f (%+.2f)\n', b, fl(b), fl(b)-medAll); end
diary off;
