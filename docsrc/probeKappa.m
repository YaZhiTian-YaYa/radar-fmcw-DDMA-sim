% DML 分离质量取证: D01(等强1.0°触发) D02(10/-4融合) D04(10/-4@1.2°触发) D05(12/0触发但拉平)
% 打印配对簇的 角度/κ/封顶与未封顶权重/权重噪声σ/RCS
diary('C:\Users\NJ\AppData\Local\Temp\probe95\probeKappa.log'); clc; diary on;
addpath('F:\Work3\CH_NEW_95_射频_DDMA');
saveFigDir = 'C:\Users\NJ\AppData\Local\Temp\probe95\cohfigs';
if ~exist(saveFigDir, 'dir'), mkdir(saveFigDir); end
T = @(gR,gD,az,el,rcs,ph) struct('gate_R',gR,'gate_D',gD,'azi',az,'ele',el, ...
    'rcsDb',rcs,'rcsPhDeg',ph,'vTan',0,'isCorner',0);
C = struct('gate_R',5,'gate_D',0,'azi',0,'ele',0,'rcsDb',-5,'rcsPhDeg',0,'vTan',0,'isCorner',1);
sc = struct('name',{},'seed',{},'tg',{});  mi = 0;
mi=mi+1; sc(mi).name='D01'; sc(mi).seed=4;  sc(mi).tg=[C, T(285,0,-0.5,0,10,0), T(285,0,0.5,0,10,90), T(257,40,-8,0,5,0)];
mi=mi+1; sc(mi).name='D02'; sc(mi).seed=4;  sc(mi).tg=[C, T(285,0,-0.5,0,10,0), T(285,0,0.5,0,-4,90), T(257,40,-8,0,5,0)];
mi=mi+1; sc(mi).name='D04'; sc(mi).seed=4;  sc(mi).tg=[C, T(214,-20,-0.6,0,10,30), T(214,-20,0.6,0,-4,150), T(357,60,-6,0,12,0)];
mi=mi+1; sc(mi).name='D05'; sc(mi).seed=4;  sc(mi).tg=[C, T(214,40,-0.5,0,12,0), T(214,40,0.5,0,0,120), T(285,-30,10,0,10,0)];
batchState.scenes = sc;  batchState.idx = 0;
for ii = 1:numel(batchState.scenes)
    batchState.idx = ii;
    targetsOverride = batchState.scenes(ii).tg;  seedOverride = batchState.scenes(ii).seed;
    fprintf('\n########## KAPPA %d [%s] ##########\n', ii, batchState.scenes(ii).name);
    run('F:\Work3\CH_NEW_95_射频_DDMA\main_sim.m');
    ci = batchState.idx;  scN = batchState.scenes(ci).name;   % clearvars 后 ii 已清, 走 batchState
    tg = batchState.scenes(ci).tg;                 % 真值参考(跳过角反)
    for g = 2:numel(tg)
        fprintf('KAPPA %s 真值: tgt%d azi%+.1f° rcs %g dBsm (R%gm)\n', ...
            scN, g-1, tg(g).azi, tg(g).rcsDb, tg(g).gate_R*dims.DeltaR);
    end
    for t = 1:numel(tgtOut)
        pt = tgtOut(t);
        if ~isfield(pt, 'dmlKappa') || isnan(pt.dmlKappa), continue; end
        if isempty(find([abs(pt.rangIdx - ([tg.gate_R] + 1)) <= 3 & ...
                         abs(pt.azimuth - tg(2).azi) <= 6], 1)), continue; end  % 只看配对簇
        fprintf(['KAPPA %s 点#%d: azi%+7.2f° tag%d pass%d sep%.2f° | ' ...
                 'ang[%+.2f,%+.2f] κ=%.1f | W封顶[%6.2f,%6.2f] W全[%6.2f,%6.2f] ' ...
                 'σw[%5.2f,%5.2f] | rcs %6.2f\n'], ...
            scN, t, pt.azimuth, pt.splitTag, pt.dmlPass, ...
            pt.dmlSepDeg, pt.dmlAng1, pt.dmlAng2, pt.dmlKappa, ...
            pt.dmlW1Db, pt.dmlW2Db, pt.dmlW1DbU, pt.dmlW2DbU, ...
            20*log10(pt.dmlSigW1 + eps), 20*log10(pt.dmlSigW2 + eps), pt.rcsDb);
    end
end
fprintf('KAPPA DONE\n');
diary off;
