% D 批重设计重跑(2026-11-04): DML 可靠窗 = 强弱对(母10~12/弱-4~4)@0.8~1.2°
% + O04 修正(5dBsm@414 违反 ncSNR>0 准则 → 371 门); 匹配器改最近方位优先, 增列 isPair/点索引
diary('C:\Users\NJ\AppData\Local\Temp\probe95\rerunD.log'); clc; diary on;
addpath('F:\Work3\CH_NEW_95_射频_DDMA');
saveFigDir = 'C:\Users\NJ\AppData\Local\Temp\probe95\cohfigs';
if ~exist(saveFigDir, 'dir'), mkdir(saveFigDir); end
T = @(gR,gD,az,el,rcs,ph) struct('gate_R',gR,'gate_D',gD,'azi',az,'ele',el, ...
    'rcsDb',rcs,'rcsPhDeg',ph,'vTan',0,'isCorner',0);
C = struct('gate_R',5,'gate_D',0,'azi',0,'ele',0,'rcsDb',-5,'rcsPhDeg',0,'vTan',0,'isCorner',1);
sc = struct('name',{},'seed',{},'tg',{});  mi = 0;
% ---- O04 修正(分离度保持 12°, 371 门≈260m 保 ncSNR>0) ----
mi=mi+1; sc(mi).name='O04'; sc(mi).seed=4;  sc(mi).tg=[C, T(371,20,-6,0,5,0), T(371,20,6,0,10,60), T(257,-50,10,0,10,90)];
% ---- D 批: 强弱对 0.8~1.2°(探针确认窗) + 等强 0.8° 融合回退演示 ----
mi=mi+1; sc(mi).name='D01'; sc(mi).seed=4;  sc(mi).tg=[C, T(285,0,-0.5,0,10,0), T(285,0,0.5,0,10,90), T(257,40,-8,0,5,0)];   % 等强1.0°(已知可触发, 对照)
mi=mi+1; sc(mi).name='D02'; sc(mi).seed=4;  sc(mi).tg=[C, T(285,0,-0.5,0,10,0), T(285,0,0.5,0,-4,90), T(257,40,-8,0,5,0)];   % 1.0° 10/-4
mi=mi+1; sc(mi).name='D03'; sc(mi).seed=4;  sc(mi).tg=[C, T(257,30,-0.4,0,12,0), T(257,30,0.4,0,4,45), T(357,-10,5,0,10,0)]; % 0.8° 12/4
mi=mi+1; sc(mi).name='D04'; sc(mi).seed=4;  sc(mi).tg=[C, T(214,-20,-0.6,0,10,30), T(214,-20,0.6,0,-4,150), T(357,60,-6,0,12,0)]; % 1.2° 10/-4
mi=mi+1; sc(mi).name='D05'; sc(mi).seed=4;  sc(mi).tg=[C, T(214,40,-0.5,0,12,0), T(214,40,0.5,0,0,120), T(285,-30,10,0,10,0)];   % 1.0° 12/0
mi=mi+1; sc(mi).name='D06'; sc(mi).seed=4;  sc(mi).tg=[C, T(314,-60,-0.5,0,10,60), T(314,-60,0.5,0,2,240), T(200,20,12,0,10,0)];% 1.0° 10/2
mi=mi+1; sc(mi).name='D07'; sc(mi).seed=11; sc(mi).tg=[C, T(285,0,-0.5,0,10,0), T(285,0,0.5,0,-4,90), T(257,40,-8,0,5,0)];   % D02 换种子
mi=mi+1; sc(mi).name='D08'; sc(mi).seed=4;  sc(mi).tg=[C, T(171,20,-0.4,0,10,0), T(171,20,0.4,0,10,30)];                      % 等强0.8° 预期融合回退
batchState.scenes = sc;  batchState.RES = zeros(0, 15);  batchState.nOk = 0;  batchState.nFail = 0;
for ii = 1:numel(batchState.scenes)
    batchState.idx = ii;
    targetsOverride = batchState.scenes(ii).tg;  seedOverride = batchState.scenes(ii).seed;
    fprintf('\n########## RERUN %d/%d [%s] seed %d ##########\n', ...
        ii, numel(batchState.scenes), batchState.scenes(ii).name, batchState.scenes(ii).seed);
    batchState.okRun = true;
    try
        run('F:\Work3\CH_NEW_95_射频_DDMA\main_sim.m');
    catch ME
        batchState.okRun = false;
        fprintf('RERUN RUNFAIL [%s]: %s\n', batchState.scenes(batchState.idx).name, ME.message);
    end
    if batchState.okRun
        r = matchRows(tgtOut, targets, dims, prm, batchState.idx, batchState.scenes(batchState.idx).seed);
        batchState.RES = [batchState.RES; r]; %#ok<AGROW>
        batchState.nOk = batchState.nOk + 1;
    else
        batchState.nFail = batchState.nFail + 1;
    end
    fprintf('RERUN PROG %d/%d (ok %d fail %d)\n', batchState.idx, numel(batchState.scenes), batchState.nOk, batchState.nFail);
end
save('C:\Users\NJ\AppData\Local\Temp\probe95\rerunD_results.mat', 'batchState');
fprintf('RERUN ALL DONE ok %d fail %d\n', batchState.nOk, batchState.nFail);
diary off;

function rows = matchRows(tgtOut, targets, dims, prm, caseId, seedVal)
nT  = numel(tgtOut);
claimed = false(1, nT);
rows = nan(numel(targets), 15);
gp = false(numel(targets), 1);                     % 同距同速成员标记
for a = 1:numel(targets)
    for b = 1:numel(targets)
        if a ~= b && targets(a).gate_R == targets(b).gate_R && targets(a).gate_D == targets(b).gate_D
            gp(a) = true;
        end
    end
end
if ~isfield(prm, 'stepDir') || isempty(prm.stepDir), sd = 1; else, sd = prm.stepDir; end
for g = 1:numel(targets)
    dopC = mod(targets(g).gate_D + sd*prm.cr*prm.B/(dims.dopBinNum-1) * ...
               (2*targets(g).gate_R*3e8/(2*prm.B))/3e8 * dims.dopBinNum, dims.dopBinNum);
    idx = find(abs([tgtOut.rangIdx] - (targets(g).gate_R+1)) <= 2 & ...
               abs([tgtOut.dopplerIdx] - (round(dopC)+1)) <= 3 & ...
               abs([tgtOut.azimuth] - targets(g).azi) <= 3 & ...
               ismember([tgtOut.peak_type], [3 4]));
    if isempty(idx), continue; end
    cand = idx(~claimed(idx));
    if isempty(cand), cand = idx; end
    nrr = cand(abs([tgtOut(cand).azimuth] - targets(g).azi) <= 1.0);
    if gp(g) && ~isempty(nrr)                       % 配对成员: ±1°内优先认领角色点(tag1-4)
        prr = nrr(ismember([tgtOut(nrr).splitTag], 1:4));
        if ~isempty(prr), nrr = prr; end
    end
    if ~isempty(nrr)
        dAz = abs([tgtOut(nrr).azimuth] - targets(g).azi);
        rv  = [tgtOut(nrr).rcsDb];
        [~, im] = min(dAz(:).' - 1e-6*rv(:).');   % 最近方位优先, 并列取强
        t = nrr(im);
    else
        [~, im] = min(abs([tgtOut(cand).azimuth] - targets(g).azi));
        t = cand(im);
    end
    claimed(t) = true;
    rows(g, :) = [caseId, g, targets(g).rcsDb, double(tgtOut(t).rcsDb), ...
        double(tgtOut(t).rcsDb) - targets(g).rcsDb, double(tgtOut(t).splitTag), ...
        double(tgtOut(t).rcsConf), double(tgtOut(t).rcsDom), ...
        targets(g).gate_R*dims.DeltaR, targets(g).azi, targets(g).ele, seedVal, ...
        double(tgtOut(t).rcsDbAlpha) - targets(g).rcsDb, double(gp(g)), t];
end
for g = 1:numel(targets)                           % 未匹配行补基础列
    if isnan(rows(g, 4))
        rows(g, :) = [caseId, g, targets(g).rcsDb, NaN, NaN, NaN, NaN, NaN, ...
            targets(g).gate_R*dims.DeltaR, targets(g).azi, targets(g).ele, seedVal, NaN, double(gp(g)), NaN];
    end
end
end
