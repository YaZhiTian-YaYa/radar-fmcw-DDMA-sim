% 文档用最终验证批(2026-11-04): 三批 = 独立 / OMP对(分离度≥OMP窗) / DML对(亚OMP分离度)
% 全部目标 ncSNR>0; 结果存 .mat 供图表与 docx 生成
diary('C:\Users\NJ\AppData\Local\Temp\probe95\finalbatch.log'); clc; diary on;
addpath('F:\Work3\CH_NEW_95_射频_DDMA');
saveFigDir = 'C:\Users\NJ\AppData\Local\Temp\probe95\cohfigs';
if ~exist(saveFigDir, 'dir'), mkdir(saveFigDir); end
T = @(gR,gD,az,el,rcs,ph) struct('gate_R',gR,'gate_D',gD,'azi',az,'ele',el, ...
    'rcsDb',rcs,'rcsPhDeg',ph,'vTan',0,'isCorner',0);
C = struct('gate_R',5,'gate_D',0,'azi',0,'ele',0,'rcsDb',-5,'rcsPhDeg',0,'vTan',0,'isCorner',1);
sc = struct('name',{},'seed',{},'tg',{});  mi = 0;
% ---- 批1 独立目标(12 场景, 多 R/azi/ele/RCS/半bin/种子) ----
mi=mi+1; sc(mi).name='I01'; sc(mi).seed=4;  sc(mi).tg=[C, T(214,20,5,0,10,0), T(285,-40,-8,0,4,30), T(357,60,12,0,7,60)];
mi=mi+1; sc(mi).name='I02'; sc(mi).seed=4;  sc(mi).tg=[C, T(414,10,0,0,12,0), T(257,-60,-15,0,5,90), T(171,30,8,0,4,20)];
mi=mi+1; sc(mi).name='I03'; sc(mi).seed=4;  sc(mi).tg=[C, T(357,-20,3,0,15,0), T(414,50,-5,0,10,120), T(285,0,15,0,4,200)];
mi=mi+1; sc(mi).name='I04'; sc(mi).seed=4;  sc(mi).tg=[C, T(143,-10,-12,0,10,0), T(414,20,5,0,10,45), T(400,-50,-3,0,8,75)];
mi=mi+1; sc(mi).name='I05'; sc(mi).seed=4;  sc(mi).tg=[C, T(142.5,15,6,0,10,0), T(200.5,-25,-10,0,4,30), T(270.5,45,14,0,7,150)];  % 半bin
mi=mi+1; sc(mi).name='I06'; sc(mi).seed=11; sc(mi).tg=[C, T(214,20,5,0,10,0), T(285,-40,-8,0,4,30), T(357,60,12,0,7,60)];
mi=mi+1; sc(mi).name='I07'; sc(mi).seed=23; sc(mi).tg=[C, T(414,10,0,0,12,0), T(257,-60,-15,0,5,90), T(171,30,8,0,4,20)];
mi=mi+1; sc(mi).name='I08'; sc(mi).seed=35; sc(mi).tg=[C, T(142.5,15,6,0,10,0), T(200.5,-25,-10,0,4,30), T(270.5,45,14,0,7,150)];
mi=mi+1; sc(mi).name='I09'; sc(mi).seed=4;  sc(mi).tg=[C, T(214,30,70,0,15,0), T(214,-30,-70,0,15,45), T(257,50,40,0,10,90)];    % 宽角
mi=mi+1; sc(mi).name='I10'; sc(mi).seed=4;  sc(mi).tg=[C, T(214,30,0,12,10,0), T(257,-30,0,-10,10,45), T(285,50,5,8,10,90)];    % 俯仰
mi=mi+1; sc(mi).name='I11'; sc(mi).seed=4;  sc(mi).tg=[C, T(357,60,-12,0,7,30), T(400,-20,8,0,9,0)];                             % 弱目标带
mi=mi+1; sc(mi).name='I12'; sc(mi).seed=47; sc(mi).tg=[C, T(285,0,0,0,10,0)];                                                    % 种子散布
% ---- 批2 OMP 对(8 场景, 分离度 2.5°~6°) ----
mi=mi+1; sc(mi).name='O01'; sc(mi).seed=4;  sc(mi).tg=[C, T(285,40,-5,0,12,45), T(285,40,5,0,2,210), T(400,30,-15,0,8,0)];
mi=mi+1; sc(mi).name='O02'; sc(mi).seed=4;  sc(mi).tg=[C, T(357,0,-3,0,10,0), T(357,0,3,0,10,90), T(214,-20,8,0,5,30)];
mi=mi+1; sc(mi).name='O03'; sc(mi).seed=4;  sc(mi).tg=[C, T(214,-30,-4,0,4,20), T(214,-30,4,0,10,140), T(414,10,0,0,10,0)];
mi=mi+1; sc(mi).name='O04'; sc(mi).seed=4;  sc(mi).tg=[C, T(414,20,-6,0,5,0), T(414,20,6,0,10,60), T(257,-50,10,0,10,90)];
mi=mi+1; sc(mi).name='O05'; sc(mi).seed=4;  sc(mi).tg=[C, T(314,60,-2.5,0,10,0), T(314,60,2.5,0,4,120), T(357,-35,-12,0,7,30)];
mi=mi+1; sc(mi).name='O06'; sc(mi).seed=4;  sc(mi).tg=[C, T(257,-45,-5,0,4,0), T(257,-45,5,0,10,90), T(371,25,6,0,8,45)];
mi=mi+1; sc(mi).name='O07'; sc(mi).seed=11; sc(mi).tg=[C, T(285,40,-5,0,12,45), T(285,40,5,0,2,210), T(400,30,-15,0,8,0)];
mi=mi+1; sc(mi).name='O08'; sc(mi).seed=4;  sc(mi).tg=[C, T(357,0,-4,0,12,0), T(357,0,4,0,0,90), T(214,30,10,0,10,0)];
% ---- 批3 DML 对(亚 OMP 分离度, 依探针窗口; 分离度 0.6°~1.2°) ----
mi=mi+1; sc(mi).name='D01'; sc(mi).seed=4;  sc(mi).tg=[C, T(285,0,-0.5,0,10,0), T(285,0,0.5,0,10,90), T(257,40,-8,0,5,0)];
mi=mi+1; sc(mi).name='D02'; sc(mi).seed=4;  sc(mi).tg=[C, T(257,30,-0.4,0,10,0), T(257,30,0.4,0,10,45), T(414,-10,5,0,7,90)];
mi=mi+1; sc(mi).name='D03'; sc(mi).seed=4;  sc(mi).tg=[C, T(285,-20,-0.6,0,10,0), T(285,-20,0.6,0,6,120), T(214,60,-6,0,12,0)];
mi=mi+1; sc(mi).name='D04'; sc(mi).seed=4;  sc(mi).tg=[C, T(414,10,-0.5,0,7,0), T(414,10,0.5,0,7,90), T(257,-30,10,0,10,30)];
mi=mi+1; sc(mi).name='D05'; sc(mi).seed=4;  sc(mi).tg=[C, T(214,40,-0.6,0,12,0), T(214,40,0.6,0,0,60), T(357,-60,-14,0,7,0)];
mi=mi+1; sc(mi).name='D06'; sc(mi).seed=4;  sc(mi).tg=[C, T(314,-60,-0.4,0,6,0), T(314,-60,0.4,0,10,90), T(400,20,12,0,10,45)];
mi=mi+1; sc(mi).name='D07'; sc(mi).seed=11; sc(mi).tg=[C, T(285,0,-0.5,0,10,0), T(285,0,0.5,0,10,90), T(257,40,-8,0,5,0)];
mi=mi+1; sc(mi).name='D08'; sc(mi).seed=4;  sc(mi).tg=[C, T(171,20,-0.5,0,10,0), T(171,20,0.5,0,6,30)];                            % 近距
batchState.scenes = sc;  batchState.RES = zeros(0, 13);  batchState.nOk = 0;  batchState.nFail = 0;
for ii = 1:numel(batchState.scenes)
    batchState.idx = ii;
    targetsOverride = batchState.scenes(ii).tg;  seedOverride = batchState.scenes(ii).seed;
    fprintf('\n########## FIN %d/%d [%s] seed %d ##########\n', ...
        ii, numel(batchState.scenes), batchState.scenes(ii).name, batchState.scenes(ii).seed);
    batchState.okRun = true;
    try
        run('F:\Work3\CH_NEW_95_射频_DDMA\main_sim.m');
    catch ME
        batchState.okRun = false;
        fprintf('FIN RUNFAIL [%s]: %s\n', batchState.scenes(batchState.idx).name, ME.message);
    end
    if batchState.okRun
        r = matchRows(tgtOut, targets, dims, prm, batchState.idx, batchState.scenes(batchState.idx).seed);
        batchState.RES = [batchState.RES; r]; %#ok<AGROW>
        batchState.nOk = batchState.nOk + 1;
    else
        batchState.nFail = batchState.nFail + 1;
    end
    fprintf('FIN PROG %d/%d (ok %d fail %d)\n', batchState.idx, numel(batchState.scenes), batchState.nOk, batchState.nFail);
end
save('C:\Users\NJ\AppData\Local\Temp\probe95\finalRCS_results.mat', 'batchState');
fprintf('FIN ALL DONE ok %d fail %d, 结果已存 finalRCS_results.mat\n', batchState.nOk, batchState.nFail);
diary off;

function rows = matchRows(tgtOut, targets, dims, prm, caseId, seedVal)
nT  = numel(tgtOut);
claimed = false(1, nT);
rows = nan(numel(targets), 13);
if ~isfield(prm, 'stepDir') || isempty(prm.stepDir), sd = 1; else, sd = prm.stepDir; end
for g = 1:numel(targets)
    dopC = mod(targets(g).gate_D + sd*prm.cr*prm.B/(dims.dopBinNum-1) * ...
               (2*targets(g).gate_R*3e8/(2*prm.B))/3e8 * dims.dopBinNum, dims.dopBinNum);
    rows(g, :) = [caseId, g, targets(g).rcsDb, NaN, NaN, NaN, NaN, NaN, ...
                  targets(g).gate_R*dims.DeltaR, targets(g).azi, targets(g).ele, seedVal, NaN];
    idx = find(abs([tgtOut.rangIdx] - (targets(g).gate_R+1)) <= 2 & ...
               abs([tgtOut.dopplerIdx] - (round(dopC)+1)) <= 3 & ...
               abs([tgtOut.azimuth] - targets(g).azi) <= 3 & ...
               ismember([tgtOut.peak_type], [3 4]));
    if isempty(idx), continue; end
    cand = idx(~claimed(idx));
    if isempty(cand), cand = idx; end
    nrr = cand(abs([tgtOut(cand).azimuth] - targets(g).azi) <= 1.0);
    if ~isempty(nrr)
        [~, im] = max([tgtOut(nrr).rcsDb]);  t = nrr(im);
    else
        [~, im] = min(abs([tgtOut(cand).azimuth] - targets(g).azi));  t = cand(im);
    end
    claimed(t) = true;
    rows(g, :) = [caseId, g, targets(g).rcsDb, double(tgtOut(t).rcsDb), ...
        double(tgtOut(t).rcsDb) - targets(g).rcsDb, double(tgtOut(t).splitTag), ...
        double(tgtOut(t).rcsConf), double(tgtOut(t).rcsDom), ...
        targets(g).gate_R*dims.DeltaR, targets(g).azi, targets(g).ele, seedVal, ...
        double(tgtOut(t).rcsDbAlpha) - targets(g).rcsDb];
end
end
