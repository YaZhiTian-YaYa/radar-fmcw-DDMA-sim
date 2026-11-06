% 单场景驱动(进程隔离): bash 循环逐场景调用 matlab -batch, 每场景新进程, 规避累积堆损坏
% 状态: rerunD_state.mat = .scenes .RES .nOk .nFail .nextIdx
addpath('F:\Work3\CH_NEW_95_射频_DDMA');
S = load('C:\Users\NJ\AppData\Local\Temp\probe95\rerunD_state.mat');
ii = S.nextIdx;
if ii > numel(S.scenes)
    fprintf('RERUN1 ALL DONE ok %d fail %d\n', S.nOk, S.nFail);
    return;
end
diary('C:\Users\NJ\AppData\Local\Temp\probe95\rerunD.log'); diary on;
saveFigDir = 'C:\Users\NJ\AppData\Local\Temp\probe95\cohfigs';
if ~exist(saveFigDir, 'dir'), mkdir(saveFigDir); end
batchState.idx = ii;  batchState.okRun = true;
batchState.scnSeed = S.scenes(ii).seed;              % clearvars 前备好(run 后 S 已被抹)
batchState.scnName = S.scenes(ii).name;
if isfield(S.scenes(ii), 'cfarScale'), cfarFarScale = S.scenes(ii).cfarScale; else, cfarFarScale = []; end
if isfield(S.scenes(ii), 'mc'), mcPairWOverride = logical(S.scenes(ii).mc); else, mcPairWOverride = false; end
if isfield(S.scenes(ii), 'azDec'), azDecoupOverride = logical(S.scenes(ii).azDec); else, azDecoupOverride = true; end
if isfield(S.scenes(ii), 'win'), winTypeOverride = S.scenes(ii).win; else, winTypeOverride = []; end
targetsOverride = S.scenes(ii).tg;  seedOverride = S.scenes(ii).seed;
fprintf('\n########## RERUN %d/%d [%s] seed %d ##########\n', ...
    ii, numel(S.scenes), S.scenes(ii).name, S.scenes(ii).seed);
try
    run('F:\Work3\CH_NEW_95_射频_DDMA\main_sim.m');
catch ME
    batchState.okRun = false;
    fprintf('RERUN RUNFAIL [%s]: %s\n', batchState.scnName, ME.message);
end
S = load('C:\Users\NJ\AppData\Local\Temp\probe95\rerunD_state.mat');   % clearvars 后重载状态
if batchState.okRun
    farN = sum([tgtOut.rangIdx] >= 300);            % 远场云点数(假警观察, bin300≈210m)
    fprintf('THRSW %s: 云点%d 远场点(rangIdx≥300)=%d\n', ...
        S.scenes(batchState.idx).name, numel(tgtOut), farN);
end
if batchState.okRun && isfield(tgtOut, 'dmlImprDb')
    ip = find(~isnan([tgtOut.dmlImprDb]));
    pv = ip([tgtOut(ip).dmlPass] == 1 | ismember([tgtOut(ip).splitTag], [3 4]));
    np = setdiff(ip, pv);
    o2 = [];
    if ~isempty(np)
        [~, o2] = maxk([tgtOut(np).dmlImprDb], min(2, numel(np)));  % 未通过中假改善最高 2 点
    end
    for t = [pv, np(o2)]
        pt = tgtOut(t);
        fprintf(['DMLDIAG %s #%d: ang[%+.2f,%+.2f] sep%.2f pass%d impr%.2f κG=%.1f resPm%.1f imprS%.2f | ' ...
                 'W[%6.2f,%6.2f] Mc[%6.2f,%6.2f] Keff%.2f\n'], S.scenes(batchState.idx).name, t, ...
            pt.dmlAng1, pt.dmlAng2, pt.dmlSepDeg, pt.dmlPass, pt.dmlImprDb, pt.dmlKappa, ...
            pt.dmlResPm, pt.dmlImprDbS, pt.dmlW1Db, pt.dmlW2Db, pt.dmlW1DbMc, pt.dmlW2DbMc, pt.dmlKeff);
        fprintf('DMLALT %s #%d: alt ang[%+.2f,%+.2f] sep%.2f W[%6.2f,%6.2f] imprS%.2f\n', ...
            S.scenes(batchState.idx).name, t, pt.dmlAltAng(1), pt.dmlAltAng(2), ...
            abs(pt.dmlAltAng(1)-pt.dmlAltAng(2)), pt.dmlAltWDb(1), pt.dmlAltWDb(2), pt.dmlAltImprS);
    end
end
if batchState.okRun
    r = matchRows(tgtOut, targets, dims, prm, batchState.idx, batchState.scnSeed);
    S.RES = [S.RES; r];
    % DOA 偏差场测量(遗留①): 真值 vs 点云 (方位, 俯仰) 逐目标打印
    for g = 1:numel(targets)
        if isnan(r(g, 15)), continue; end
        pt = tgtOut(r(g, 15));
        fprintf('DOADIAG %s tgt%d: 真(az%+6.1f, el%+5.1f) -> 估(az%+7.2f, el%+6.2f)  d(az%+6.2f, el%+5.2f)\n', ...
            S.scenes(batchState.idx).name, g, targets(g).azi, targets(g).ele, ...
            pt.azimuth, pt.elevation, pt.azimuth - targets(g).azi, pt.elevation - targets(g).ele);
    end
    S.nOk = S.nOk + 1;
else
    S.nFail = S.nFail + 1;
end
S.nextIdx = batchState.idx + 1;
save('C:\Users\NJ\AppData\Local\Temp\probe95\rerunD_state.mat', '-struct', 'S');
fprintf('RERUN PROG %d/%d (ok %d fail %d)\n', batchState.idx, numel(S.scenes), S.nOk, S.nFail);
diary off;

function rows = matchRows(tgtOut, targets, dims, prm, caseId, seedVal)
nT  = numel(tgtOut);
claimed = false(1, nT);
rows = nan(numel(targets), 15);
gp = false(numel(targets), 1);
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
               ismember([tgtOut.peak_type], [3 4]) & ...
               [tgtOut.nonCohSNR] > 0);                % 噪声点出局(与批验证口径一致; 否则最近方位规则会被 az 恰更近的噪声点偷走匹配, I09 实证)
    if isempty(idx), continue; end
    cand = idx(~claimed(idx));
    if isempty(cand), cand = idx; end
    nrr = cand(abs([tgtOut(cand).azimuth] - targets(g).azi) <= 1.0);
    if gp(g) && ~isempty(nrr)
        prr = nrr(ismember([tgtOut(nrr).splitTag], 1:4));
        if ~isempty(prr), nrr = prr; end
    end
    if ~isempty(nrr)                                % 强度过滤: 剔除比候选最强者弱 6dB 的病态副本
        rvW = [tgtOut(nrr).rcsDb];                  % (I09 实证: 同目标 straddle 副本 rcs 可差 21dB)
        nrrW = nrr(rvW >= max(rvW) - 6);
        if ~isempty(nrrW), nrr = nrrW; end
    end
    if ~isempty(nrr)
        dAz = abs([tgtOut(nrr).azimuth] - targets(g).azi);
        rv  = [tgtOut(nrr).rcsDb];
        if gp(g)
            [~, im] = min(dAz(:).' - 1e-6*rv(:).'); % 配对成员: 最近方位(角色点已由 prr 优先)
        else
            [~, im] = max(rv(:).' - 1e-3*dAz(:).'); % 独立目标: 最强副本(straddle 副本 rcs 可
        end                                          % 发散 21dB, 近距规则会被弱副本偷走)
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
for g = 1:numel(targets)
    if isnan(rows(g, 4))
        rows(g, :) = [caseId, g, targets(g).rcsDb, NaN, NaN, NaN, NaN, NaN, ...
            targets(g).gate_R*dims.DeltaR, targets(g).azi, targets(g).ele, seedVal, NaN, double(gp(g)), NaN];
    end
end
end
