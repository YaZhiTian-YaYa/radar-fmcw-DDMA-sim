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
targetsOverride = S.scenes(ii).tg;  seedOverride = S.scenes(ii).seed;
fprintf('\n########## RERUN %d/%d [%s] seed %d ##########\n', ...
    ii, numel(S.scenes), S.scenes(ii).name, S.scenes(ii).seed);
try
    run('F:\Work3\CH_NEW_95_射频_DDMA\main_sim.m');
catch ME
    batchState.okRun = false;
    fprintf('RERUN RUNFAIL [%s]: %s\n', S.scenes(batchState.idx).name, ME.message);
end
S = load('C:\Users\NJ\AppData\Local\Temp\probe95\rerunD_state.mat');   % clearvars 后重载状态
if batchState.okRun
    r = matchRows(tgtOut, targets, dims, prm, batchState.idx, batchState.scnSeed);
    S.RES = [S.RES; r];
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
               ismember([tgtOut.peak_type], [3 4]));
    if isempty(idx), continue; end
    cand = idx(~claimed(idx));
    if isempty(cand), cand = idx; end
    nrr = cand(abs([tgtOut(cand).azimuth] - targets(g).azi) <= 1.0);
    if gp(g) && ~isempty(nrr)
        prr = nrr(ismember([tgtOut(nrr).splitTag], 1:4));
        if ~isempty(prr), nrr = prr; end
    end
    if ~isempty(nrr)
        dAz = abs([tgtOut(nrr).azimuth] - targets(g).azi);
        rv  = [tgtOut(nrr).rcsDb];
        [~, im] = min(dAz(:).' - 1e-6*rv(:).');
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
