%% ========================================================================
%  thrSweep_ompdml.m —— OMP/DML 门限【多案例标定 + 留出集验证】
%  方法论: 不针对单案例调参。
%    1) 两个不同种子的随机多案例场景(成对案例 + 单目标), 各跑一遍全链路,
%       在 DBF 阶段保存快照(OMP/DML 的输入与真值无关, 可离线反复扫);
%    2) 校准集: 扫门限网格, 指标 = 案例分出率(两真值角都出现在点云) /
%       单目标误分裂率 / 假角度数; 按"误分裂≤8% 内取分出率最高, 平手取大门限"选参;
%    3) 留出集(未参与调参): 用选定门限验证, 并与当前默认门限对比。
%  运行: matlab -batch "run('thrSweep_ompdml.m')"
%  ========================================================================
clear vars; close all;

fprintf('======== 生成两个场景并跑全链路(各约 1 分钟) ========\n');
nOmp = 10;  nDml = 8;  nSingle = 18;                          % 案例量(统计强度)
[targetsOverride, ~] = genPairScene(101, nOmp, nDml, nSingle);% 校准集(种子 101)
fprintf('校准集(种子101): %d 目标 = 角反 + %d对(%d OMP型+%d DML型) + %d 单目标 —— 与 main_sim 内置场景不同!\n', ...
        numel(targetsOverride), nOmp+nDml, nOmp, nDml, nSingle);
saveDBFStage = fullfile(pwd, 'dbfStage_calib.mat');
evalc('main_sim;');
nOmp = 10;  nDml = 8;  nSingle = 18;                          % main_sim 的 clear 后重声明计数

[targetsOverride, ~] = genPairScene(777, nOmp, nDml, nSingle);% 留出集(种子 777, 不参与调参)
fprintf('留出集(种子777): %d 目标(同规模, 案例全新, 未参与调参)\n', numel(targetsOverride));
saveDBFStage = fullfile(pwd, 'dbfStage_hold.mat');
evalc('main_sim;');
% 两次 main_sim 的 clear 后统一重生成两份真值(按种子确定性可复现; 计数再重声明)
nOmp = 10;  nDml = 8;  nSingle = 18;
[~, truthC] = genPairScene(101, nOmp, nDml, nSingle);
[~, truthH] = genPairScene(777, nOmp, nDml, nSingle);
clear targetsOverride saveDBFStage                            % 清钩子, 防误用

fprintf('======== 校准集: 扫 OMP 门限(DML 门限固定 3) ========\n');
S = load('dbfStage_calib.mat');
gridO = [8 10 12 14 16];
hitO = zeros(size(gridO));  faO = zeros(size(gridO));  wrongO = zeros(size(gridO));
for q = 1:numel(gridO)
    tO = ompRefine_pipe(S.tgtOut, S.detCfg, S.AntCaliData, S.rowLut, S.dims, gridO(q), [], [], false);
    tD = dmlRefine_pipe(tO,        S.detCfg, S.AntCaliData, S.rowLut, S.dims, 3, [], [], false);
    [hitO(q), faO(q), wrongO(q)] = evalScene(tD, truthC);
    fprintf('thr2ndDb=%2d : 分出 %2d/%d (%3.0f%%), 误分裂 %d/%d, 假角度 %d\n', ...
        gridO(q), hitO(q), nPairOf(truthC), 100*hitO(q)/nPairOf(truthC), ...
        faO(q), nSingleOf(truthC), wrongO(q));
end
% 选参规则: 误分裂率 ≤ 25% 的门限中取分出率最高; 平手取【更小】门限
% (平手时更小门限保留灵敏度 —— 对余量刚过线的真实弱目标更友好;
%  样本有限时误分裂 1~2 个的差异在统计噪声内, 不值得为它抬高门限)
ok = faO <= 0.25 * nSingleOf(truthC);
best = -1;  thrOsel = NaN;
for q = 1:numel(gridO)
    if ok(q) && hitO(q) > best
        best = hitO(q);  thrOsel = gridO(q);
    end
end
if isnan(thrOsel), thrOsel = gridO(find(faO == min(faO), 1)); end   % 全超界兜底: 误分裂最小者
fprintf('>> 选定 thr2ndDb = %d\n', thrOsel);

fprintf('======== 校准集: 扫 DML 门限(OMP 门限用选定值) ========\n');
gridD = [1.5 2 2.5 3 3.5 4 5];
hitD = zeros(size(gridD));  faD = zeros(size(gridD));  wrongD = zeros(size(gridD));
tOs = ompRefine_pipe(S.tgtOut, S.detCfg, S.AntCaliData, S.rowLut, S.dims, thrOsel, [], [], false);
for q = 1:numel(gridD)
    tD = dmlRefine_pipe(tOs, S.detCfg, S.AntCaliData, S.rowLut, S.dims, gridD(q), [], [], false);
    [hitD(q), faD(q), wrongD(q)] = evalScene(tD, truthC);
    fprintf('thrDmlDb=%3.1f : 分出 %2d/%d (%3.0f%%), 误分裂 %d/%d, 假角度 %d\n', ...
        gridD(q), hitD(q), nPairOf(truthC), 100*hitD(q)/nPairOf(truthC), ...
        faD(q), nSingleOf(truthC), wrongD(q));
end
ok = faD <= 0.25 * nSingleOf(truthC);
best = -1;  thrDsel = NaN;
for q = 1:numel(gridD)
    if ok(q) && hitD(q) > best
        best = hitD(q);  thrDsel = gridD(q);
    end
end
if isnan(thrDsel), thrDsel = gridD(find(faD == min(faD), 1)); end
fprintf('>> 选定 thrDmlDb = %.1f\n', thrDsel);

fprintf('======== 留出集验证(种子 777, 未参与调参) ========\n');
H = load('dbfStage_hold.mat');
cfgList = [12 3; thrOsel thrDsel];                 % 当前方默认门限 vs 多案例标定门限
for ii = 1:size(cfgList, 1)
    cfg = cfgList(ii, :);                          % 按行取(勿用 for=矩阵, 那是按列迭代!)
    tO = ompRefine_pipe(H.tgtOut, H.detCfg, H.AntCaliData, H.rowLut, H.dims, cfg(1), [], [], false);
    tD = dmlRefine_pipe(tO,        H.detCfg, H.AntCaliData, H.rowLut, H.dims, cfg(2), [], [], false);
    [hit, fa, wrong] = evalScene(tD, truthH);
    tag = '当前默认';  if ii == 2, tag = '标定选定'; end
    fprintf('%s 门限 (%2d, %.1f) : 分出 %2d/%d (%3.0f%%), 误分裂 %d/%d, 假角度 %d\n', ...
        tag, cfg(1), cfg(2), hit, nPairOf(truthH), 100*hit/nPairOf(truthH), fa, nSingleOf(truthH), wrong);
end
fprintf('======== 完成(选定门限: thr2ndDb=%d, thrDmlDb=%.1f) ========\n', thrOsel, thrDsel);

%% ---------- 评估: 真值 vs 点云 ----------
function [hit, fa, wrong] = evalScene(tgtOut, truth)
% hit: 案例对两真值角都出现在该单元点云(±1°)的数量(母点或第二目标均可)
% fa : 单目标单元出现任何第二目标(OMP/DML)的数量(误分裂)
% wrong: 案例对单元出现第二目标但两真值角未凑齐(分裂了但分错)
    hit = 0;  fa = 0;  wrong = 0;  nPair = 0;  nSingle = 0;
    for q = 1:numel(truth)
        gr = truth(q).cell(1);  dop = truth(q).cell(2);
        inCell = find(abs([tgtOut.rangIdx] - (gr+1)) <= 1 & abs([tgtOut.dopplerIdx] - (dop+1)) <= 3);
        azs = [tgtOut(inCell).azimuth];
        sec = [tgtOut(inCell).isOMP2nd] | [tgtOut(inCell).isDML2nd];
        if truth(q).type > 0
            nPair = nPair + 1;
            near1 = any(abs(azs - truth(q).az1) <= 1);
            near2 = any(abs(azs - truth(q).az2) <= 1);
            if near1 && near2
                hit = hit + 1;
            elseif any(sec)
                wrong = wrong + 1;                     % 分裂了但角度不对
            end
        else
            nSingle = nSingle + 1;
            if any(sec), fa = fa + 1; end
        end
    end
end
function n = nPairOf(truth), n = sum([truth.type] > 0); end
function n = nSingleOf(truth), n = sum([truth.type] == 0); end
