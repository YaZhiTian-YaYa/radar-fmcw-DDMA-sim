%% ========================================================================
%  phaseSweepProbe.m —— 控制变量: 只扫 RCS 初相差的双近目标(多角度差 × 多初相)
%  前提(已验证): Δφ=0° 时 DBF 融合不可分 + DML 可分 → 基线成立。
%  本版: 3 档角差(1.5/2/3°) × 8 档初相(0~180°), 共 24 案例。
%  评分改进: O2 角度与真值逐个对比(不再只判"有没有", 还判"角度对不对")。
%  产物: "初相扫描双近目标\angX_phYYY"(图 + rawdata + 评分)。
%  运行: matlab -batch "run('phaseSweepProbe.m')"
%  ========================================================================
close all;
outRoot = fullfile(pwd, '初相扫描双近目标');

angList = [1.5, 2.0, 3.0];                        % 角差列表(deg, DML 域)
phList  = [0, 15, 30, 45, 60, 90, 135, 180];      % 初相差列表(deg)
nCase = numel(angList) * numel(phList);
gR = 180;  gD = 200;  azM = 3.0;  rcsPair = 30;   % 固定量

bs = struct('idx', 1, 'results', zeros(nCase, 5), 'outRoot', outRoot, ...
            'angList', angList, 'phList', phList, 'nCase', nCase, ...
            'gR', gR, 'gD', gD, 'azM', azM, 'rcsPair', rcsPair, 'phSweep', phList);
batchState = bs;

while bs.idx <= bs.nCase
    % ★ 每轮开头统一从 bs 恢复全部状态(main_sim 的 clearvars 会清一切局部变量)
    outRoot = bs.outRoot;  azM = bs.azM;  gR = bs.gR;  gD = bs.gD;
    rcsPair = bs.rcsPair;  nCase = bs.nCase;  phList = bs.phList;  angList = bs.angList;
    c = bs.idx;
    ia = ceil(c / numel(bs.phList));  ip = mod(c-1, numel(bs.phList)) + 1;
    dAzi = bs.angList(ia);  phB = bs.phList(ip);
    azA = bs.azM - dAzi/2;  azB = bs.azM + dAzi/2;

    targets2 = struct('gate_R', 5, 'gate_D', 0, 'azi', 0, 'ele', 0, ...
                      'rcsDb', -5, 'rcsPhDeg', 0, 'vTan', 0, 'isCorner', 1);
    targets2(2) = struct('gate_R', 428.57, 'gate_D', 0, 'azi', 0, 'ele', 0, ...
                         'rcsDb', 20, 'rcsPhDeg', 0, 'vTan', 0, 'isCorner', 0);
    targets2(3) = struct('gate_R', bs.gR, 'gate_D', bs.gD, 'azi', azA, 'ele', 0, ...
                         'rcsDb', bs.rcsPair, 'rcsPhDeg', 0, 'vTan', 0, 'isCorner', 0);
    targets2(4) = struct('gate_R', bs.gR, 'gate_D', bs.gD, 'azi', azB, 'ele', 0, ...
                         'rcsDb', bs.rcsPair, 'rcsPhDeg', phB, 'vTan', 0, 'isCorner', 0);

    caseDir = fullfile(outRoot, sprintf('ang%.1f_ph%03d', dAzi, phB));
    if ~exist(caseDir, 'dir'), mkdir(caseDir); end
    targetsOverride = targets2;
    saveFigDir = caseDir;
    fprintf('======== %d/%d: Δazi=%.1f° Δφ=%d° ========\n', c, nCase, dAzi, phB);
    evalc('main_sim;');
    % main_sim 后从 batchState 恢复(下面评分用的所有真值从 bs 派生)
    bs = batchState;
    c = bs.idx;
    dAzi = bs.angList(ceil(c / numel(bs.phList)));  phB = bs.phList(mod(c-1, numel(bs.phList)) + 1);
    azA = bs.azM - dAzi/2;  azB = bs.azM + dAzi/2;
    outRoot = bs.outRoot;
    caseDir = fullfile(outRoot, sprintf('ang%.1f_ph%03d', dAzi, phB));

    % ---- rawdata ----
    raw = struct('adcRawData', adcRawData, 'prm', prm, 'dims', dims, ...
                 'targets', targets, 'encodeSeq', encodeSeq, 'ddmShift', ddmShift, ...
                 'TX_A', TX_A, 'RX_A', RX_A, 'antenna_t', antenna_t, 'antenna_r', antenna_r, ...
                 'FFT_d0', FFT_d0, 'lambda', lambda, 'rfInfo', rfInfo, ...
                 'truth', struct('gR', bs.gR, 'gD', bs.gD, 'azA', azA, 'azB', azB, ...
                                 'dAzi', dAzi, 'phB', phB, 'rcsPair', bs.rcsPair));
    raw.readme = sprintf('角差 %.1f° 初相差 %d° 同 RD(%d,%d) RCS %g\n', dAzi, phB, bs.gR, bs.gD, bs.rcsPair);
    save(fullfile(caseDir, 'rawdata.mat'), '-struct', 'raw');

    % ---- 评分 ----
    dopC = mod(bs.gD + (-1)*2*214e6/511*(2*bs.gR*3e8/(2*214e6))/3e8*512, 512);
    inC = find(abs([tgtOut.rangIdx]-(bs.gR+1))<=1 & abs([tgtOut.dopplerIdx]-(round(dopC)+1))<=3);
    azs = [tgtOut(inC).azimuth];  isO = [tgtOut(inC).isOMP2nd];  isD = [tgtOut(inC).isDML2nd];

    % 各阶段: DBF(母点快照)能否分出两点分别靠近 A/B
    dbfSep = any(abs(stageDbf.azi - azA) <= 0.8) && any(abs(stageDbf.azi - azB) <= 0.8);
    % OMP: O2 点角度与真值 B 的偏差
    o2az = azs(isO == 1);
    if ~isempty(o2az)
        o2mean = mean(o2az); o2errB = abs(o2mean - azB);
    else
        o2mean = NaN; o2errB = NaN;
    end
    % DML: D2 点角度
    d2az = azs(isD == 1);
    if ~isempty(d2az), d2mean = mean(d2az); else, d2mean = NaN; end
    dmlSep = ~isempty(d2az);
    % 母点角度(DML 替换后)
    mAz = mean(azs(isO == 0 & isD == 0));
    % 最终: 两真值各有点
    hA = sum(abs(azs - azA) <= 1);  hB = sum(abs(azs - azB) <= 1);

    bs.results(c,:) = [hA, hB, o2errB, dbfSep, dmlSep];

    fid = fopen(fullfile(caseDir, '评分.txt'), 'w', 'n', 'UTF-8');
    fprintf(fid, '角差 %.1f° 初相差 %d°: A=%.1f° B=%.1f° RCS %g dBsm\n', dAzi, phB, azA, azB, bs.rcsPair);
    fprintf(fid, 'DBF 可分: %s | OMP O2 角: %s(误差B %s°) | DML D2 角: %s | 最终 A:%d B:%d\n\n', ...
        string(dbfSep), num2str(o2mean, '%.2f'), num2str(o2errB, '%.2f'), num2str(d2mean, '%.2f'), hA, hB);
    for q = 1:numel(inC)
        tag = 'M';  if isD(q)==1, tag='D2'; end;  if isO(q)==1, tag='O2'; end
        fprintf(fid, '  %s %7.2f°\n', tag, azs(q));
    end
    fclose(fid);
    fprintf('  DBF:%s OMP O2@%.2f°(errB %.2f°) DML D2@%.2f° | A:%d B:%d\n', ...
        string(dbfSep), o2mean, o2errB, d2mean, hA, hB);

    batchState = bs;
    close all;
    bs.idx = bs.idx + 1;
    batchState = bs;
end
clear targetsOverride saveFigDir batchState

fprintf('\n======== 汇总 ========\n');
fprintf('%6s %6s | %5s %5s | %8s %8s | %5s %5s %5s\n', ...
    'Dazi', 'Dphi', 'A', 'B', 'O2ang', 'O2errB', 'DBF', 'OMPok', 'DML');
for c = 1:bs.nCase
    ia = ceil(c / numel(bs.phList));  ip = mod(c-1, numel(bs.phList)) + 1;
    r = bs.results(c,:);
    fprintf('%6.1f %6d | %5d %5d | %8s %8s | %5s %5s %5s\n', ...
        angList(ia), phList(ip), r(1), r(2), ...
        num2str(r(3), '%.2f'), num2str(r(3), '%.2f'), ...
        string(r(4)==1), string(r(5)>0 || (r(3)<=1 && ~isnan(r(3)))), string(r(5)==1));
end
fprintf('产物: %s\n', outRoot);
