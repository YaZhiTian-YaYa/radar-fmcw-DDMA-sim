%% ========================================================================
%  phasePairProbe.m —— 同距同速双目标(角度近 + RCS 初相差大): OMP 分离能力验证
%  条件: 两个目标角度相近(4~8°, OMP 域), RCS 幅度相同, 但初相差距大(90°/180°等)。
%  物理含义: 实际目标 RCS 是复数(σ∠φ), 两个同距同速目标即使幅度相同,
%  初相不同 → 阵列上的干涉图样完全不同 → DBF 融合峰偏移/残差结构变化。
%  场景: 保留角反(bin5) + 300m 目标(bin428.57), 只加这一对新目标。
%  每案例一组随机参数 + 随机初相组合, 跑全链路, 统计两真值角各几个点命中。
%  产物: "同距同速初相双目标\案例N"(图 + rawdata + 评分)。
%  运行: matlab -batch "run('phasePairProbe.m')"
%  ========================================================================
close all;
outRoot = fullfile(pwd, '同距同速初相双目标');
nCase = 5;                                          % 案例数(每案例一对目标)
seed0 = 4201;
results = zeros(nCase, 2);
c = 1;
while c <= nCase
    rng(seed0 + c - 1);
    % ---- 双目标随机参数 ----
    gR = 120 + round(rand*180);                     % 距离 bin(120~300)
    gD = round(rand*400);                           % 速度 bin
    dAzi = 4 + rand*4;                              % 角度差 4~8°(OMP 域, 与 DML 重叠)
    azM  = (2*rand-1) * 8;                          % 中心角(±8°, ROI 内)
    azA  = azM - dAzi/2;  azB = azM + dAzi/2;
    rcsPair = 28 + rand*8;                          % RCS(28~36 dBsm, 两目标相同)
    phList = [0, 90, 180, -90, 45, -135, 120, -60]; % 初相池(deg)
    phA = 0;                                        % 目标A 初相固定 0°(参考)
    phB = phList(randi(numel(phList)));             % 目标B 初相从池中随机
    if phB == 0, phB = 180; end                     % 避免 0/0(同相退化为常规场景)

    targets2 = struct('gate_R', 5, 'gate_D', 0, 'azi', 0, 'ele', 0, ...
                      'rcsDb', -5, 'rcsPhDeg', 0, 'vTan', 0, 'isCorner', 1);
    targets2(2) = struct('gate_R', 428.57, 'gate_D', 0, 'azi', 0, 'ele', 0, ...
                         'rcsDb', 20, 'rcsPhDeg', 0, 'vTan', 0, 'isCorner', 0);
    targets2(3) = struct('gate_R', gR, 'gate_D', gD, 'azi', azA, 'ele', 0, ...
                         'rcsDb', rcsPair, 'rcsPhDeg', phA, 'vTan', 0, 'isCorner', 0);
    targets2(4) = struct('gate_R', gR, 'gate_D', gD, 'azi', azB, 'ele', 0, ...
                         'rcsDb', rcsPair, 'rcsPhDeg', phB, 'vTan', 0, 'isCorner', 0);

    caseDir = fullfile(outRoot, sprintf('案例%d', c));
    if ~exist(caseDir, 'dir'), mkdir(caseDir); end
    targetsOverride = targets2;
    saveFigDir = caseDir;
    fprintf('======== 案例 %d/%d: RD(%.0f,%.0f) A=%.1f°/∠%d° B=%.1f°/∠%d° / %g dBsm ========\n', ...
        c, nCase, gR, gD, azA, phA, azB, phB, rcsPair);
    batchState = struct('caseIdx', c, 'nCase', nCase, 'seed0', seed0, 'results', results, ...
        'outRoot', outRoot, ...
        'truthCase', struct('gR', gR, 'gD', gD, 'azA', azA, 'azB', azB, ...
                            'phA', phA, 'phB', phB, 'rcsPair', rcsPair, 'dAzi', dAzi));
    evalc('main_sim;');
    c = batchState.caseIdx;  nCase = batchState.nCase;  seed0 = batchState.seed0;
    results = batchState.results;  outRoot = batchState.outRoot;  truthCase = batchState.truthCase;
    caseDir = fullfile(outRoot, sprintf('案例%d', c));

    % ---- rawdata.mat ----
    raw = struct();
    raw.adcRawData = adcRawData;  raw.prm = prm;  raw.dims = dims;
    raw.targets = targets;  raw.truth = truthCase;
    raw.encodeSeq = encodeSeq;  raw.ddmShift = ddmShift;
    raw.TX_A = TX_A;  raw.RX_A = RX_A;
    raw.antenna_t = antenna_t;  raw.antenna_r = antenna_r;
    raw.FFT_d0 = FFT_d0;  raw.lambda = lambda;  raw.rfInfo = rfInfo;
    raw.readme = sprintf([ ...
        '同距同速双目标(初相差大)案例 rawdata\n' ...
        '场景: 角反 + 300m 目标 + 同距同速双目标 A/B。\n' ...
        '  A=%.1f°/∠%d°, B=%.1f°/∠%d°, 角差 %.1f°(OMP 域), RCS %g dBsm(相同)。\n' ...
        '  初相差 = %d°(复数 RCS: 回波干涉图样随初相变)\n' ...
        '  谱线 = mod(gate_D − 2·gate_R, 512)(下跳频 cr=2)。\n'], ...
        truthCase.azA, truthCase.phA, truthCase.azB, truthCase.phB, ...
        truthCase.dAzi, truthCase.rcsPair, truthCase.phB - truthCase.phA);
    save(fullfile(caseDir, 'rawdata.mat'), '-struct', 'raw');

    % ---- 评分: 两真值角在点云中的命中点数 ----
    sd = -1;  B_ = 214e6;  cr_ = 2;  N_ = 512;
    dopC = mod(truthCase.gD + sd*cr_*B_/(N_-1)*(2*truthCase.gR*3e8/(2*B_))/3e8*N_, N_);
    inCell = find(abs([tgtOut.rangIdx] - (truthCase.gR+1)) <= 1 & ...
                  abs([tgtOut.dopplerIdx] - (round(dopC)+1)) <= 3);
    azs = [tgtOut(inCell).azimuth];
    isO = [tgtOut(inCell).isOMP2nd];
    isD = [tgtOut(inCell).isDML2nd];
    truth = [truthCase.azA, truthCase.azB];
    nm = {'A(强)', 'B(强)'};
    fid = fopen(fullfile(caseDir, '评分.txt'), 'w', 'n', 'UTF-8');
    fprintf(fid, '案例 %d: RD(%.0f,%.0f) A=%.1f°/∠%d° B=%.1f°/∠%d° / %g dBsm | 初相差 %d°\n', ...
        c, truthCase.gR, truthCase.gD, truthCase.azA, truthCase.phA, ...
        truthCase.azB, truthCase.phB, truthCase.rcsPair, truthCase.phB - truthCase.phA);
    for q = 1:2
        hit = abs(azs - truth(q)) <= 1;
        fprintf('  %s 真%.1f°: %d 点命中', nm{q}, truth(q), sum(hit));
        fprintf(fid, '  %s 真%.1f°: %d 点命中', nm{q}, truth(q), sum(hit));
        for r = find(hit)
            if     isD(r) == 1, tag = 'D2';
            elseif isO(r) == 1, tag = 'O2';
            else,               tag = 'M ';
            end
            fprintf(' [%s %.1f]', tag, azs(r));
            fprintf(fid, ' [%s %.1f]', tag, azs(r));
        end
        fprintf('\n');  fprintf(fid, '\n');
        results(c, q) = sum(hit);
    end
    fclose(fid);
    fprintf('案例 %d 完成: %s\n', c, caseDir);
    batchState.results = results;
    close all;
    c = c + 1;
end
clear targetsOverride saveFigDir batchState
fprintf('======== 汇总(%d 案例 × 2 目标命中点数) ========\n', nCase);
for c = 1:nCase
    fprintf('  案例%d: A=%d B=%d\n', c, results(c,:));
end
fprintf('双目标全命中案例: %d/%d\n', sum(results(:,1) > 0 & results(:,2) > 0), nCase);
fprintf('全部产物已存: %s\n', outRoot);
