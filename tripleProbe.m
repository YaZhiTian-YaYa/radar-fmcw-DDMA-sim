%% ========================================================================
%  tripleProbe.m —— 同距同速 3 目标: OMP+DML 联合性能验证(结果落盘版)
%  场景: 角反(bin5) + 300m 目标(不动) + 新增同距同速 3 目标:
%    目标A/B: 角度差小(1.5~4°), RCS 高且差距小(预期 DML 分出)
%    目标C  : 角度相隔远(≥6°), RCS 小(预期 OMP 分出 —— 仿真 0/6, 见汇总)
%  每案例输出到 "同距同速3目标仿真\案例N": 37 张结果图 + rawdata.mat(ADC 数据
%  + 场景真值 + 全链参数, 含 readme) + 评分.txt(三目标命中明细)。
%  运行: matlab -batch "run('tripleProbe.m')"  (每案例 ~1.2 分钟)
%  ========================================================================
close all;
outRoot = fullfile(pwd, '同距同速3目标仿真');
nCase = 6;                                          % 案例数
seed0 = 3101;                                       % 随机种子起点(每案例 +1)
results = zeros(nCase, 3);                          % [A命中, B命中, C命中] 点云点数
c = 1;
while c <= nCase
    rng(seed0 + c - 1);
    % ---- 三目标随机参数(在约定框架内) ----
    gR = 100 + round(rand*200);                     % 距离 bin(100~300, 中远距)
    gD = round(rand*400);                           % 速度 bin(0~399, 折叠内)
    dAzi = 1.5 + rand*2.5;                          % A/B 角度差 1.5~4°(DML 域)
    azCsign = (rand < 0.5)*2 - 1;                   % C 在 A/B 一侧远处
    sepC = 7 + rand*15;                             % C 与 A/B 中心差 ≥7°(OMP 域)
    azM  = (2*rand-1) * 8;                          % A/B 中心角(±8°, ROI 内)
    azA  = azM - dAzi/2;  azB = azM + dAzi/2;
    azCfull = azM + azCsign*sepC;                   % C 的方位角
    rcsAB = 28 + rand*8;                            % A/B RCS 高(28~36 dBsm)
    dRcs  = rand*3 - 1.5;                           % A/B 差距小(±1.5 dB)
    rcsC  = rcsAB - (8 + rand*6);                   % C RCS 小(比 A/B 低 8~14 dB)

    targets3 = struct('gate_R', 5, 'gate_D', 0, 'azi', 0, 'ele', 0, ...
                      'rcsDb', -5, 'vTan', 0, 'isCorner', 1);
    targets3(2) = struct('gate_R', 428.57, 'gate_D', 0, 'azi', 0, 'ele', 0, ...
                         'rcsDb', 20, 'vTan', 0, 'isCorner', 0);
    targets3(3) = struct('gate_R', gR, 'gate_D', gD, 'azi', azA, 'ele', 0, ...
                         'rcsDb', rcsAB, 'vTan', 0, 'isCorner', 0);
    targets3(4) = struct('gate_R', gR, 'gate_D', gD, 'azi', azB, 'ele', 0, ...
                         'rcsDb', rcsAB + dRcs, 'vTan', 0, 'isCorner', 0);
    targets3(5) = struct('gate_R', gR, 'gate_D', gD, 'azi', azCfull, 'ele', 0, ...
                         'rcsDb', rcsC, 'vTan', 0, 'isCorner', 0);

    caseDir = fullfile(outRoot, sprintf('案例%d', c));
    if ~exist(caseDir, 'dir'), mkdir(caseDir); end
    targetsOverride = targets3;
    saveFigDir = caseDir;                            % 图存到案例文件夹
    fprintf('======== 案例 %d/%d: RD(%.0f,%.0f) A=%.1f°/%.0f B=%.1f°/%.0f C=%.1f°/%.0f dBsm ========\n', ...
        c, nCase, gR, gD, azA, rcsAB, azB, rcsAB+dRcs, azCfull, rcsC);
    % 全部状态存 batchState(main_sim 的 clearvars 只保留钩子变量)
    batchState = struct('caseIdx', c, 'nCase', nCase, 'seed0', seed0, 'results', results, ...
        'outRoot', outRoot, ...
        'truthCase', struct('gR', gR, 'gD', gD, 'azA', azA, 'azB', azB, 'azC', azCfull, ...
                            'rcsA', rcsAB, 'rcsB', rcsAB+dRcs, 'rcsC', rcsC, 'dAzi', dAzi, ...
                            'azCfull', azCfull, 'sepC', sepC));
    evalc('main_sim;');
    c = batchState.caseIdx;  nCase = batchState.nCase;  seed0 = batchState.seed0;
    results = batchState.results;  outRoot = batchState.outRoot;  truthCase = batchState.truthCase;
    caseDir = fullfile(outRoot, sprintf('案例%d', c));  % 从恢复的状态重建(main_sim 清掉了局部变量)

    % ---- rawdata.mat(ADC + 场景 + 参数 + readme) ----
    raw = struct();
    raw.adcRawData = adcRawData;  raw.prm = prm;  raw.dims = dims;
    raw.targets = targets;  raw.truth = truthCase;
    raw.encodeSeq = encodeSeq;  raw.ddmShift = ddmShift;
    raw.TX_A = TX_A;  raw.RX_A = RX_A;
    raw.antenna_t = antenna_t;  raw.antenna_r = antenna_r;
    raw.FFT_d0 = FFT_d0;  raw.lambda = lambda;  raw.rfInfo = rfInfo;
    raw.readme = sprintf([ ...
        '同距同速 3 目标案例 rawdata 存放规则(详细版见工程目录《案例RawData数据说明.docx》)\n' ...
        '1. adcRawData [1024×512×8] int16 口径实采样 counts(ADC 60 MHz, 16bit 量化+左移)。\n' ...
        '2. 场景: 角反(bin5) + 300m 目标(bin428.57) + 同距同速 3 目标(A/B/C)。\n' ...
        '   A/B 角差 %.2f°(DML 域), C 与中心隔 %.1f°(OMP 域); A %.1f / B %.1f / C %.1f dBsm。\n' ...
        '3. 三目标同 RD 单元: 谱线 = mod(gate_D − 2·gate_R, 512)(下跳频 cr=2, 每 bin 耦合 −2)。\n' ...
        '4. 其余字段与《案例RawData数据说明.docx》一致。\n'], ...
        truthCase.dAzi, truthCase.sepC, truthCase.rcsA, truthCase.rcsB, truthCase.rcsC);
    save(fullfile(caseDir, 'rawdata.mat'), '-struct', 'raw');

    % ---- 评分: 三真值角在点云中的命中点数(该单元内 |az−truth| ≤ 1°) ----
    sd = -1;  B_ = 214e6;  cr_ = 2;  N_ = 512;          % 与当前 main_sim 参数一致
    dopC = mod(truthCase.gD + sd*cr_*B_/(N_-1)*(2*truthCase.gR*3e8/(2*B_))/3e8*N_, N_);
    inCell = find(abs([tgtOut.rangIdx] - (truthCase.gR+1)) <= 1 & ...
                  abs([tgtOut.dopplerIdx] - (round(dopC)+1)) <= 3);
    azs = [tgtOut(inCell).azimuth];
    isO = [tgtOut(inCell).isOMP2nd];
    isD = [tgtOut(inCell).isDML2nd];
    truth = [truthCase.azA, truthCase.azB, truthCase.azC];
    nm = {'A(DML主)', 'B(DML二)', 'C(OMP二)'};
    fid = fopen(fullfile(caseDir, '评分.txt'), 'w', 'n', 'UTF-8');
    fprintf(fid, '案例 %d: RD(%.0f,%.0f) A=%.1f°/%.0f B=%.1f°/%.0f C=%.1f°/%.0f dBsm\n', ...
        c, truthCase.gR, truthCase.gD, truthCase.azA, truthCase.rcsA, ...
        truthCase.azB, truthCase.rcsB, truthCase.azCfull, truthCase.rcsC);
    for q = 1:3
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
    batchState.results = results;                    % 回写 batchState(下一轮 main_sim 后取回)
    close all;
    c = c + 1;
end
clear targetsOverride saveFigDir batchState
fprintf('======== 汇总(%d 案例 × 3 目标命中点数) ========\n', nCase);
for c = 1:nCase
    fprintf('  案例%d: A=%d B=%d C=%d\n', c, results(c,:));
end
fprintf('A/B(DML 域) 全命中案例: %d/%d | C(OMP 域) 全命中: %d/%d\n', ...
    sum(results(:,1) > 0 & results(:,2) > 0), nCase, sum(results(:,3) > 0), nCase);
fprintf('全部产物已存: %s\n', outRoot);
