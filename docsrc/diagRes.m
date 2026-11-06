% 导向模型残差解剖(头号专项第一步): 残差 r = x − a(az_true)·α 的结构判定
% 通道结构(cCal 残差, 两目标 r 模式相同) vs 角结构(泄漏/插值, 随角度变)
diary('C:\Users\NJ\AppData\Local\Temp\probe95\diagRes.log'); clc; diary on;
addpath('F:\Work3\CH_NEW_95_射频_DDMA');
saveFigDir = 'C:\Users\NJ\AppData\Local\Temp\probe95\cohfigs';
if ~exist(saveFigDir, 'dir'), mkdir(saveFigDir); end
T = @(gR,gD,az,el,rcs,ph) struct('gate_R',gR,'gate_D',gD,'azi',az,'ele',el, ...
    'rcsDb',rcs,'rcsPhDeg',ph,'vTan',0,'isCorner',0);
C = struct('gate_R',5,'gate_D',0,'azi',0,'ele',0,'rcsDb',-5,'rcsPhDeg',0,'vTan',0,'isCorner',1);
targetsOverride = [C, T(214,-20,-0.6,0,10,30), T(357,60,-6,0,12,0)];  % C2 强目标 + 参照
seedOverride = 4;
batchState.idx = 1;  batchState.okRun = true;
run('F:\Work3\CH_NEW_95_射频_DDMA\main_sim.m');
% ---- 解剖(clearvars 后这些量都在 main_sim 工作区) ----
aziGrid = AntCaliData.aziSteerVecRange(1):AntCaliData.aziSteerVecRange(3):AntCaliData.aziSteerVecRange(2);
[~, aziRowBest] = max(sum(rowLut, 2));
chanAzi = find(rowLut(aziRowBest, :));
A_dml = conj(AntCaliData.aziSteerVec(chanAzi, :));
trA = [targets(1).azi, targets(2).azi];  trR = [targets(1).gate_R, targets(2).gate_R];
fprintf('RES解剖: 真值 az%.2f(R%d) 与 az%.2f(R%d)\n', trA(1), trR(1), trA(2), trR(2));
fprintf('RES维度: A_dml %dx%d, aziGrid %d 点[%g..%g], chanAzi %d, 云点 %d\n', ...
    size(A_dml,1), size(A_dml,2), numel(aziGrid), aziGrid(1), aziGrid(end), numel(chanAzi), numel(tgtOut));
rPats = [];
for g = 1:2
    % 找该目标的点(角度最近且 ncSNR>0)
    dAz = abs([tgtOut.azimuth] - trA(g));
    cand = find(dAz < 1.0 & [tgtOut.nonCohSNR] > 0 & (1:numel(tgtOut)) <= size(detCfg.doaChanIn, 1));
    [~, im] = min(dAz(cand));  t = cand(im);
    x = double(detCfg.doaChanIn(t, chanAzi)).';
    aT = steerFitAt(A_dml, aziGrid, trA(g));            % 真角度原子
    al = (aT' * x) / (aT' * aT);
    r = x - aT * al;
    rNrm = norm(r) / norm(x);
    % 残差角度谱(前3峰): 看残差能量在哪些角度
    cR = 20*log10(abs(A_dml' * r) / max(abs(A_dml' * r) + eps));
    [pkV, pkI] = maxk(cR, 3);
    % 残差通道模式(归一化复向量, 供两目标对比)
    rP = r / (norm(r) + eps);
    rPats = [rPats, rP];
    % 投影损失: 用真角度原子 vs 用拟合角度原子 的 |α|
    aE = steerFitAt(A_dml, aziGrid, tgtOut(t).azimuth);
    alE = abs((aE' * x) / (aE' * aE));
    fprintf(['RES tgt%d: 点#%d az%+6.2f | 残差能量占比 %.3f | |α(真角)|%.4f |α(拟合角)|%.4f | ' ...
             '残差谱峰 @%+.1f°/%.1f dB, @%+.1f°/%.1f dB\n'], g, t, tgtOut(t).azimuth, rNrm^2, ...
             abs(al), alE, aziGrid(pkI(1)), pkV(1), aziGrid(pkI(2)), pkV(2));
    fprintf('     残差通道幅(dB, 相对r最大): %s\n', mat2str(round(20*log10(abs(r) / max(abs(r) + eps)), 0)).');
end
% 两目标残差通道模式的相关性: |<r1,r2>| —— 高=通道结构(cCal残差), 低=角结构
c12 = abs(rPats(:,1)' * rPats(:,2));
fprintf('RES 残差模式相关 |<r1,r2>| = %.3f (1=纯通道结构, 0=无关/角结构)\n', c12);
% 噪声底参考: 噪声在36通道上的能量 vs 残差能量
snrC = tgtOut(t).nonCohSNR;    % 最后目标的 ncSNR(参考量级)
fprintf('RES 参照: 最后点 ncSNR %.1f dB → 噪声占比预期 ~%.3f\n', snrC, 10^(-snrC/10));
diary off;
