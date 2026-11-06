% Oracle-J 检验(诊断层): 强弱对的真配置 J 是否低于搜索最优 bestJ —— 检验"起点遗漏"假设
% 若 J_true < bestJ: 五次证伪的地貌论断推翻, 病根=候选起点生成(网格量化残差峰被主瓣余项霸占)
diary('C:\Users\NJ\AppData\Local\Temp\probe95\diagJtrue.log'); clc; diary on;
addpath('F:\Work3\CH_NEW_95_射频_DDMA');
saveFigDir = 'C:\Users\NJ\AppData\Local\Temp\probe95\cohfigs';
if ~exist(saveFigDir, 'dir'), mkdir(saveFigDir); end
T = @(gR,gD,az,el,rcs,ph) struct('gate_R',gR,'gate_D',gD,'azi',az,'ele',el, ...
    'rcsDb',rcs,'rcsPhDeg',ph,'vTan',0,'isCorner',0);
C = struct('gate_R',5,'gate_D',0,'azi',0,'ele',0,'rcsDb',-5,'rcsPhDeg',0,'vTan',0,'isCorner',1);
targetsOverride = [C, T(214,-20,-0.6,0,10,30), T(214,-20,0.6,0,-4,150), T(357,60,-6,0,12,0)];  % C4
seedOverride = 4;
batchState.idx = 1;  batchState.okRun = true;
run('F:\Work3\CH_NEW_95_射频_DDMA\main_sim.m');
aziGrid = AntCaliData.aziSteerVecRange(1):AntCaliData.aziSteerVecRange(3):AntCaliData.aziSteerVecRange(2);
[~, aziRowBest] = max(sum(rowLut, 2));
chanAzi = find(rowLut(aziRowBest, :));
A_dml = conj(AntCaliData.aziSteerVec(chanAzi, :));
fitFun = @(th) steerFitAt(A_dml, aziGrid, th);
% 配对簇主点: 靠近真值 az=-0.6 的最强点(有 doaChanIn 行)
dAz = abs([tgtOut.azimuth] - (-0.6));
nRow = size(detCfg.doaChanIn, 1);
cand = find(dAz < 1.5 & (1:numel(tgtOut)) <= nRow);
[~, im] = max([tgtOut(cand).doaMag]);  t = cand(im);
x = double(detCfg.doaChanIn(t, chanAzi)).';
C0 = real(x' * x);
% Jtrue: 真角度对 [-0.6, +0.6](无细化, 真值即精确)
At = [fitFun(-0.6), fitFun(0.6)];
aT = inv2x2(At'*At) * (At'*x);
Jtrue = C0 - real(x' * At * aT);
% Jsearch: 搜索最优解的 J(用 dmlAng)
As = [fitFun(tgtOut(t).dmlAng1), fitFun(tgtOut(t).dmlAng2)];
aS = inv2x2(As'*As) * (As'*x);
Jsearch = C0 - real(x' * As * aS);
% J1(单原子, 亚bin)
bb = abs(A_dml' * x); [~, p1] = max(bb);
d1 = 0; if p1 > 1 && p1 < numel(bb), b3 = bb(p1-1:p1+1); d1 = max(min(0.5*(b3(1)-b3(3))/(b3(1)-2*b3(2)+b3(3)+eps), 0.5), -0.5); end
a1 = fitFun(aziGrid(p1) + d1);
J1s = C0 - abs(a1' * x)^2 / real(a1' * a1);
fprintf('JTRUE 点#%d: C0=%.3e | J1s(单原子亚bin)=%.3e | Jsearch[%+.2f,%+.2f]=%.3e | Jtrue[-0.6,+0.6]=%.3e\n', ...
    t, C0, J1s, tgtOut(t).dmlAng1, tgtOut(t).dmlAng2, Jsearch, Jtrue);
if Jtrue < Jsearch
    fprintf('JTRUE 结论: 真配置 J 更低 %.3e —— 地貌论断推翻, 病根=候选起点(真配置未被搜索尝试)\n', Jsearch - Jtrue);
    fprintf('JTRUE 对应净改善 imprS_true = %.2f dB (对比门限 1.0)\n', 10*log10(J1s / Jtrue));
    aTw = inv2x2(At'*At) * (At'*x);
    fprintf('JTRUE 真配置权重: [%6.2f, %6.2f] dB (真值 10/-4 → 两成员差应 ≈14dB)\n', ...
        20*log10(abs(aTw(1)) + eps), 20*log10(abs(aTw(2)) + eps));
else
    fprintf('JTRUE 结论: 搜索解仍更优 (ΔJ=%.3e) —— 地貌论断维持\n', Jsearch - Jtrue);
end
diary off;
