% 链路预算审计对账抓取: 场景由 auditScene 选择(1=默认50dBsm@300m, 2=10dBsm@200m)
diary('C:\Users\NJ\AppData\Local\Temp\probe95\auditCap.log'); clc; diary on;
addpath('F:\Work3\CH_NEW_95_射频_DDMA');
if auditScene == 2
    T = @(gR,gD,az,el,rcs,ph) struct('gate_R',gR,'gate_D',gD,'azi',az,'ele',el, ...
        'rcsDb',rcs,'rcsPhDeg',ph,'vTan',0,'isCorner',0);
    C = struct('gate_R',5,'gate_D',0,'azi',0,'ele',0,'rcsDb',-5,'rcsPhDeg',0,'vTan',0,'isCorner',1);
    targetsOverride = [C, T(285,0,0,0,10,0)];
    seedOverride = 4;
else
    clear targetsOverride seedOverride
end
saveFigDir = 'C:\Users\NJ\AppData\Local\Temp\probe95\cohfigs';
if ~exist(saveFigDir,'dir'), mkdir(saveFigDir); end
batchState.idx = 1;  batchState.okRun = true;
run('F:\Work3\CH_NEW_95_射频_DDMA\main_sim.m');
diary off;
