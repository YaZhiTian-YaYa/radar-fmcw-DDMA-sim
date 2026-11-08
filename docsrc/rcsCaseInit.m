% RCS 案例复查初始化(3 场景: I01 多目标 / I12 单目标种子 / O04 配对, 与今晨批同定义)
T = @(gR,gD,az,el,rcs,ph) struct('gate_R',gR,'gate_D',gD,'azi',az,'ele',el, ...
    'rcsDb',rcs,'rcsPhDeg',ph,'vTan',0,'isCorner',0);
C = struct('gate_R',5,'gate_D',0,'azi',0,'ele',0,'rcsDb',-5,'rcsPhDeg',0,'vTan',0,'isCorner',1);
sc = struct('name',{},'seed',{},'tg',{});  mi = 0;
mi=mi+1; sc(mi).name='I01'; sc(mi).seed=4;  sc(mi).tg=[C, T(214,20,5,0,10,0), T(285,-40,-8,0,4,30), T(357,60,12,0,7,60)];
mi=mi+1; sc(mi).name='I12'; sc(mi).seed=47; sc(mi).tg=[C, T(285,0,0,0,10,0)];
mi=mi+1; sc(mi).name='O04'; sc(mi).seed=4;  sc(mi).tg=[C, T(414,20,-6,0,5,0), T(414,20,6,0,10,60), T(257,-50,10,0,10,90)];
scenes = sc;  RES = zeros(0, 15);  nOk = 0;  nFail = 0;  nextIdx = 1;
save('C:\Users\NJ\AppData\Local\Temp\probe95\rerunD_state.mat', 'scenes', 'RES', 'nOk', 'nFail', 'nextIdx');
fprintf('RCSCASE INIT %d scenes\n', numel(scenes));
