% 初始化进程隔离批跑状态(9 场景: O04 + D01~D08)
T = @(gR,gD,az,el,rcs,ph) struct('gate_R',gR,'gate_D',gD,'azi',az,'ele',el, ...
    'rcsDb',rcs,'rcsPhDeg',ph,'vTan',0,'isCorner',0);
C = struct('gate_R',5,'gate_D',0,'azi',0,'ele',0,'rcsDb',-5,'rcsPhDeg',0,'vTan',0,'isCorner',1);
sc = struct('name',{},'seed',{},'tg',{});  mi = 0;
mi=mi+1; sc(mi).name='O04'; sc(mi).seed=4;  sc(mi).tg=[C, T(371,20,-6,0,5,0), T(371,20,6,0,10,60), T(257,-50,10,0,10,90)];
mi=mi+1; sc(mi).name='D01'; sc(mi).seed=4;  sc(mi).tg=[C, T(285,0,-0.5,0,10,0), T(285,0,0.5,0,10,90), T(257,40,-8,0,5,0)];
mi=mi+1; sc(mi).name='D02'; sc(mi).seed=4;  sc(mi).tg=[C, T(285,0,-0.5,0,10,0), T(285,0,0.5,0,-4,90), T(257,40,-8,0,5,0)];
mi=mi+1; sc(mi).name='D03'; sc(mi).seed=4;  sc(mi).tg=[C, T(257,30,-0.4,0,12,0), T(257,30,0.4,0,4,45), T(357,-10,5,0,10,0)];
mi=mi+1; sc(mi).name='D04'; sc(mi).seed=4;  sc(mi).tg=[C, T(214,-20,-0.6,0,10,30), T(214,-20,0.6,0,-4,150), T(357,60,-6,0,12,0)];
mi=mi+1; sc(mi).name='D05'; sc(mi).seed=4;  sc(mi).tg=[C, T(214,40,-0.5,0,12,0), T(214,40,0.5,0,0,120), T(285,-30,10,0,10,0)];
mi=mi+1; sc(mi).name='D06'; sc(mi).seed=4;  sc(mi).tg=[C, T(314,-60,-0.5,0,10,60), T(314,-60,0.5,0,2,240), T(200,20,12,0,10,0)];
mi=mi+1; sc(mi).name='D07'; sc(mi).seed=11; sc(mi).tg=[C, T(285,0,-0.5,0,10,0), T(285,0,0.5,0,-4,90), T(257,40,-8,0,5,0)];
mi=mi+1; sc(mi).name='D08'; sc(mi).seed=4;  sc(mi).tg=[C, T(171,20,-0.4,0,10,0), T(171,20,0.4,0,10,30)];
scenes = sc;  RES = zeros(0, 15);  nOk = 0;  nFail = 0;  nextIdx = 1;
save('C:\Users\NJ\AppData\Local\Temp\probe95\rerunD_state.mat', 'scenes', 'RES', 'nOk', 'nFail', 'nextIdx');
fprintf('RERUN INIT %d scenes\n', numel(scenes));
