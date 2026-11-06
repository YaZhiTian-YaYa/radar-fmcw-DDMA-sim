% 远旁瓣自投影实验场景(同 C 系列): C1/C2 单目标对照, C3/C4 强弱对, C5 等强对回归
T = @(gR,gD,az,el,rcs,ph) struct('gate_R',gR,'gate_D',gD,'azi',az,'ele',el, ...
    'rcsDb',rcs,'rcsPhDeg',ph,'vTan',0,'isCorner',0);
C = struct('gate_R',5,'gate_D',0,'azi',0,'ele',0,'rcsDb',-5,'rcsPhDeg',0,'vTan',0,'isCorner',1);
sc = struct('name',{},'seed',{},'tg',{},'cfarScale',{},'mc',{},'rob',{});  mi = 0;
mi=mi+1; sc(mi).name='C1'; sc(mi).seed=4;  sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=true;
sc(mi).tg=[C, T(214,40,-0.5,0,12,0), T(285,-30,10,0,10,0)];
mi=mi+1; sc(mi).name='C2'; sc(mi).seed=4;  sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=true;
sc(mi).tg=[C, T(214,-20,-0.6,0,10,30), T(357,60,-6,0,12,0)];
mi=mi+1; sc(mi).name='C3'; sc(mi).seed=4;  sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=true;
sc(mi).tg=[C, T(214,40,-0.5,0,12,0), T(214,40,0.5,0,0,120), T(285,-30,10,0,10,0)];
mi=mi+1; sc(mi).name='C4'; sc(mi).seed=4;  sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=true;
sc(mi).tg=[C, T(214,-20,-0.6,0,10,30), T(214,-20,0.6,0,-4,150), T(357,60,-6,0,12,0)];
mi=mi+1; sc(mi).name='C5'; sc(mi).seed=4;  sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=true;
sc(mi).tg=[C, T(285,0,-0.5,0,10,0), T(285,0,0.5,0,10,90), T(257,40,-8,0,5,0)];
scenes = sc;  RES = zeros(0, 15);  nOk = 0;  nFail = 0;  nextIdx = 1;
save('C:\Users\NJ\AppData\Local\Temp\probe95\rerunD_state.mat', 'scenes', 'RES', 'nOk', 'nFail', 'nextIdx');
fprintf('ROB INIT %d scenes\n', numel(scenes));
