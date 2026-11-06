% 角色排序修复验证: D01(等强分裂) D13(强弱融合) O13(弱母强二) O01(正常对照)
T = @(gR,gD,az,el,rcs,ph) struct('gate_R',gR,'gate_D',gD,'azi',az,'ele',el, ...
    'rcsDb',rcs,'rcsPhDeg',ph,'vTan',0,'isCorner',0);
C = struct('gate_R',5,'gate_D',0,'azi',0,'ele',0,'rcsDb',-5,'rcsPhDeg',0,'vTan',0,'isCorner',1);
sc = struct('name',{},'seed',{},'tg',{},'cfarScale',{},'mc',{},'rob',{});  mi = 0;
mi=mi+1; sc(mi).name='D01'; sc(mi).seed=4;  sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(285,0,-0.5,0,10,0), T(285,0,0.5,0,10,90), T(257,40,-8,0,5,0)];
mi=mi+1; sc(mi).name='D13'; sc(mi).seed=4;  sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(285,0,-0.5,0,10,0), T(285,0,0.5,0,0,90), T(257,40,-8,0,5,0)];
mi=mi+1; sc(mi).name='O13'; sc(mi).seed=4;  sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(414,20,-6,0,5,0), T(414,20,6,0,10,60), T(257,-50,10,0,10,90)];
mi=mi+1; sc(mi).name='O01'; sc(mi).seed=4;  sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(285,40,-5,0,12,45), T(285,40,5,0,2,210), T(400,30,-15,0,8,0)];
scenes = sc;  RES = zeros(0, 15);  nOk = 0;  nFail = 0;  nextIdx = 1;
save('C:\Users\NJ\AppData\Local\Temp\probe95\rerunD_state.mat', 'scenes', 'RES', 'nOk', 'nFail', 'nextIdx');
fprintf('ROLE INIT %d scenes\n', numel(scenes));
