% I批抽查(当前代码假通过普遍性): I01/I04/I11 单目标常规几何
T = @(gR,gD,az,el,rcs,ph) struct('gate_R',gR,'gate_D',gD,'azi',az,'ele',el, ...
    'rcsDb',rcs,'rcsPhDeg',ph,'vTan',0,'isCorner',0);
C = struct('gate_R',5,'gate_D',0,'azi',0,'ele',0,'rcsDb',-5,'rcsPhDeg',0,'vTan',0,'isCorner',1);
sc = struct('name',{},'seed',{},'tg',{},'cfarScale',{});  mi = 0;
mi=mi+1; sc(mi).name='I01'; sc(mi).seed=4;  sc(mi).cfarScale=[];
sc(mi).tg=[C, T(214,20,5,0,10,0), T(285,-40,-8,0,4,30), T(357,60,12,0,7,60)];
mi=mi+1; sc(mi).name='I04'; sc(mi).seed=4;  sc(mi).cfarScale=[];
sc(mi).tg=[C, T(143,-10,-12,0,10,0), T(414,20,5,0,10,45), T(400,-50,-3,0,8,75)];
mi=mi+1; sc(mi).name='I11'; sc(mi).seed=4;  sc(mi).cfarScale=[];
sc(mi).tg=[C, T(357,60,-12,0,7,30), T(400,-20,8,0,9,0)];
scenes = sc;  RES = zeros(0, 15);  nOk = 0;  nFail = 0;  nextIdx = 1;
save('C:\Users\NJ\AppData\Local\Temp\probe95\rerunD_state.mat', 'scenes', 'RES', 'nOk', 'nFail', 'nextIdx');
fprintf('ISPOT INIT %d scenes\n', numel(scenes));
