% 遗留③ thrSweep 门限灵敏度扫参: 远场粗门限缩放 × {1.00, 0.85, 0.72}
% 目标: 12dBsm@330m + 10dBsm@300m 检出; 近距-only 对照看远场假警
T = @(gR,gD,az,el,rcs,ph) struct('gate_R',gR,'gate_D',gD,'azi',az,'ele',el, ...
    'rcsDb',rcs,'rcsPhDeg',ph,'vTan',0,'isCorner',0);
C = struct('gate_R',5,'gate_D',0,'azi',0,'ele',0,'rcsDb',-5,'rcsPhDeg',0,'vTan',0,'isCorner',1);
sc = struct('name',{},'seed',{},'tg',{},'cfarScale',{});  mi = 0;
% 远弱目标场景 × 3 档门限
mi=mi+1; sc(mi).name='T100'; sc(mi).seed=4;  sc(mi).cfarScale=1.00;
sc(mi).tg=[C, T(471,30,10,0,12,0), T(428,-20,-8,0,10,30)];
mi=mi+1; sc(mi).name='T085'; sc(mi).seed=4;  sc(mi).cfarScale=0.85;
sc(mi).tg=[C, T(471,30,10,0,12,0), T(428,-20,-8,0,10,30)];
mi=mi+1; sc(mi).name='T072'; sc(mi).seed=4;  sc(mi).cfarScale=0.72;
sc(mi).tg=[C, T(471,30,10,0,12,0), T(428,-20,-8,0,10,30)];
mi=mi+1; sc(mi).name='T085s11'; sc(mi).seed=11; sc(mi).cfarScale=0.85;
sc(mi).tg=[C, T(471,30,10,0,12,0), T(428,-20,-8,0,10,30)];
% 近距-only 对照(远场无真值, 远场点=假警) × 3 档
mi=mi+1; sc(mi).name='E100'; sc(mi).seed=4;  sc(mi).cfarScale=1.00;
sc(mi).tg=[C, T(143,-10,-12,0,10,0)];
mi=mi+1; sc(mi).name='E085'; sc(mi).seed=4;  sc(mi).cfarScale=0.85;
sc(mi).tg=[C, T(143,-10,-12,0,10,0)];
mi=mi+1; sc(mi).name='E072'; sc(mi).seed=4;  sc(mi).cfarScale=0.72;
sc(mi).tg=[C, T(143,-10,-12,0,10,0)];
scenes = sc;  RES = zeros(0, 15);  nOk = 0;  nFail = 0;  nextIdx = 1;
save('C:\Users\NJ\AppData\Local\Temp\probe95\rerunD_state.mat', 'scenes', 'RES', 'nOk', 'nFail', 'nextIdx');
fprintf('THR INIT %d scenes\n', numel(scenes));
