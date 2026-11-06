% 遗留② 多cell合成 A/B 验证: D01 几何多种子 + 分离度变体(等强对, 唯一可信分离 regime)
% 单 cell 权重已入生产(rcsDb), MC 版按 rcsDb + (W_Mc − W_Sc) 离线折算对比
T = @(gR,gD,az,el,rcs,ph) struct('gate_R',gR,'gate_D',gD,'azi',az,'ele',el, ...
    'rcsDb',rcs,'rcsPhDeg',ph,'vTan',0,'isCorner',0);
C = struct('gate_R',5,'gate_D',0,'azi',0,'ele',0,'rcsDb',-5,'rcsPhDeg',0,'vTan',0,'isCorner',1);
sc = struct('name',{},'seed',{},'tg',{},'cfarScale',{},'mc',{});  mi = 0;
mi=mi+1; sc(mi).name='M1'; sc(mi).seed=4;  sc(mi).cfarScale=[]; sc(mi).mc=false;
sc(mi).tg=[C, T(285,0,-0.50,0,10,0), T(285,0,0.50,0,10,90), T(257,40,-8,0,5,0)];
mi=mi+1; sc(mi).name='M2'; sc(mi).seed=11; sc(mi).cfarScale=[]; sc(mi).mc=false;
sc(mi).tg=[C, T(285,0,-0.50,0,10,0), T(285,0,0.50,0,10,90), T(257,40,-8,0,5,0)];
mi=mi+1; sc(mi).name='M3'; sc(mi).seed=23; sc(mi).cfarScale=[]; sc(mi).mc=false;
sc(mi).tg=[C, T(285,0,-0.50,0,10,0), T(285,0,0.50,0,10,90), T(257,40,-8,0,5,0)];
mi=mi+1; sc(mi).name='M4'; sc(mi).seed=35; sc(mi).cfarScale=[]; sc(mi).mc=false;
sc(mi).tg=[C, T(285,0,-0.50,0,10,0), T(285,0,0.50,0,10,90), T(257,40,-8,0,5,0)];
mi=mi+1; sc(mi).name='M5'; sc(mi).seed=47; sc(mi).cfarScale=[]; sc(mi).mc=false;
sc(mi).tg=[C, T(285,0,-0.50,0,10,0), T(285,0,0.50,0,10,90), T(257,40,-8,0,5,0)];
mi=mi+1; sc(mi).name='M6'; sc(mi).seed=4;  sc(mi).cfarScale=[]; sc(mi).mc=false;
sc(mi).tg=[C, T(285,0,-0.45,0,10,0), T(285,0,0.45,0,10,90), T(257,40,-8,0,5,0)];   % 0.9°
mi=mi+1; sc(mi).name='M7'; sc(mi).seed=4;  sc(mi).cfarScale=[]; sc(mi).mc=false;
sc(mi).tg=[C, T(285,0,-0.55,0,10,0), T(285,0,0.55,0,10,90), T(257,40,-8,0,5,0)];   % 1.1°
scenes = sc;  RES = zeros(0, 15);  nOk = 0;  nFail = 0;  nextIdx = 1;
save('C:\Users\NJ\AppData\Local\Temp\probe95\rerunD_state.mat', 'scenes', 'RES', 'nOk', 'nFail', 'nextIdx');
fprintf('MCVAL INIT %d scenes\n', numel(scenes));
