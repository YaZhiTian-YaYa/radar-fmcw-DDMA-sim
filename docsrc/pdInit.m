% 亚分辨率 RD 邻近对探针: RD 不完全同 cell 但分辨率分不开, 只能靠角度分
% 等强 10/10 @±0.5°(sep1.0°) 变 RD 偏移 + 强弱 10/-4 加多普勒偏移(正交性悬念)
T = @(gR,gD,az,el,rcs,ph) struct('gate_R',gR,'gate_D',gD,'azi',az,'ele',el, ...
    'rcsDb',rcs,'rcsPhDeg',ph,'vTan',0,'isCorner',0);
C = struct('gate_R',5,'gate_D',0,'azi',0,'ele',0,'rcsDb',-5,'rcsPhDeg',0,'vTan',0,'isCorner',1);
sc = struct('name',{},'seed',{},'tg',{},'cfarScale',{},'mc',{},'rob',{});  mi = 0;
mi=mi+1; sc(mi).name='P0'; sc(mi).seed=4;  sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(285,0,-0.5,0,10,0), T(285,0,+0.5,0,10,90), T(257,40,-8,0,5,0)];          % 基线: 同RD(=D01)
mi=mi+1; sc(mi).name='P1'; sc(mi).seed=4;  sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(285,0,-0.5,0,10,0), T(285,0.5,+0.5,0,10,90), T(257,40,-8,0,5,0)];        % ΔD=0.5bin
mi=mi+1; sc(mi).name='P2'; sc(mi).seed=4;  sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(285,0,-0.5,0,10,0), T(285,1.0,+0.5,0,10,90), T(257,40,-8,0,5,0)];        % ΔD=1bin
mi=mi+1; sc(mi).name='P3'; sc(mi).seed=4;  sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(285,0,-0.5,0,10,0), T(285.5,0,+0.5,0,10,90), T(257,40,-8,0,5,0)];        % ΔR=0.5bin
mi=mi+1; sc(mi).name='P4'; sc(mi).seed=4;  sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(285,0,-0.5,0,10,0), T(285.5,0.5,+0.5,0,10,90), T(257,40,-8,0,5,0)];      % ΔR+ΔD 各半bin
mi=mi+1; sc(mi).name='P5'; sc(mi).seed=4;  sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(285,0,-0.5,0,10,0), T(285,0.5,+0.5,0,-4,90), T(257,40,-8,0,5,0)];        % 强弱 ΔD=0.5
mi=mi+1; sc(mi).name='P6'; sc(mi).seed=4;  sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(285,0,-0.5,0,10,0), T(285,1.0,+0.5,0,-4,90), T(257,40,-8,0,5,0)];        % 强弱 ΔD=1
scenes = sc;  RES = zeros(0, 15);  nOk = 0;  nFail = 0;  nextIdx = 1;
save('C:\Users\NJ\AppData\Local\Temp\probe95\rerunD_state.mat', 'scenes', 'RES', 'nOk', 'nFail', 'nextIdx');
fprintf('PD INIT %d scenes\n', numel(scenes));
