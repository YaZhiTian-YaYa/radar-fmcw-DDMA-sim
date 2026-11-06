% 大批量验证场景(66): I 批 30(独立, 多几何多种子) / O 批 20(OMP 对, sep3~12) / D 批 16(等强多分离度多种子+强弱融合)
T = @(gR,gD,az,el,rcs,ph) struct('gate_R',gR,'gate_D',gD,'azi',az,'ele',el, ...
    'rcsDb',rcs,'rcsPhDeg',ph,'vTan',0,'isCorner',0);
C = struct('gate_R',5,'gate_D',0,'azi',0,'ele',0,'rcsDb',-5,'rcsPhDeg',0,'vTan',0,'isCorner',1);
sc = struct('name',{},'seed',{},'tg',{},'cfarScale',{},'mc',{},'rob',{});  mi = 0;
% ================= I 批 30 =================
mi=mi+1; sc(mi).name='I01'; sc(mi).seed=4;  sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(214,20,5,0,10,0), T(285,-40,-8,0,4,30), T(357,60,12,0,7,60)];
mi=mi+1; sc(mi).name='I02'; sc(mi).seed=4;  sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(414,10,0,0,12,0), T(257,-60,-15,0,5,90), T(171,30,8,0,4,20)];
mi=mi+1; sc(mi).name='I03'; sc(mi).seed=4;  sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(357,-20,3,0,15,0), T(414,50,-5,0,10,120), T(285,0,15,0,4,200)];
mi=mi+1; sc(mi).name='I04'; sc(mi).seed=4;  sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(143,-10,-12,0,10,0), T(414,20,5,0,10,45), T(400,-50,-3,0,8,75)];
mi=mi+1; sc(mi).name='I05'; sc(mi).seed=4;  sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(142.5,15,6,0,10,0), T(200.5,-25,-10,0,4,30), T(270.5,45,14,0,7,150)];
mi=mi+1; sc(mi).name='I06'; sc(mi).seed=4;  sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(314,35,25,0,10,0), T(257,-45,-35,0,7,60), T(371,10,0,0,10,120)];
mi=mi+1; sc(mi).name='I07'; sc(mi).seed=4;  sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(200,25,60,0,12,0), T(200,-25,-60,0,10,45), T(257,50,45,0,7,90)];
mi=mi+1; sc(mi).name='I08'; sc(mi).seed=4;  sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(220,30,70,0,15,0), T(220,-30,-70,0,12,45), T(250,45,40,0,10,90)];
mi=mi+1; sc(mi).name='I09'; sc(mi).seed=4;  sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(214,30,0,12,10,0), T(257,-30,0,-10,10,45), T(285,50,5,8,10,90)];
mi=mi+1; sc(mi).name='I10'; sc(mi).seed=4;  sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(285,15,0,15,10,0), T(314,-25,10,-15,7,60), T(357,40,-12,10,7,120)];
mi=mi+1; sc(mi).name='I11'; sc(mi).seed=4;  sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(171,20,-3,0,7,0), T(229,40,9,0,4,30), T(400,-30,-6,0,12,60)];
mi=mi+1; sc(mi).name='I12'; sc(mi).seed=4;  sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(285,0,0,0,10,0)];
mi=mi+1; sc(mi).name='I13'; sc(mi).seed=4;  sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(245,35,18,0,10,0), T(335,-15,-22,0,7,45), T(429,25,30,0,10,90)];
mi=mi+1; sc(mi).name='I14'; sc(mi).seed=4;  sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(190,-50,45,0,7,0), T(300,55,-42,0,10,30), T(380,-5,38,0,7,60)];
mi=mi+1; sc(mi).name='I15'; sc(mi).seed=4;  sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(205.5,10,-8,0,10,0), T(260.5,-30,12,0,4,45), T(315.5,50,-16,0,7,90)];
mi=mi+1; sc(mi).name='I16'; sc(mi).seed=4;  sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(155,45,55,0,10,0), T(210,-40,-58,0,10,30), T(265,20,62,0,7,60)];
mi=mi+1; sc(mi).name='I17'; sc(mi).seed=4;  sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(230,5,2,0,15,0), T(340,-10,-4,0,15,45), T(425,40,6,0,10,90)];
mi=mi+1; sc(mi).name='I18'; sc(mi).seed=4;  sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(175,-25,-28,0,12,0), T(280,30,35,0,4,30), T(420,-45,-55,0,10,60)];
mi=mi+1; sc(mi).name='I19'; sc(mi).seed=11; sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(214,20,5,0,10,0), T(285,-40,-8,0,4,30), T(357,60,12,0,7,60)];
mi=mi+1; sc(mi).name='I20'; sc(mi).seed=11; sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(414,10,0,0,12,0), T(257,-60,-15,0,5,90), T(171,30,8,0,4,20)];
mi=mi+1; sc(mi).name='I21'; sc(mi).seed=11; sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(142.5,15,6,0,10,0), T(200.5,-25,-10,0,4,30), T(270.5,45,14,0,7,150)];
mi=mi+1; sc(mi).name='I22'; sc(mi).seed=11; sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(220,30,70,0,15,0), T(220,-30,-70,0,12,45), T(250,45,40,0,10,90)];
mi=mi+1; sc(mi).name='I23'; sc(mi).seed=11; sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(214,30,0,12,10,0), T(257,-30,0,-10,10,45), T(285,50,5,8,10,90)];
mi=mi+1; sc(mi).name='I24'; sc(mi).seed=11; sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(285,0,0,0,10,0)];
mi=mi+1; sc(mi).name='I25'; sc(mi).seed=23; sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(214,20,5,0,10,0), T(285,-40,-8,0,4,30), T(357,60,12,0,7,60)];
mi=mi+1; sc(mi).name='I26'; sc(mi).seed=23; sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(143,-10,-12,0,10,0), T(414,20,5,0,10,45), T(400,-50,-3,0,8,75)];
mi=mi+1; sc(mi).name='I27'; sc(mi).seed=23; sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(314,35,25,0,10,0), T(257,-45,-35,0,7,60), T(371,10,0,0,10,120)];
mi=mi+1; sc(mi).name='I28'; sc(mi).seed=23; sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(200,25,60,0,12,0), T(200,-25,-60,0,10,45), T(257,50,45,0,7,90)];
mi=mi+1; sc(mi).name='I29'; sc(mi).seed=23; sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(285,15,0,15,10,0), T(314,-25,10,-15,7,60), T(357,40,-12,10,7,120)];
mi=mi+1; sc(mi).name='I30'; sc(mi).seed=23; sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(245,35,18,0,10,0), T(335,-15,-22,0,7,45), T(429,25,30,0,10,90)];
% ================= O 批 20 =================
mi=mi+1; sc(mi).name='O01'; sc(mi).seed=4;  sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(285,40,-1.5,0,10,0), T(285,40,1.5,0,6,90), T(257,-20,8,0,10,45)];         % sep3
mi=mi+1; sc(mi).name='O02'; sc(mi).seed=11; sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(214,-30,-2.0,0,8,45), T(214,-30,2.0,0,8,135), T(357,35,-10,0,7,0)];       % sep4
mi=mi+1; sc(mi).name='O03'; sc(mi).seed=4;  sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(357,0,-2.0,0,12,0), T(357,0,2.0,0,4,90), T(214,-20,8,0,10,0)];            % sep4 弱二
mi=mi+1; sc(mi).name='O04'; sc(mi).seed=23; sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(257,20,-2.5,0,10,0), T(257,20,2.5,0,10,60), T(314,45,-12,0,7,30)];        % sep5 等强
mi=mi+1; sc(mi).name='O05'; sc(mi).seed=4;  sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(314,60,-2.5,0,10,0), T(314,60,2.5,0,4,120), T(357,-35,-12,0,7,30)];      % sep5
mi=mi+1; sc(mi).name='O06'; sc(mi).seed=11; sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(285,0,-2.5,0,7,0), T(285,0,2.5,0,7,90), T(414,20,10,0,10,45)];            % sep5 中距
mi=mi+1; sc(mi).name='O07'; sc(mi).seed=4;  sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(414,10,-3,0,12,0), T(414,10,3,0,0,90), T(257,-50,10,0,10,30)];            % sep6 远距极弱二
mi=mi+1; sc(mi).name='O08'; sc(mi).seed=4;  sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(200,-40,-3,0,10,0), T(200,-40,3,0,6,120), T(285,30,-8,0,10,60)];          % sep6 近距
mi=mi+1; sc(mi).name='O09'; sc(mi).seed=4;  sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(285,40,-4,0,10,0), T(285,40,4,0,2,210), T(400,30,-15,0,8,0)];             % sep8
mi=mi+1; sc(mi).name='O10'; sc(mi).seed=23; sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(357,-20,-4,0,12,0), T(357,-20,4,0,12,90), T(214,30,10,0,10,0)];           % sep8 等强
mi=mi+1; sc(mi).name='O11'; sc(mi).seed=4;  sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(214,30,-5,0,8,0), T(214,30,5,0,0,90), T(335,-15,-14,0,7,30)];            % sep10
mi=mi+1; sc(mi).name='O12'; sc(mi).seed=4;  sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(371,-60,-5,0,10,0), T(371,-60,5,0,10,60), T(257,40,12,0,7,0)];            % sep10 宽角
mi=mi+1; sc(mi).name='O13'; sc(mi).seed=4;  sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(414,20,-6,0,5,0), T(414,20,6,0,10,60), T(257,-50,10,0,10,90)];            % sep12 弱母强二
mi=mi+1; sc(mi).name='O14'; sc(mi).seed=11; sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(257,50,-6,0,12,0), T(257,50,6,0,4,150), T(314,-35,-18,0,10,30)];          % sep12
mi=mi+1; sc(mi).name='O15'; sc(mi).seed=4;  sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(285,0,-1.25,0,10,0), T(285,0,1.25,0,10,90), T(357,40,-8,0,7,0)];          % sep2.5 OMP下边界
mi=mi+1; sc(mi).name='O16'; sc(mi).seed=4;  sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(335,15,-1.25,0,10,30), T(335,15,1.25,0,4,120), T(214,-30,9,0,10,0)];      % sep2.5
mi=mi+1; sc(mi).name='O17'; sc(mi).seed=11; sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(285,40,-5,0,12,45), T(285,40,5,0,2,210), T(400,30,-15,0,8,0)];            % sep10 大差
mi=mi+1; sc(mi).name='O18'; sc(mi).seed=23; sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(357,0,-3,0,10,0), T(357,0,3,0,10,90), T(214,-20,8,0,5,30)];               % sep6 等强
mi=mi+1; sc(mi).name='O19'; sc(mi).seed=4;  sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(229,25,-4.5,0,10,0), T(229,25,4.5,0,7,150), T(300,-40,-20,0,10,60)];      % sep9
mi=mi+1; sc(mi).name='O20'; sc(mi).seed=11; sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(300,55,-7,0,10,0), T(300,55,7,0,6,240), T(245,-20,15,0,10,30)];           % sep14 宽角
% ================= D 批 16 =================
mi=mi+1; sc(mi).name='D01'; sc(mi).seed=4;  sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(285,0,-0.5,0,10,0), T(285,0,0.5,0,10,90), T(257,40,-8,0,5,0)];            % sep1.0 等强
mi=mi+1; sc(mi).name='D02'; sc(mi).seed=11; sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(285,0,-0.5,0,10,0), T(285,0,0.5,0,10,90), T(257,40,-8,0,5,0)];
mi=mi+1; sc(mi).name='D03'; sc(mi).seed=35; sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(285,0,-0.5,0,10,0), T(285,0,0.5,0,10,90), T(257,40,-8,0,5,0)];
mi=mi+1; sc(mi).name='D04'; sc(mi).seed=4;  sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(285,0,-0.45,0,10,0), T(285,0,0.45,0,10,120), T(257,40,-8,0,5,0)];         % sep0.9
mi=mi+1; sc(mi).name='D05'; sc(mi).seed=11; sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(285,0,-0.45,0,10,0), T(285,0,0.45,0,10,120), T(257,40,-8,0,5,0)];
mi=mi+1; sc(mi).name='D06'; sc(mi).seed=35; sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(285,0,-0.45,0,10,0), T(285,0,0.45,0,10,120), T(257,40,-8,0,5,0)];
mi=mi+1; sc(mi).name='D07'; sc(mi).seed=4;  sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(285,0,-0.55,0,10,0), T(285,0,0.55,0,10,45), T(257,40,-8,0,5,0)];          % sep1.1
mi=mi+1; sc(mi).name='D08'; sc(mi).seed=11; sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(285,0,-0.55,0,10,0), T(285,0,0.55,0,10,45), T(257,40,-8,0,5,0)];
mi=mi+1; sc(mi).name='D09'; sc(mi).seed=35; sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(285,0,-0.55,0,10,0), T(285,0,0.55,0,10,45), T(257,40,-8,0,5,0)];
mi=mi+1; sc(mi).name='D10'; sc(mi).seed=4;  sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(214,30,-0.5,0,10,0), T(214,30,0.5,0,10,90), T(357,-40,-6,0,10,0)];        % sep1.0 近距
mi=mi+1; sc(mi).name='D11'; sc(mi).seed=4;  sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(357,-20,-0.55,0,10,0), T(357,-20,0.55,0,10,90), T(214,50,8,0,12,0)];      % sep1.1 远距
mi=mi+1; sc(mi).name='D12'; sc(mi).seed=4;  sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(314,45,-0.45,0,10,0), T(314,45,0.45,0,10,60), T(257,-35,-12,0,7,30)];    % sep0.9 宽角簇
mi=mi+1; sc(mi).name='D13'; sc(mi).seed=4;  sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(285,0,-0.5,0,10,0), T(285,0,0.5,0,0,90), T(257,40,-8,0,5,0)];             % 强弱 10/0 → 预期融合
mi=mi+1; sc(mi).name='D14'; sc(mi).seed=11; sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(285,0,-0.5,0,10,0), T(285,0,0.5,0,-4,90), T(257,40,-8,0,5,0)];            % 强弱 10/-4 → 融合
mi=mi+1; sc(mi).name='D15'; sc(mi).seed=23; sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(214,40,-0.6,0,12,0), T(214,40,0.6,0,0,120), T(285,-30,10,0,10,0)];       % 强弱 12/0 → 融合
mi=mi+1; sc(mi).name='D16'; sc(mi).seed=4;  sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).rob=false;
sc(mi).tg=[C, T(357,0,-0.5,0,10,0), T(357,0,0.5,0,4,90), T(214,60,-6,0,12,0)];           % 10/4 中间差 → 看走向
scenes = sc;  RES = zeros(0, 15);  nOk = 0;  nFail = 0;  nextIdx = 1;
save('C:\Users\NJ\AppData\Local\Temp\probe95\rerunD_state.mat', 'scenes', 'RES', 'nOk', 'nFail', 'nextIdx');
fprintf('BIG INIT %d scenes\n', numel(scenes));
