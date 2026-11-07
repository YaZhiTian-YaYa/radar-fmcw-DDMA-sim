% 初始化进程隔离批跑状态: 28 场景(I01-I12/O01-O08/D01-D08, 与 docsrc/finalBatch.m 场景表一致)
T = @(gR,gD,az,el,rcs,ph) struct('gate_R',gR,'gate_D',gD,'azi',az,'ele',el, ...
    'rcsDb',rcs,'rcsPhDeg',ph,'vTan',0,'isCorner',0);
C = struct('gate_R',5,'gate_D',0,'azi',0,'ele',0,'rcsDb',-5,'rcsPhDeg',0,'vTan',0,'isCorner',1);
sc = struct('name',{},'seed',{},'tg',{});  mi = 0;
% ---- 批1 独立目标 ----
mi=mi+1; sc(mi).name='I01'; sc(mi).seed=4;  sc(mi).tg=[C, T(214,20,5,0,10,0), T(285,-40,-8,0,4,30), T(357,60,12,0,7,60)];
mi=mi+1; sc(mi).name='I02'; sc(mi).seed=4;  sc(mi).tg=[C, T(414,10,0,0,12,0), T(257,-60,-15,0,5,90), T(171,30,8,0,4,20)];
mi=mi+1; sc(mi).name='I03'; sc(mi).seed=4;  sc(mi).tg=[C, T(357,-20,3,0,15,0), T(414,50,-5,0,10,120), T(285,0,15,0,4,200)];
mi=mi+1; sc(mi).name='I04'; sc(mi).seed=4;  sc(mi).tg=[C, T(143,-10,-12,0,10,0), T(414,20,5,0,10,45), T(400,-50,-3,0,8,75)];
mi=mi+1; sc(mi).name='I05'; sc(mi).seed=4;  sc(mi).tg=[C, T(142.5,15,6,0,10,0), T(200.5,-25,-10,0,4,30), T(270.5,45,14,0,7,150)];
mi=mi+1; sc(mi).name='I06'; sc(mi).seed=11; sc(mi).tg=[C, T(214,20,5,0,10,0), T(285,-40,-8,0,4,30), T(357,60,12,0,7,60)];
mi=mi+1; sc(mi).name='I07'; sc(mi).seed=23; sc(mi).tg=[C, T(414,10,0,0,12,0), T(257,-60,-15,0,5,90), T(171,30,8,0,4,20)];
mi=mi+1; sc(mi).name='I08'; sc(mi).seed=35; sc(mi).tg=[C, T(142.5,15,6,0,10,0), T(200.5,-25,-10,0,4,30), T(270.5,45,14,0,7,150)];
mi=mi+1; sc(mi).name='I09'; sc(mi).seed=4;  sc(mi).tg=[C, T(214,30,70,0,15,0), T(214,-30,-70,0,15,45), T(257,50,40,0,10,90)];
mi=mi+1; sc(mi).name='I10'; sc(mi).seed=4;  sc(mi).tg=[C, T(214,30,0,12,10,0), T(257,-30,0,-10,10,45), T(285,50,5,8,10,90)];
mi=mi+1; sc(mi).name='I11'; sc(mi).seed=4;  sc(mi).tg=[C, T(357,60,-12,0,7,30), T(400,-20,8,0,9,0)];
mi=mi+1; sc(mi).name='I12'; sc(mi).seed=47; sc(mi).tg=[C, T(285,0,0,0,10,0)];
% ---- 批2 OMP 对 ----
mi=mi+1; sc(mi).name='O01'; sc(mi).seed=4;  sc(mi).tg=[C, T(285,40,-5,0,12,45), T(285,40,5,0,2,210), T(400,30,-15,0,8,0)];
mi=mi+1; sc(mi).name='O02'; sc(mi).seed=4;  sc(mi).tg=[C, T(357,0,-3,0,10,0), T(357,0,3,0,10,90), T(214,-20,8,0,5,30)];
mi=mi+1; sc(mi).name='O03'; sc(mi).seed=4;  sc(mi).tg=[C, T(214,-30,-4,0,4,20), T(214,-30,4,0,10,140), T(414,10,0,0,10,0)];
mi=mi+1; sc(mi).name='O04'; sc(mi).seed=4;  sc(mi).tg=[C, T(414,20,-6,0,5,0), T(414,20,6,0,10,60), T(257,-50,10,0,10,90)];
mi=mi+1; sc(mi).name='O05'; sc(mi).seed=4;  sc(mi).tg=[C, T(314,60,-2.5,0,10,0), T(314,60,2.5,0,4,120), T(357,-35,-12,0,7,30)];
mi=mi+1; sc(mi).name='O06'; sc(mi).seed=4;  sc(mi).tg=[C, T(257,-45,-5,0,4,0), T(257,-45,5,0,10,90), T(371,25,6,0,8,45)];
mi=mi+1; sc(mi).name='O07'; sc(mi).seed=11; sc(mi).tg=[C, T(285,40,-5,0,12,45), T(285,40,5,0,2,210), T(400,30,-15,0,8,0)];
mi=mi+1; sc(mi).name='O08'; sc(mi).seed=4;  sc(mi).tg=[C, T(357,0,-4,0,12,0), T(357,0,4,0,0,90), T(214,30,10,0,10,0)];
% ---- 批3 DML 对 ----
mi=mi+1; sc(mi).name='D01'; sc(mi).seed=4;  sc(mi).tg=[C, T(285,0,-0.5,0,10,0), T(285,0,0.5,0,10,90), T(257,40,-8,0,5,0)];
mi=mi+1; sc(mi).name='D02'; sc(mi).seed=4;  sc(mi).tg=[C, T(257,30,-0.4,0,10,0), T(257,30,0.4,0,10,45), T(414,-10,5,0,7,90)];
mi=mi+1; sc(mi).name='D03'; sc(mi).seed=4;  sc(mi).tg=[C, T(285,-20,-0.6,0,10,0), T(285,-20,0.6,0,6,120), T(214,60,-6,0,12,0)];
mi=mi+1; sc(mi).name='D04'; sc(mi).seed=4;  sc(mi).tg=[C, T(414,10,-0.5,0,7,0), T(414,10,0.5,0,7,90), T(257,-30,10,0,10,30)];
mi=mi+1; sc(mi).name='D05'; sc(mi).seed=4;  sc(mi).tg=[C, T(214,40,-0.6,0,12,0), T(214,40,0.6,0,0,60), T(357,-60,-14,0,7,0)];
mi=mi+1; sc(mi).name='D06'; sc(mi).seed=4;  sc(mi).tg=[C, T(314,-60,-0.4,0,6,0), T(314,-60,0.4,0,10,90), T(400,20,12,0,10,45)];
mi=mi+1; sc(mi).name='D07'; sc(mi).seed=11; sc(mi).tg=[C, T(285,0,-0.5,0,10,0), T(285,0,0.5,0,10,90), T(257,40,-8,0,5,0)];
mi=mi+1; sc(mi).name='D08'; sc(mi).seed=4;  sc(mi).tg=[C, T(171,20,-0.5,0,10,0), T(171,20,0.5,0,6,30)];
scenes = sc;  RES = zeros(0, 15);  nOk = 0;  nFail = 0;  nextIdx = 1;
save('C:\Users\NJ\AppData\Local\Temp\probe95\rerunD_state.mat', 'scenes', 'RES', 'nOk', 'nFail', 'nextIdx');
fprintf('RERUN INIT ALL %d scenes\n', numel(scenes));
