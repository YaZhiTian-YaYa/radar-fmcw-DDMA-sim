% 初始化进程隔离批跑状态(9 场景: O04 + D01~D08)
% 2026-11-08 按新 RF 态(线功率+2~3dB)重设 D 系裕量: 只用 probeDwin 探针证过的稳定点。
%   可分×3: 1.2°@214m 12+12(+2dB 裕量) / 1.6°@285m 10+10(跨种子稳, 角色2 稳定低读
%           5.7dB = 分裂边界态如实保留) / 1.2°@214m 10+10(D04 原版, 两时代稳);
%           双稳带(1.0°@285m 等强, W01/W03 坍缩↔分裂翻态)弃用。
%   对不分离×3: 0.8° 12+0(W06 原点) / 0.5° 12+0(D05 原版) / 1.0° 10+0(W09 原点)
%           —— 弱成员≤0dBsm 稳定读强线(单点读数=强成员, splitTag=0)。
T = @(gR,gD,az,el,rcs,ph) struct('gate_R',gR,'gate_D',gD,'azi',az,'ele',el, ...
    'rcsDb',rcs,'rcsPhDeg',ph,'vTan',0,'isCorner',0);
C = struct('gate_R',5,'gate_D',0,'azi',0,'ele',0,'rcsDb',-5,'rcsPhDeg',0,'vTan',0,'isCorner',1);
sc = struct('name',{},'seed',{},'tg',{});  mi = 0;
mi=mi+1; sc(mi).name='O04'; sc(mi).seed=4;  sc(mi).tg=[C, T(371,20,-6,0,5,0), T(371,20,6,0,10,60), T(257,-50,10,0,10,90)];
mi=mi+1; sc(mi).name='D01'; sc(mi).seed=4;  sc(mi).tg=[C, T(214,-20,-0.6,0,12,30), T(214,-20,0.6,0,12,150), T(357,60,-6,0,7,0)];
mi=mi+1; sc(mi).name='D02'; sc(mi).seed=4;  sc(mi).tg=[C, T(285,0,-0.8,0,10,0), T(285,0,0.8,0,10,90), T(257,40,-8,0,5,0)];
mi=mi+1; sc(mi).name='D03'; sc(mi).seed=4;  sc(mi).tg=[C, T(257,30,-0.4,0,12,0), T(257,30,0.4,0,0,45), T(357,-10,5,0,10,0)];
mi=mi+1; sc(mi).name='D04'; sc(mi).seed=4;  sc(mi).tg=[C, T(214,-20,-0.6,0,10,30), T(214,-20,0.6,0,10,150), T(357,60,-6,0,12,0)];
mi=mi+1; sc(mi).name='D05'; sc(mi).seed=4;  sc(mi).tg=[C, T(214,40,-0.5,0,12,0), T(214,40,0.5,0,0,120), T(285,-30,10,0,10,0)];
mi=mi+1; sc(mi).name='D06'; sc(mi).seed=4;  sc(mi).tg=[C, T(314,-60,-0.5,0,10,60), T(314,-60,0.5,0,0,240), T(200,20,12,0,10,0)];
mi=mi+1; sc(mi).name='D07'; sc(mi).seed=11; sc(mi).tg=[C, T(214,-20,-0.6,0,12,30), T(214,-20,0.6,0,12,150), T(357,60,-6,0,7,0)];
mi=mi+1; sc(mi).name='D08'; sc(mi).seed=4;  sc(mi).tg=[C, T(171,20,-0.4,0,10,0), T(171,20,0.4,0,10,30)];
scenes = sc;  RES = zeros(0, 15);  nOk = 0;  nFail = 0;  nextIdx = 1;
save('C:\Users\NJ\AppData\Local\Temp\probe95\rerunD_state.mat', 'scenes', 'RES', 'nOk', 'nFail', 'nextIdx');
fprintf('RERUN INIT %d scenes\n', numel(scenes));
