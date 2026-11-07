% D 系触发窗口探针(2026-11-08 新 RF 态重设裕量): 分离度×成员电平网格 + 跨种子
% 几何对照 docsrc/rerunInit.m(Nov-7 精修表): D01=±0.5°/10+10, D03=±0.4°/12+4, D06=±0.5°/10+2
% 稳定模式观察: 弱成员 ≤0dBsm(D02/-4, D05/0, D07/-4)两时代都稳定"对不分离";
% 双稳带: 等强 10/10@1.0°(D01) 与 弱成员 2/4dBsm(D06/D03) 在新 RF(线功率+2~3dB)下翻态。
T = @(gR,gD,az,el,rcs,ph) struct('gate_R',gR,'gate_D',gD,'azi',az,'ele',el, ...
    'rcsDb',rcs,'rcsPhDeg',ph,'vTan',0,'isCorner',0);
C = struct('gate_R',5,'gate_D',0,'azi',0,'ele',0,'rcsDb',-5,'rcsPhDeg',0,'vTan',0,'isCorner',1);
sc = struct('name',{},'seed',{},'tg',{});  mi = 0;
% ---- D01 型(等强对 1.0°, 双稳带) ----
mi=mi+1; sc(mi).name='W01'; sc(mi).seed=4;  sc(mi).tg=[C, T(285,0,-0.5,0,10,0),  T(285,0,0.5,0,10,90),  T(257,40,-8,0,5,0)];  % =D01 基线复测
mi=mi+1; sc(mi).name='W02'; sc(mi).seed=4;  sc(mi).tg=[C, T(285,0,-0.8,0,10,0),  T(285,0,0.8,0,10,90),  T(257,40,-8,0,5,0)];  % 1.6° 10/10
mi=mi+1; sc(mi).name='W03'; sc(mi).seed=4;  sc(mi).tg=[C, T(285,0,-0.5,0,12,0),  T(285,0,0.5,0,12,90),  T(257,40,-8,0,5,0)];  % 1.0° 12/12
mi=mi+1; sc(mi).name='W04'; sc(mi).seed=4;  sc(mi).tg=[C, T(285,0,-0.8,0,12,0),  T(285,0,0.8,0,12,90),  T(257,40,-8,0,5,0)];  % 1.6° 12/12
% ---- D03 型(强/弱 0.8°, 双稳带) ----
mi=mi+1; sc(mi).name='W05'; sc(mi).seed=4;  sc(mi).tg=[C, T(257,30,-0.4,0,12,0), T(257,30,0.4,0,4,45),  T(357,-10,5,0,10,0)];  % =D03 基线复测
mi=mi+1; sc(mi).name='W06'; sc(mi).seed=4;  sc(mi).tg=[C, T(257,30,-0.4,0,12,0), T(257,30,0.4,0,0,45),  T(357,-10,5,0,10,0)];  % 弱 0dBsm(明确不可分)
mi=mi+1; sc(mi).name='W07'; sc(mi).seed=4;  sc(mi).tg=[C, T(257,30,-0.4,0,12,0), T(257,30,0.4,0,8,45),  T(357,-10,5,0,10,0)];  % 弱 8dBsm(明确可分)
% ---- D06 型(强/弱 1.0°, 双稳带) ----
mi=mi+1; sc(mi).name='W08'; sc(mi).seed=4;  sc(mi).tg=[C, T(314,-60,-0.5,0,10,60), T(314,-60,0.5,0,2,240), T(200,20,12,0,10,0)]; % =D06 基线复测
mi=mi+1; sc(mi).name='W09'; sc(mi).seed=4;  sc(mi).tg=[C, T(314,-60,-0.5,0,10,60), T(314,-60,0.5,0,0,240), T(200,20,12,0,10,0)]; % 弱 0dBsm(明确不可分)
% ---- 跨种子稳定性 ----
mi=mi+1; sc(mi).name='W10'; sc(mi).seed=11; sc(mi).tg=[C, T(285,0,-0.8,0,10,0),  T(285,0,0.8,0,10,90),  T(257,40,-8,0,5,0)];   % W02 @seed11
mi=mi+1; sc(mi).name='W11'; sc(mi).seed=11; sc(mi).tg=[C, T(257,30,-0.4,0,12,0), T(257,30,0.4,0,8,45),  T(357,-10,5,0,10,0)]; % W07 @seed11
scenes = sc;  RES = zeros(0, 15);  nOk = 0;  nFail = 0;  nextIdx = 1;
save('C:\Users\NJ\AppData\Local\Temp\probe95\rerunD_state.mat', 'scenes', 'RES', 'nOk', 'nFail', 'nextIdx');
fprintf('DWIN INIT %d scenes\n', numel(scenes));
