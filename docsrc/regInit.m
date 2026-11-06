% 遗留① DOA 二维耦合偏差场测量: az×el 网格单目标(每场景4个不同RD目标, 省场次)
% az∈{0,15,30,50} × el∈{0,5,10,15} + az=70 列 + 负俯仰抽查; RCS 10dBsm(强目标, 测系统性)
T = @(gR,gD,az,el,rcs,ph) struct('gate_R',gR,'gate_D',gD,'azi',az,'ele',el, ...
    'rcsDb',rcs,'rcsPhDeg',ph,'vTan',0,'isCorner',0);
C = struct('gate_R',5,'gate_D',0,'azi',0,'ele',0,'rcsDb',-5,'rcsPhDeg',0,'vTan',0,'isCorner',1);
sc = struct('name',{},'seed',{},'tg',{},'cfarScale',{},'mc',{},'azDec',{},'win',{});  mi = 0;
mi=mi+1; sc(mi).name='I01'; sc(mi).seed=4; sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).azDec=true; sc(mi).win=[];
sc(mi).tg=[C, T(214,20,5,0,10,0), T(285,-40,-8,0,4,30), T(357,60,12,0,7,60)];
mi=mi+1; sc(mi).name='I09'; sc(mi).seed=4; sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).azDec=true; sc(mi).win=[];
sc(mi).tg=[C, T(214,30,70,0,15,0), T(214,-30,-70,0,15,45), T(257,50,40,0,10,90)];
mi=mi+1; sc(mi).name='I10'; sc(mi).seed=4; sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).azDec=true; sc(mi).win=[];
sc(mi).tg=[C, T(214,30,0,12,10,0), T(257,-30,0,-10,10,45), T(285,50,5,8,10,90)];
scenes = sc;  RES = zeros(0, 15);  nOk = 0;  nFail = 0;  nextIdx = 1;
save('C:\Users\NJ\AppData\Local\Temp\probe95\rerunD_state.mat', 'scenes', 'RES', 'nOk', 'nFail', 'nextIdx');
fprintf('BIAS INIT %d scenes\n', numel(scenes));
