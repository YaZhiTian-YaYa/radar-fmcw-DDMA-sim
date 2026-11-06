% 遗留① DOA 二维耦合偏差场测量: az×el 网格单目标(每场景4个不同RD目标, 省场次)
% az∈{0,15,30,50} × el∈{0,5,10,15} + az=70 列 + 负俯仰抽查; RCS 10dBsm(强目标, 测系统性)
T = @(gR,gD,az,el,rcs,ph) struct('gate_R',gR,'gate_D',gD,'azi',az,'ele',el, ...
    'rcsDb',rcs,'rcsPhDeg',ph,'vTan',0,'isCorner',0);
C = struct('gate_R',5,'gate_D',0,'azi',0,'ele',0,'rcsDb',-5,'rcsPhDeg',0,'vTan',0,'isCorner',1);
sc = struct('name',{},'seed',{},'tg',{},'cfarScale',{},'mc',{});  mi = 0;
mi=mi+1; sc(mi).name='B0';  sc(mi).seed=4; sc(mi).cfarScale=[]; sc(mi).mc=false;
sc(mi).tg=[C, T(143,10,0,0,10,0), T(214,-20,15,0,10,45), T(285,40,30,0,10,90), T(357,-50,50,0,10,135)];
mi=mi+1; sc(mi).name='B5';  sc(mi).seed=4; sc(mi).cfarScale=[]; sc(mi).mc=false;
sc(mi).tg=[C, T(143,10,0,5,10,0), T(214,-20,15,5,10,45), T(285,40,30,5,10,90), T(357,-50,50,5,10,135)];
mi=mi+1; sc(mi).name='B10'; sc(mi).seed=4; sc(mi).cfarScale=[]; sc(mi).mc=false;
sc(mi).tg=[C, T(143,10,0,10,10,0), T(214,-20,15,10,10,45), T(285,40,30,10,10,90), T(357,-50,50,10,10,135)];
mi=mi+1; sc(mi).name='B15'; sc(mi).seed=4; sc(mi).cfarScale=[]; sc(mi).mc=false;
sc(mi).tg=[C, T(143,10,0,15,10,0), T(214,-20,15,15,10,45), T(285,40,30,15,10,90), T(357,-50,50,15,10,135)];
mi=mi+1; sc(mi).name='Bm10'; sc(mi).seed=4; sc(mi).cfarScale=[]; sc(mi).mc=false;
sc(mi).tg=[C, T(143,10,0,-10,10,0), T(214,-20,15,-10,10,45), T(285,40,30,-10,10,90), T(357,-50,50,-10,10,135)];
mi=mi+1; sc(mi).name='B70'; sc(mi).seed=4; sc(mi).cfarScale=[]; sc(mi).mc=false;
sc(mi).tg=[C, T(143,10,70,0,10,0), T(214,-20,70,5,10,45), T(285,40,70,10,10,90), T(357,-50,70,15,10,135)];
mi=mi+1; sc(mi).name='Bm15'; sc(mi).seed=4; sc(mi).cfarScale=[]; sc(mi).mc=false;
sc(mi).tg=[C, T(143,10,0,-15,10,0), T(214,-20,15,-15,10,45), T(285,40,30,-15,10,90), T(357,-50,50,-15,10,135)];
scenes = sc;  RES = zeros(0, 15);  nOk = 0;  nFail = 0;  nextIdx = 1;
save('C:\Users\NJ\AppData\Local\Temp\probe95\rerunD_state.mat', 'scenes', 'RES', 'nOk', 'nFail', 'nextIdx');
fprintf('BIAS INIT %d scenes\n', numel(scenes));
