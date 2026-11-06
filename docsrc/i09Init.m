% 遗留④ SVA 档 winStat 重标定: 多目标 2 种子(8 样本), 从审计 差=仿真−解析 中位数刷新
T = @(gR,gD,az,el,rcs,ph) struct('gate_R',gR,'gate_D',gD,'azi',az,'ele',el, ...
    'rcsDb',rcs,'rcsPhDeg',ph,'vTan',0,'isCorner',0);
C = struct('gate_R',5,'gate_D',0,'azi',0,'ele',0,'rcsDb',-5,'rcsPhDeg',0,'vTan',0,'isCorner',1);
sc = struct('name',{},'seed',{},'tg',{},'cfarScale',{},'mc',{},'azDec',{},'win',{});  mi = 0;
mi=mi+1; sc(mi).name='I09'; sc(mi).seed=4; sc(mi).cfarScale=[]; sc(mi).mc=false; sc(mi).azDec=true; sc(mi).win=[];
sc(mi).tg=[C, T(214,30,70,0,15,0), T(214,-30,-70,0,15,45), T(257,50,40,0,10,90)];
scenes = sc;  RES = zeros(0, 15);  nOk = 0;  nFail = 0;  nextIdx = 1;
save('C:\Users\NJ\AppData\Local\Temp\probe95\rerunD_state.mat', 'scenes', 'RES', 'nOk', 'nFail', 'nextIdx');
fprintf('SVA INIT %d scenes\n', numel(scenes));
