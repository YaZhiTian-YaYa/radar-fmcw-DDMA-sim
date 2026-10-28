function [targets, truth] = genPairScene(seed, nOmp, nDml, nSingle, weakFrac, snrOffDb)
%GENPAIRSCENE 生成稀疏多目标随机场景(批量案例用, 与 main_sim 的 B/cr 一致)
%   [targets, truth] = genPairScene(seed, nOmp, nDml, nSingle, weakFrac, snrOffDb)
%   场景风格(十字路口): 前方车流(Y 向行驶) + 横穿车(X 向, 单目标里 30%) + 角反。
%   【多而不密】: 任意两个 RD 单元(距离bin, 谱线bin)至少满足其一 ——
%     距离差 ≥ 16 bin(≈9 m) 或 谱线差 ≥ 30 bin(速度差 ≈1.9 m/s) ——
%     目标总数可以多, 但 RD 平面上彼此散开, 不出现"堵车式"扎堆。
%   结构: 角反(目标1) + 第1对 OMP 型(2/3, OMP 种子) + 第2对 DML 型(4/5, DML 种子)
%         + 其余 OMP/DML 型对 + 单目标(考核误分裂)。每对两成员同距同速bin(独占该单元)。
%   案例维度随机化:
%     OMP 型: 角度差 5~28°, RCS 差 {0,5,10,15} dB;  DML 型: 角度差 1.5~6°, RCS 差 ±2.5 dB
%     距离 60~440 bin 随机(量程/LPF 带内); 30% 对为 n=1 同速bin(gate_D+512, 低速组合;
%       真速差 = 512×Vres, 随 Tc 联动 —— 帧内走动恒 ≈1 m, 与 Tc 无关)
%     全部对约束在 ROI 内(|R·sin az|≤25 m 留 1° 裕量)
%   单目标: 随机角度/距离/速度; 运动方向 dir: 70% Y 向(沿路行驶), 30% X 向(横穿, |az|≥8°),
%     vTan 按 (gate_D, azi, ele) 矢量分解自动算出(切向不产生多普勒, 物理自洽)
%   weakFrac(默认 0): 对的弱成员再压 5~15 dB(概率 weakFrac), 单目标压 5~10 dB
%     (概率 weakFrac/2) —— 制造低信噪比案例考核 OMP/DML 下限(允许漏检)
%   snrOffDb(默认 0): 全场景 RCS 压低量(dB, 角反除外) —— 整场景低信噪比案例:
%     0=常规可检测基准; 8~15 dB 逐级压低, 大量目标接近/低于检测门限(考核下限)
%   truth: struct('cell',[gate_R, dopC], 'type'(1=OMP 2=DML 0=单), 'az1','az2','rcs1','rcs2')

    rng(seed);
    if nargin < 5, weakFrac = 0; end
    if nargin < 6, snrOffDb = 0; end
    B_ = 256.844e6;  cr_ = 1;  dopBinNum_ = 512;      % 与 main_sim 一致
    stepDir_ = 1;                                      % 与 main_sim 的 stepDir 联动(下跳频时同步改 −1)
    Tc_ = 37e-6;                                       % 与 main_sim 的 Tc 一致(改 Tc 时同步这里)
    fC_ = 76.2e9 + 2*B_/2;                             % 全帧中心频率(与 main_sim 的 fCenter 同式)
    Vres_ = 3e8/(2*fC_*Tc_*dopBinNum_);                % 速度分辨率(与 dims.Vres 同式, vTan 分解用)
    DeltaR_ = 3e8/(2*B_);
    rcsBase = @(gr) 8 + 40*log10(gr/64) - snrOffDb;   % 按距离补偿的基准(dBsm) − 全场压低量
    dopOf   = @(gr, gd) mod(gd + stepDir_*cr_*B_/(dopBinNum_-1) * (2*gr*3e8/(2*B_))/3e8 * dopBinNum_, dopBinNum_);

    % ---- 角反(校准) ----
    targets = struct('gate_R',5,'gate_D',0,'azi',0,'ele',0,'rcsDb',-5,'vTan',0,'isCorner',1);
    truth = repmat(struct('cell',[],'type',[],'az1',[],'az2',[],'rcs1',[],'rcs2',[]), 0, 1);

    % ---- 稀疏化: RD 单元二维最小间隔(距离 16 bin 或 谱线 30 bin 环形差) ----
    minGapR = 16;  minGapD = 30;
    cells = [5, dopOf(5, 0)];                         % 已占用单元(角反先占)
    poolR = 60:4:440;  poolR = poolR(randperm(numel(poolR)));   % 距离候选(洗牌不重复)
    ip = 0;

    % ---- 案例对 ----
    for k = 1:(nOmp + nDml)
        if k == 1,      typeOmp = true;               % 第1对 OMP 型(OMP 过程图种子)
        elseif k == 2,  typeOmp = false;              % 第2对 DML 型(DML 过程图种子)
        else,           typeOmp = (k <= nOmp + 1);    % 其余: 先补齐 OMP 型, 再 DML 型
        end
        n1 = rand < 0.3;                               % 30% 对做 n=1 同速bin(先定, 决定抽速上限)
        [gr, gd] = drawCellInner(200*n1 + 480*~n1);   % 抽满足间隔的单元(n=1 限低速 gd≤200 防走动)
        R = gr * DeltaR_;
        lim = asind(min(0.99, 25/R)) - 1;             % ROI 内最大 |az|(留 1° 裕量)
        if typeOmp, sep = 5 + rand*23;  else, sep = 1.5 + rand*4.5;  end
        if sep/2 > lim, sep = 2*lim*0.8;  end         % 远距时收缩间隔(类型语义不变)
        azC = (2*rand-1) * max(lim - sep/2, 0.5);     % 对中心角随机, 保证两角都在 ROI
        az1 = azC - sep/2;  az2 = azC + sep/2;
        if rand < 0.5, t_ = az1; az1 = az2; az2 = t_; end
        if typeOmp, dr = [0 5 10 15]; dr = dr(randi(4));  else, dr = 5*rand - 2.5;  end
        rcs1 = rcsBase(gr) + 2 + dr/2;   rcs2 = rcs1 - dr;
        if rand < weakFrac, rcs2 = rcs2 - (5 + 10*rand);  end   % 弱成员再压(低信噪比案例)
        ele = (2*rand-1)*8;
        n1 = rand < 0.3;                              % 30% 概率 n=1 同速bin(低速组合)
        gd2 = gd + 512*n1;                            %   (n=1 的低速约束已在抽单元时给出)
        targets(end+1) = struct('gate_R',gr,'gate_D',gd ,'azi',az1,'ele',ele,'rcsDb',rcs1,'vTan',0,'isCorner',0);
        targets(end+1) = struct('gate_R',gr,'gate_D',gd2,'azi',az2,'ele',ele,'rcsDb',rcs2,'vTan',0,'isCorner',0);
        truth(end+1) = struct('cell',[gr, dopOf(gr,gd)], 'type', double(~typeOmp) + 1, ...
                              'az1',az1, 'az2',az2, 'rcs1',rcs1, 'rcs2',rcs2);  % type: 1=OMP型, 2=DML型
    end

    % ---- 单目标(考核误分裂; 70% Y 向行驶 / 30% X 向横穿) ----
    for k = 1:nSingle
        [gr, gd] = drawCellInner();
        R = gr * DeltaR_;
        lim = asind(min(0.99, 25/R)) - 1;
        dirX = rand < 0.3;                            % X 向横穿比例
        if dirX, az = sign(2*rand-1) * (8 + rand*(lim-8));  else, az = (2*rand-1)*lim;  end
        ele = (2*rand-1)*10;
        rcs = rcsBase(gr) + 6*rand - 2;
        if rand < weakFrac/2, rcs = rcs - (5 + 5*rand);  end    % 部分单目标也压(低信噪比)
        vRad = gd * Vres_;                            % 法向速度(m/s; Vres 与主链路同式联动)
        if dirX
            vTan = vRad / (tand(az) * cosd(ele));     % X 向横穿: 切向 = 法向/(tan·cos)
        else
            vTan = -vRad * tand(az) / cosd(ele);      % Y 向行驶: 切向 = −法向·tan/cos
        end
        targets(end+1) = struct('gate_R',gr,'gate_D',gd,'azi',az,'ele',ele,'rcsDb',rcs,'vTan',vTan,'isCorner',0);
        truth(end+1) = struct('cell',[gr, dopOf(gr,gd)], 'type', 0, 'az1',az, 'az2',NaN, ...
                              'rcs1',rcs, 'rcs2',NaN);
    end

    % ---- 局部函数: 抽一个满足二维间隔的 (gr, gd), 并登记占用 ----
    function [gr, gd] = drawCellInner(maxGd)
        gr = 0;  gd = 0;
        if nargin < 1, maxGd = 480;  end
        for try_ = 1:80
            if ip >= numel(poolR), ip = 0;  poolR = poolR(randperm(numel(poolR)));  end %#ok 距离候选循环使用
            ip = ip + 1;
            grTry = poolR(ip);
            gdTry = randi([0 maxGd]);
            dTry = dopOf(grTry, gdTry);
            % 间隔检查必须在此处做(嵌套函数直接读活动工作区的 cells;
            % 匿名函数按创建时刻的值捕获, 看不到后续登记的占用 —— 曾因此失效)
            ok = all( abs(grTry - cells(:,1)) >= minGapR | ...
                      abs(mod(dTry - cells(:,2) + 256, 512) - 256) >= minGapD );
            if ok
                gr = grTry;  gd = gdTry;
                cells(end+1, :) = [gr, dTry];          %#ok<AGROW> 登记占用(嵌套函数直接可见)
                return;
            end
        end
        % 80 次未抽到(极小概率): 放宽接受, 保证函数总能返回
        ip = ip + 1;  if ip > numel(poolR), ip = 1;  end
        gr = poolR(ip);  gd = randi([0 maxGd]);
        cells(end+1, :) = [gr, dopOf(gr, gd)];
    end
end
