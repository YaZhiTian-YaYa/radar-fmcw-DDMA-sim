function [tgtOut, detCfg, dmlCfg] = dmlRefine_pipe(tgtOut, detCfg, AntCaliData, rowLut, dims, thrDmlDb, seedR, seedD, drawFigs)
%DMLREFINE_PIPE 对全部点云(不含 OMP 第二目标)逐点做 DML 角度超分辨(最多两目标)
%   [tgtOut, detCfg, dmlCfg] = dmlRefine_pipe(tgtOut, detCfg, AntCaliData, rowLut, dims, thrDmlDb, seedR, seedD)
%   输入:
%     tgtOut      OMP 模块输出的点云列表(OMP 第二目标也在内 —— 本模块跳过它不做 DML)
%     detCfg      doaChanIn: 各点 64 通道快照(每行一点, 与 tgtOut 原始行序一致)
%     AntCaliData 虚拟导向矢量表(1° 间隔, 工程采集格式)
%     rowLut      通道俯仰行查找表(取最长行通道做方位维)
%     dims        公共维度
%     thrDmlDb    【第二目标判定门限】残差改善超过此值(dB)才输出第二目标
%     seedR/seedD 指定做 DML 过程图的目标的 距离bin/多普勒bin(1 起, 容差 ±1 bin;
%                 注意步进频耦合会使实际 D 序号偏离放置值, 按 RD 图实际位置填)
%   处理(逐点依次, 算法逻辑与已验证的 dmlAngle 完全一致):
%     0) 跳过 OMP 第二目标; 门限检查(|DBF方位|≤20° 且 |x|≤25 m, 旋钮见下)超出跳过;
%     1) K=1 粗搜: 1° 表上最佳单角度(强目标);
%     2) 残差匹配曲线取候选峰(限于 DBF 主峰 ±dmlSpanDeg 内, 破 39° 边缘假峰)
%        → 各起点做交替坐标细分(±1° 内 dmlStepDeg 步进 ×3 轮, 表格插值原子);
%     3) 取最终代价最小的起点为解;
%     4) 判定: 残差改善 imprDb = 10log10(J1/J2) > thrDmlDb 且 0.2° < 间隔 ≤ sepLimDeg(6°, 与 OMP 衔接)
%        → 通过: 第一目标方位替换为 DML 第一角度, 第二目标追加点云(继承俯仰等);
%        不通过: 该点保留 DML 之前的 DOA 结果(方位不动, 不追加)。
%   输出:
%     tgtOut  更新后的点云列表(全部原点新增 isDML2nd/dmlGate/dmlPass/dmlAng1/
%             dmlAng2/dmlW1Db/dmlW2Db/dmlW1/dmlW2 字段; 通过判定的第二目标已追加,
%             isDML2nd=1; dmlW1/dmlW2 为线性 LS 幅度, 供 RCS 推算)
%     detCfg  原样返回
%     dmlCfg  .figSeed   指定点的 DML 过程图(4 格; 未找到指定点时为 [])
%             .figCloud  更新点云的 6 子图 DOA 图(同 fig14 口径, 含 DML 第二目标)
%             .nDML2nd   判定出 DML 第二目标的点数
%             .nGate     因门限(角度/x位置)跳过 DML 的点数
%             .nDone     实际做了 DML 细分的点数
%             .seedFound 指定点是否找到

    % ---------- 方位字典(最长行 36 通道, 共轭域; 与 DOA/OMP 方位维一致) ----------
    aziRange = AntCaliData.aziSteerVecRange;
    aziGrid  = aziRange(1) : aziRange(3) : aziRange(2);
    [~, aziRowBest] = max(sum(rowLut, 2));
    chanAzi  = find(rowLut(aziRowBest, :));
    A_dml    = conj(AntCaliData.aziSteerVec(chanAzi, :));    % 信号模型域字典
    fitFun   = @(th) steerFitAt(A_dml, aziGrid, th);         % 表格插值原子(按需单向量)
    DeltaR   = 3e8/(2*dims.B);                          % 距离分辨率 (m)

    % ---------- DML 分工限定(与 OMP 无缝衔接: DML 管 1.5~6°, OMP 管 ≥5°, 5~6° 双保险) ----------
    dmlSpanDeg = 8;                     % [旋钮] DML 搜索半跨度: 第二角度候选限于 DBF 主峰 ±此值内 (deg, 覆盖 6° 分离+裕量)
    sepLimDeg  = 6;                     % [旋钮] 两目标角度间隔上限 (deg): 超出判不通过(曾为 4°, 与 OMP 的 6° 之间留 4~6° 空档, 已改 6° 衔接)
    dmlStepDeg = 0.05;                  % [旋钮] 交替坐标细分步进 (deg): ±1° 内按此步进搜极小
    nSweep     = 3;                     % [旋钮] 交替坐标细分轮数
    nStart     = 20;                    % [旋钮] 残差匹配先取的峰数(再按 ±dmlSpanDeg 筛选)

    % ---------- DML 处理门限(降低虚假点: 只对近距离/小角度点做 DML) ----------
    aziLimDeg = 20;                     % [旋钮] 方位角门限: |DBF拟合方位角| 超过则跳过
    xLimM     = 25;                     % [旋钮] x 位置门限 (m): |x| 超过则跳过

    nTgtAll = numel(tgtOut);            % 输入点云总数(含 OMP 第二目标)
    nDML2nd = 0;  nGate = 0;  nDone = 0;
    seedFound = false;  seed = [];
    if nargin < 8, seedD = []; end
    if nargin < 7, seedR = []; end
    if nargin < 9 || isempty(drawFigs), drawFigs = true; end   % false = 批量扫描不画图
    haveSeed = ~isempty(seedR) && ~isempty(seedD);

    % 统一新增标记字段(先建字段, 再用 deal 对【全部元素】填默认值; 必须在循环前 ——
    % 循环内要读 isOMP2nd。注意: 只写 tgtOut(1).字段 = 值时, 其余元素该字段是空 [],
    % 空值会在 [tgtOut.字段] 拼接时丢元素, 云图三色圈掩码就会错位)
    tgtOut(1).isDML2nd = 0;    tgtOut(1).dmlGate = 1;   tgtOut(1).dmlPass = 0;
    tgtOut(1).dmlAng1  = NaN;  tgtOut(1).dmlAng2 = NaN;
    tgtOut(1).dmlW1Db  = NaN;  tgtOut(1).dmlW2Db = NaN;
    tgtOut(1).dmlW1    = [];   tgtOut(1).dmlW2  = [];            % 线性幅度(RCS 推算用)
    [tgtOut.isDML2nd] = deal(0);
    [tgtOut.dmlGate]  = deal(1);
    [tgtOut.dmlPass]  = deal(0);
    [tgtOut.dmlAng1]  = deal(NaN);  [tgtOut.dmlAng2] = deal(NaN);
    [tgtOut.dmlW1Db]  = deal(NaN);  [tgtOut.dmlW2Db] = deal(NaN);
    [tgtOut.dmlW1]    = deal([]);   [tgtOut.dmlW2]   = deal([]);

    % ---------- 逐点 DML(跳过 OMP 第二目标; 门限内的点才做) ----------
    for t = 1:nTgtAll
        if tgtOut(t).isOMP2nd == 1
            tgtOut(t).dmlGate = 0;                       % OMP 第二目标不做 DML
            continue;
        end
        azi1 = double(tgtOut(t).azimuth);                % DBF 拟合方位角(该点第一目标)
        rngM = double(tgtOut(t).rangIdx) * DeltaR;       % 距离 (m)
        eleM = double(tgtOut(t).elevation);
        xPosT = rngM * cosd(eleM) * sind(azi1);          % x 位置 (m)

        % ---- 门限检查: |方位| ≤ 20° 且 |x| ≤ 25 m; 超出 → 保留 DOA 原结果 ----
        if abs(azi1) > aziLimDeg || abs(xPosT) > xLimM
            tgtOut(t).dmlGate = 0;
            nGate = nGate + 1;
            continue;
        end

        x  = double(detCfg.doaChanIn(t, :)).';
        x  = x(chanAzi);                                 % 方位维 36 通道
        C0 = x' * x;                                     % 快照总功率
        b  = A_dml' * x;                                 % 各角度匹配向量 [M×1]
        nDone = nDone + 1;

        % ---- K=1: 最佳单角度(残差 = C0 − |b(p)|²/36, 原子单位模) ----
        [J1s, p1] = max(abs(b).^2 / 36);
        J1 = C0 - J1s;
        thS = aziGrid(p1);                               % 强目标角度(所有起点的 θ1)

        % ---- 残差匹配曲线 + 候选峰(第二角度起点, 限于 DBF 主峰 ±dmlSpanDeg) ----
        a1ls = (A_dml(:, p1)' * x) / (A_dml(:, p1)' * A_dml(:, p1));
        res  = x - A_dml(:, p1) * a1ls;
        c2   = A_dml' * res;
        c2Db = 20*log10(abs(c2) + eps);
        [~, pkIall] = searchPeaksPwr_pipe(c2Db, nStart);       % 多取一些峰再按范围筛选
        pkI = pkIall(abs(aziGrid(pkIall) - thS) <= dmlSpanDeg);% 限于 DBF 主峰 ±dmlSpanDeg
        if isempty(pkI), pkI = p1; end                         % 兜底: 无合格峰则用强目标自身
        nCand = numel(pkI);                              % 实际可用峰数(可能少于 nStart)

        % ---- 多起点交替坐标细分(每起点: 两角度联合交替 ±1° 内 dmlStepDeg ×nSweep 轮) ----
        bestJ = inf;  bTh1 = [];  bTh2 = [];  bW = [];  bSt = 0;
        for st = 1:nCand
            th1 = thS;  th2 = aziGrid(pkI(st));          % 起点: 强目标 + 第 st 个候选峰
            for sw = 1:nSweep                            % 交替坐标细分 ×nSweep 轮
                f1  = th1 + (-1 : dmlStepDeg : 1);
                Jf1 = zeros(size(f1));
                for j1 = 1:numel(f1)
                    A2  = [fitFun(f1(j1)), fitFun(th2)];
                    aS  = inv2x2(A2' * A2) * (A2' * x);
                    Jf1(j1) = C0 - real(x' * A2 * aS);
                end
                [~, j1m] = min(Jf1);  th1 = f1(j1m);
                f2  = th2 + (-1 : dmlStepDeg : 1);
                Jf2 = zeros(size(f2));
                for j2 = 1:numel(f2)
                    A2  = [fitFun(th1), fitFun(f2(j2))];
                    aS  = inv2x2(A2' * A2) * (A2' * x);
                    Jf2(j2) = C0 - real(x' * A2 * aS);
                end
                [~, j2m] = min(Jf2);  th2 = f2(j2m);
            end
            A2f = [fitFun(th1), fitFun(th2)];
            aS  = inv2x2(A2f' * A2f) * (A2f' * x);
            Jf  = C0 - real(x' * A2f * aS);
            if Jf < bestJ
                bestJ = Jf;  bTh1 = th1;  bTh2 = th2;  bW = aS;  bSt = st;
            end
        end

        % ---- 角度升序整理(权重与升序角度对齐) ----
        ang1 = min(bTh1, bTh2);  ang2 = max(bTh1, bTh2);
        ampW = bW;
        if bTh1 > bTh2, ampW = ampW([2 1]); end

        % ---- 第二目标判定(残差改善 > 门限 且 间隔可分且限内) ----
        imprDb = 10 * log10(max(J1, eps) / max(bestJ, eps));
        dmlPass = (imprDb > thrDmlDb) && (abs(ang1 - ang2) > 0.2) && (abs(ang1 - ang2) <= sepLimDeg);

        % ---- 写回逐点结果 ----
        tgtOut(t).dmlAng1 = ang1;  tgtOut(t).dmlAng2 = ang2;
        tgtOut(t).dmlW1Db = 20*log10(abs(ampW(1)) + eps);
        tgtOut(t).dmlW2Db = 20*log10(abs(ampW(2)) + eps);
        tgtOut(t).dmlW1   = ampW(1);  tgtOut(t).dmlW2 = ampW(2);   % 线性幅度(RCS 推算用)
        tgtOut(t).dmlPass = dmlPass;

        % ---- 通过: 第一目标方位替换为 DML 第一角度, 追加第二目标进点云 ----
        if dmlPass
            nDML2nd = nDML2nd + 1;
            tgtOut(t).azimuth = ang1;                        % 第一目标: DBF角 → DML角
            newT          = tgtOut(t);                       % 继承同距同速/俯仰等字段
            newT.ID       = max([tgtOut.ID]) + 1;
            newT.azimuth  = ang2;                            % 第二目标方位
            newT.isDML2nd = 1;
            newT.dmlPass  = 0;                               % 追加点不是"被替换的母点", 清标记(否则计数翻倍)
            tgtOut(end+1) = newT;                            % #ok<AGROW> 追加段在云图用 is2 标记
        end

        % ---- 指定点(R/D 序号, 容差 ±1 bin) → 保存 DML 过程数据 ----
        if haveSeed && ~seedFound && ...
                abs(tgtOut(t).rangIdx - seedR) <= 1 && abs(tgtOut(t).dopplerIdx - seedD) <= 1
            seedFound = true;
            seed = struct('t', t, 'rngIdx', tgtOut(t).rangIdx, 'dopIdx', tgtOut(t).dopplerIdx, ...
                          'azi1', azi1, 'ang1', ang1, 'ang2', ang2, 'pass', dmlPass, ...
                          'imprDb', imprDb, 'nCand', nCand, 'bSt', bSt, ...
                          'bDb', 20*log10(abs(b)/max(abs(b)) + eps), 'c2Db', c2Db, ...
                          'pkI', pkI, 'x', x, 'W', ampW);
        end
    end
    if haveSeed && ~seedFound
        fprintf('DML 提示: RD(%g, %g) 处未找到确认点, 过程图跳过\n', seedR, seedD);
    end

    % ---------- 图窗公共尺寸(与 OMP/DOA 模块一致) ----------
    screenSize = get(0, 'ScreenSize');
    figW = screenSize(3) * 0.98;  figH = screenSize(4) * 0.80;
    figPos = [(screenSize(3)-figW)/2, (screenSize(4)-figH)/2, figW, figH];

    % ---------- 指定点的 DML 过程图(图24, 4 格; 未指定/未找到则为空) ----------
    figSeed = [];
    if seedFound && drawFigs
        vTxt = '不通过';  if seed.pass, vTxt = '通过'; end
        figSeed = figure('Position', figPos, 'Name', 'fig31_DML_process');
        tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

        nexttile;   % DBF 谱(单峰难分) + DML 两角度竖线
        plot(aziGrid, seed.bDb, 'b-'); hold on;
        xline(seed.ang1, 'b--', sprintf('%.2f°', seed.ang1));
        xline(seed.ang2, 'r--', sprintf('%.2f°', seed.ang2));
        grid on; xlim([min(aziGrid), max(aziGrid)]); ylim([-40 2]);
        xlabel('方位角 (deg)'); ylabel('归一化幅度 (dB)');
        title(sprintf('常规 DBF 谱 (DML 反演: %.2f° / %.2f°)', seed.ang1, seed.ang2));

        nexttile;   % 残差匹配曲线 + 全部候选起点 + 胜出起点
        plot(aziGrid, seed.c2Db - max(seed.c2Db), 'r-'); hold on;
        plot(aziGrid(seed.pkI), seed.c2Db(seed.pkI) - max(seed.c2Db), 'co', 'MarkerSize', 6);
        plot(aziGrid(seed.pkI(seed.bSt)), seed.c2Db(seed.pkI(seed.bSt)) - max(seed.c2Db), ...
            'gv', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
        yline(median(seed.c2Db) - max(seed.c2Db), 'k:');
        grid on; ylim([-40 10]);
        xlabel('方位角 (deg)'); ylabel('归一化幅度 (dB)');
        legend({'残差匹配曲线', '候选起点', '胜出起点'}, 'Location', 'best');
        title(sprintf('残差匹配: DBF主峰 ±%g°内 %d 个候选起点, 胜出 %.2f°', ...
            dmlSpanDeg, seed.nCand, aziGrid(seed.pkI(seed.bSt))));

        nexttile;   % DML 二维代价面(胜出解附近 ±2°)
        xS  = seed.x;  C0S = real(xS' * xS);
        fA  = seed.ang1 + (-2 : 0.1 : 2);  fB = seed.ang2 + (-2 : 0.1 : 2);
        Jm  = zeros(numel(fB), numel(fA));
        for ib = 1:numel(fB)
            for ia = 1:numel(fA)
                A2  = [fitFun(fA(ia)), fitFun(fB(ib))];
                aS  = inv2x2(A2' * A2) * (A2' * xS);
                Jm(ib, ia) = C0S - real(xS' * A2 * aS);
            end
        end
        JmDb = 10*log10(Jm / max(min(Jm(:)), eps) + eps);   % NaN/全零守卫(退化数据防 imagesc 崩溃)
        JmDb(~isfinite(JmDb)) = 0;
        imagesc(fA, fB, JmDb); axis xy; colorbar;
        hold on; plot(seed.ang1, seed.ang2, 'wo', 'MarkerSize', 10);
        plot(xlim, xlim, 'w:', 'LineWidth', 1);          % θ1=θ2 退化线(两原子重合, 恒为脊)
        xlabel('θ1 (deg)'); ylabel('θ2 (deg)');
        title(sprintf('DML 二维代价面 (谷底 %.2f° / %.2f°, 白虚线=θ1=θ2退化线)', seed.ang1, seed.ang2));

        nexttile;   % 重建检验
        xRef = max(abs(xS)) + eps;
        atom1R = steerFitAt(A_dml, aziGrid, seed.ang1);
        atom2R = steerFitAt(A_dml, aziGrid, seed.ang2);
        plot(20*log10(abs(xS)/xRef + eps), 'b-o', 'MarkerSize', 4); hold on;
        plot(20*log10(abs([atom1R, atom2R] * seed.W)/xRef + eps), 'r-*', 'MarkerSize', 4);
        grid on; ylim([-60 5]);
        xlabel('通道 (最长行 36)'); ylabel('相对快照峰值 (dB)');
        legend({'快照(归一化)', 'DML 双原子重建(同一参考)'}, 'Location', 'best');
        title('重建检验: 快照 vs 双原子重建');
        sgtitle(sprintf('DML 多起点超分辨: RD(%g,%g) → %.2f° / %.2f° | 残差改善 %.1f dB | 判定: %s', ...
            seedR, seedD, seed.ang1, seed.ang2, seed.imprDb, vTxt));

        fprintf('DML 指定点      : RD(%g,%g) → %.2f° / %.2f° | 残差改善 %.1f dB | 判定: %s\n', ...
            seedR, seedD, seed.ang1, seed.ang2, seed.imprDb, vTxt);
        fprintf('DML 幅度信息   : 两目标 LS 幅度 %.3g / %.3g (积累快照量纲, 供RCS推算) | 相对差 %.1f dB\n', ...
            abs(seed.W(1)), abs(seed.W(2)), 20*log10(abs(seed.W(1))/abs(seed.W(2))));
    end

    if ~drawFigs
        figCloud = [];                                   % 批量扫描: 不画图
    end
    if drawFigs
    % ---------- 更新点云的 6 子图 DBF 点迹图(同 fig28_DBF 口径, 含三种角度圈) ----------
    rngEst = cell2mat(arrayfun(@(q) double(q.range), tgtOut, 'UniformOutput', false));
    dopEst = cell2mat(arrayfun(@(q) double(q.dopplerIdx) - 1 - double(q.doppler_alpha), tgtOut, 'UniformOutput', false)); % 0基小数bin
    aziEst = cell2mat(arrayfun(@(q) double(q.azimuth),  tgtOut, 'UniformOutput', false));
    eleEst = cell2mat(arrayfun(@(q) double(q.elevation), tgtOut, 'UniformOutput', false));
    zEst   = rngEst .* sind(eleEst);
    xPos   = rngEst .* cosd(eleEst) .* sind(aziEst);
    yPos   = rngEst .* cosd(eleEst) .* cosd(aziEst);
    % 三种圈: 黑圈=OMP第二目标 / 红圈=DML第一目标(DBF角已被DML角替换) / 青圈=DML第二目标
    isOMP2 = logical([tgtOut.isOMP2nd]);                         % OMP 第二目标
    isD1   = false(1, numel(tgtOut));
    isD1(1:nTgtAll) = logical([tgtOut(1:nTgtAll).dmlPass]);      % DML 第一目标(原点段, 角度已替换)
    isD2   = logical([tgtOut.isDML2nd]);                         % DML 第二目标(追加段)
    rMax10 = ceil(dims.rngUse * 3e8/(2*dims.B) / 10) * 10;
    sym10  = @(v) max(ceil(max(abs(v)) / 10) * 10, 10);   % 全 0 时保底 ±10, 防 xlim([0,0]) 非法

    figCloud = figure('Position', figPos, 'Name', 'fig32_DBF_OMP_DML');
    tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    nexttile;
    scatter(xPos, yPos, 24, zEst, 'filled'); hold on;
    scatter(xPos(isD1), yPos(isD1), 40, 'ro', 'LineWidth', 1.5); hold on;
    scatter(xPos(isD2), yPos(isD2), 40, 'co', 'LineWidth', 1.5); hold on;
    scatter(xPos(isOMP2), yPos(isOMP2), 40, 'ko', 'LineWidth', 1.5); grid on; colorbar;
    legend({'OMP第二目标(黑圈)','DML第一目标(红圈)','DML第二目标(青圈)'}, 'Location', 'best');
    xlim([-sym10(xPos), sym10(xPos)]); ylim([0, rMax10]);
    xlabel('x (m)'); ylabel('y (m)'); title('x-y 位置');
    nexttile;
    scatter(aziEst, rngEst, 24, zEst, 'filled'); hold on;
    scatter(aziEst(isD1), rngEst(isD1), 40, 'ro', 'LineWidth', 1.5); hold on;
    scatter(aziEst(isD2), rngEst(isD2), 40, 'co', 'LineWidth', 1.5); hold on;
    scatter(aziEst(isOMP2), rngEst(isOMP2), 40, 'ko', 'LineWidth', 1.5); grid on; colorbar;
    xlim([-sym10(aziEst), sym10(aziEst)]); ylim([0, rMax10]);
    xlabel('方位角 (deg)'); ylabel('距离 (m)'); title('方位角-距离');
    nexttile;
    scatter(dopEst, rngEst/dims.DeltaR, 24, zEst, 'filled'); grid on; colorbar;
    xlim([0, dims.dopBinNum-1]); ylim([0, dims.rngUse]);
    xlabel('多普勒 bin'); ylabel('距离 bin'); title('多普勒-距离 (满bin)');
    nexttile;
    scatter(zEst, yPos, 24, zEst, 'filled'); hold on;
    scatter(zEst(isD1), yPos(isD1), 40, 'ro', 'LineWidth', 1.5); hold on;
    scatter(zEst(isD2), yPos(isD2), 40, 'co', 'LineWidth', 1.5); hold on;
    scatter(zEst(isOMP2), yPos(isOMP2), 40, 'ko', 'LineWidth', 1.5); grid on; colorbar;
    xlim([-sym10(zEst), sym10(zEst)]); ylim([0, rMax10]);
    xlabel('z (m)'); ylabel('y (m)'); title('z-y 位置');
    nexttile;
    scatter(eleEst, rngEst, 24, zEst, 'filled'); hold on;
    scatter(eleEst(isD1), rngEst(isD1), 40, 'ro', 'LineWidth', 1.5); hold on;
    scatter(eleEst(isD2), rngEst(isD2), 40, 'co', 'LineWidth', 1.5); hold on;
    scatter(eleEst(isOMP2), rngEst(isOMP2), 40, 'ko', 'LineWidth', 1.5); grid on; colorbar;
    xlim([-sym10(eleEst), sym10(eleEst)]); ylim([0, rMax10]);
    xlabel('俯仰角 (deg)'); ylabel('距离 (m)'); title('俯仰角-距离');
    nexttile;
    scatter3(xPos, yPos, zEst, 24, zEst, 'filled'); hold on;
    scatter3(xPos(isD1), yPos(isD1), zEst(isD1), 40, 'ro', 'LineWidth', 1.5); hold on;
    scatter3(xPos(isD2), yPos(isD2), zEst(isD2), 40, 'co', 'LineWidth', 1.5); hold on;
    scatter3(xPos(isOMP2), yPos(isOMP2), zEst(isOMP2), 40, 'ko', 'LineWidth', 1.5); grid on;
    xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)'); view(3);
    title('点云 xyz (可旋转)');
    sgtitle(sprintf(['DML 后点云: 原 %d 点 + DML第二目标 %d 点 | ' ...
        '黑圈=OMP第二目标, 红圈=DML第一目标(DBF角已替换), 青圈=DML第二目标'], nTgtAll, nDML2nd));
    end                                             % if drawFigs(云图段)

    dmlCfg = struct('figSeed', figSeed, 'figCloud', figCloud, 'nDML2nd', nDML2nd, ...
                    'nGate', nGate, 'nDone', nDone, 'seedFound', seedFound);
end
