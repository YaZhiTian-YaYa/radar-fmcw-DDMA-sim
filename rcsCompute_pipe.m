function [tgtOut, rcsCfg, figRcs] = rcsCompute_pipe(tgtOut, detCfg, AntCaliData, rowLut, dims, prm, targets, rcsCfg, winType, accumCfg, fftCfg, zRowV, lambdaV, rfInfo)
%RCSCOMPUTE_PIPE RCS 估计: 对全部点云逐点算 RCS, 写回点云并绘图(双域, 2026-11-04)
%   [tgtOut, rcsCfg, figRcs] = rcsCompute_pipe(tgtOut, detCfg, AntCaliData, rowLut, ...
%       dims, prm, targets, rcsCfg, winType, accumCfg, fftCfg, zRowV, lambdaV)
%   【双域口径(用户 2026-11-04: 弱目标 RCS 用 TXRX 相参累积图)】
%     - 独立点(非 OMP/DML 配对成员): 相参累积图读数 —— 方位谱(检测峰 bin±1 邻域)×
%       距离维 5 点峰-插值 × 多普勒恢复比 / 相参噪底(标定遍 cohFloorPow+拟合滤形),
%       俯仰≠0 的 z 行失配损失按 cohGainDbf(θ̂el)−cohGainDbf(0) 补回(用点云自身俯仰
%       估计, 生产合法)。动因: 非相参图无积累增益(8TX 对齐+8RX 功率平均不改均值 SNR),
%       弱目标峰单元噪声主导、RCS 起伏大; 相参图(64 元方位 FFT)+~18 dB 弱目标可读。
%       标定常数 K0CohSeg(角反 ele=0 无损失, 读数同法锚定)。
%     - 配对成员(OMP/DML 第二目标与母点): α 域联合 LS 权重(ompW1/2, dmlW1/2) ——
%       ±1° 级 DML 对在方位谱 1.4°/bin 分辨率内互染, 图读数不可分; 联合权重才有分离。
%       幅度来源(均为"目标复幅度 α"同一量纲, 积累快照单位):
%       标定常数 K0Seg(α 域角反锚定)。
%   RCS 公式(SNR 比值, 分段):
%     σ_dB = SNR_dB + 40log10(R) + K0Seg/K0CohSeg(段) + [G(0°) − G(θ)]
%     独立点 SNR = 相参读数dB − (cohFloorDb + 拟合滤形) − 俯仰失配补偿
%     配对点 SNR = 2·α_dB − 10log10(该距离门噪底参考)
%   信号与噪底同 bin 同滤形 → SNR 与 IF 滤波器形状无关; 分段按 R 与 segR 比较;
%   写回字段: rcsDb / rcsSeg / rcsAmpDb / rcsSnrDb / rcsConf / rcsDom(1=相参图, 0=α域)。
%   绘图(图35, 2×2); 打印逐真值目标对比表(含域与置信标志)。

    % ---------- 方位字典(最长行 36 通道, 共轭域; α 域分支用) ----------
    aziRange = AntCaliData.aziSteerVecRange;
    aziGrid  = aziRange(1) : aziRange(3) : aziRange(2);
    [~, aziRowBest] = max(sum(rowLut, 2));
    chanAzi  = find(rowLut(aziRowBest, :));
    A_dbf    = conj(AntCaliData.aziSteerVec(chanAzi, :));    % 信号模型域字典

    % ---------- 相参域可用性(标定量齐备 + 本遍数据在) ----------
    useCohAvail = isfield(rcsCfg, 'K0CohSeg') && ~any(isnan(rcsCfg.K0CohSeg)) && ...
                  isfield(rcsCfg, 'cohFloorDb') && nargin >= 13 && ~isempty(lambdaV) && ...
                  isfield(accumCfg, 'aziFftSpec');
    if useCohAvail
        gCoh0 = cohGainDbf(zRowV, lambdaV, 0);                 % θel=0 相参增益 = 10lg(64)
        fprintf(['RCS 双域口径     : 独立点 = TXRX 相参累积图(方位峰bin±1×距离峰插值×多普勒恢复, 俯仰损失补 cohGainDbf) | 配对点 = α 域联合权重' newline]);
    else
        fprintf(['RCS 双域口径     : 相参标定量缺失, 全部回退 α 域 LS' newline]);
    end

    % ---------- PLL 线损补偿(2026-11-04, 两分支同补): 逐 chirp 相位噪声把谱线功率按
    % e^{−σφ²} 劈走(多普勒去相关), σφ = 2π·σ_δf·τ, τ̂=2R̂/c —— σ_δf = pnTrackCoef×B/Ts
    % (芯片 PLL 跟踪项) + hopErrHz(跳频项), 均 ∝R²: 角反 τ≈23ns 免疫(K0 干净)、远目标
    % 被压。消融实测@200m 代价 1.13dB ≈ 线损式 4.343σφ²=1.22dB(自扩散项 0.5 属模型
    % 高估, 不补)。生产合法: 芯片 PLL 数据 + 自估距离; snrAnaCalc 同式(sig2) ----
    pllCoefHz = 0;  hopHz = 0;
    if nargin >= 14 && ~isempty(rfInfo)
        if isfield(rfInfo, 'pllTrackOn') && rfInfo.pllTrackOn && isfield(rfInfo, 'pnTrackCoef')
            pllCoefHz = rfInfo.pnTrackCoef * prm.B / prm.Ts;   % σ_δf 跟踪项 (Hz)
        end
        if isfield(rfInfo, 'hopErrOn') && rfInfo.hopErrOn && isfield(rfInfo, 'hopErrHz')
            hopHz = rfInfo.hopErrHz;                           % 跳频误差 σ (Hz)
        end
    end
    pllCompDb = @(Rq) 4.343 * ( (2*pi*pllCoefHz*2*Rq/3e8).^2 + (2*pi*hopHz*2*Rq/3e8).^2 );
    if pllCoefHz + hopHz > 0
        fprintf(['RCS PLL 补偿    : 线损 4.343σφ²(σ_δf=%.1f kHz + 跳频 %.3g Hz), @200m %+.2f dB, @300m %+.2f dB' newline], ...
            pllCoefHz/1e3, hopHz, pllCompDb(200), pllCompDb(300));
    end

    % ---------- 亚bin窗损补偿(α 域幅度: 配对点生产 + 独立点对照; SVA 无固定核不补) ----------
    % 2026-11-04 配对优化: 旧口径 winScal(检测插值距离)+compDop(−0.37 统计) 弃用 —— 前者
    % 有幻影(插值散布 ±0.75bin, 实测虚补 +2.3dB), 后者只对均值。改为逐点精确恢复(见
    % 循环内 leafRDb/leafDDb): 距离维相参图峰剖面 hannPk3 恢复比 + 多普勒维 8 码线逐线恢复
    doWinComp = nargin >= 9 && strcmpi(winType, 'Hanning');
    if doWinComp
        fprintf(['RCS 窗损补偿   : (α域幅度) 逐点精确恢复 = 距离维相参图峰剖面 hannPk3 + 多普勒维 8 码线逐线均值(旧 winScal/compDop 统计口径已弃)' newline]);
    else
        fprintf(['RCS 窗损补偿   : (α域幅度) %s 无固定核不补偿(口径如实)' newline], 'SVA 自适应窗');
    end

    rcsLowConfSepDeg = 3;    % [旋钮] RCS 低置信角间隔阈值(deg): OMP/DML 对内间隔小于此值
                             %   (≈1.6×波束宽, 36 元半波长 ULA 波束宽 ≈1.84°)时联合 LS 分离
                             %   病态(κ≈3.5@1°), 单实现 RCS 散布 ±1.5dB —— 打低置信标志供下游加权
    nTgt = numel(tgtOut);
    if nTgt == 0                                   % 空点云(无检测点/目标全出带): 跳过估计
        warning('RCS 估计: 点云为空, 跳过逐点估计');
        screenSize = get(0, 'ScreenSize');
        figRcs = figure('Position', [screenSize(3)/2-300, screenSize(4)/2-150, 600, 300], ...
                        'Name', 'fig35_RCS');
        text(0.5, 0.5, '无检测点, RCS 估计跳过', 'HorizontalAlignment', 'center', 'FontSize', 14);
        return;
    end
    rcsDb  = NaN(nTgt, 1);
    rcsSeg = zeros(nTgt, 1);
    ampDb  = NaN(nTgt, 1);
    snrDb  = NaN(nTgt, 1);
    domFlg = zeros(nTgt, 1);
    rcsAlp = NaN(nTgt, 1);                       % α 域对照值(独立点=旧口径, 仅供打印对比)

    for t = 1:nTgt
        % 配对成员判定(四类): OMP第二/DML第二/OMP母点/DML母点 —— 这些点保持 α 域
        isPairMem = tgtOut(t).isOMP2nd == 1 || tgtOut(t).isDML2nd == 1 || ...
                    (tgtOut(t).ompHas2nd && ~isnan(tgtOut(t).ompW1Db)) || ...
                    (tgtOut(t).dmlPass && ~isnan(tgtOut(t).dmlW1Db));
        kOk = tgtOut(t).rangIdx >= 3 && tgtOut(t).rangIdx <= dims.rngUse - 2;
        R = double(tgtOut(t).range);               % 插值距离 (m)
        seg = find(R >= rcsCfg.segR(1:end-1) & R < rcsCfg.segR(2:end), 1);
        if isempty(seg), seg = numel(rcsCfg.K0Seg); end
        az = double(tgtOut(t).azimuth);  el = double(tgtOut(t).elevation);
        % ---- α 域逐点精确叶损因子(2026-11-04 配对优化): α 幅度取自单 (距离,被抓多普勒) bin
        % 的快照, 带两维窗叶损 —— 旧口径 winScal(检测插值距离) 有幻影(插值散布 ±0.75bin,
        % 实测虚补 +2.3), compDop=−0.37 只是统计均值。改为与相参读数同机制的精确恢复:
        % 距离维 = 相参图该 cell 峰bin 距离剖面 hannPk3 峰/本 bin 功率比; 多普勒维 = 8 根
        % DDMA 码线逐线恢复等权平均。对内两成员共享同 cell 因子(距离偏移与角 bin 无关)。
        leafRDb = 0;  leafDDb = 0;
        if doWinComp && kOk && useCohAvail
            sbL = mod(double(tgtOut(t).dopplerIdx) - 1, dims.subBins) + 1;
            if isfield(tgtOut, 'cohPeakIdx') && ~isempty(tgtOut(t).cohPeakIdx)
                b0L = double(tgtOut(t).cohPeakIdx);
            else
                b0L = double(accumCfg.cohPeakPos(tgtOut(t).rangIdx, sbL));  % 追加点的 cell 峰
            end
            c5L = double(abs(squeeze(accumCfg.aziFftSpec(b0L, tgtOut(t).rangIdx-2:tgtOut(t).rangIdx+2, sbL))).^2);
            [~, pjL] = max(c5L(2:4));  pjL = pjL + 1;
            [~, pkL, okL] = hannPk3(c5L(pjL-1:pjL+1), dims.rngPts);
            if okL
                leafRDb = 10*log10(pkL / c5L(3));       % 本 bin → 真峰 功率恢复比
            end
            if isfield(accumCfg, 'chirpPosMap')
                binsL = double(accumCfg.chirpPosMap(tgtOut(t).rangIdx, sbL, :));
            else
                binsL = double(tgtOut(t).dopplerIdx);
            end
            recsL = ones(numel(binsL), 1);
            for ttL = 1:numel(binsL)
                dbL = double(binsL(ttL));
                dJL = mod((dbL-3):(dbL+1), dims.dopBinNum) + 1;
                q5L = squeeze(sum(abs(fftCfg.data_sva(tgtOut(t).rangIdx, :, dJL)).^2, 2));
                [~, qjL] = max(q5L(2:4));  qjL = qjL + 1;
                [~, QpkL, okDL] = hannPk3(q5L(qjL-1:qjL+1), dims.dopBinNum);
                if okDL, recsL(ttL) = QpkL / q5L(3); end
            end
            leafDDb = 10*log10(mean(recsL));
        end
        if useCohAvail && ~isPairMem && kOk && isfield(tgtOut, 'cohPeakIdx') && ~isempty(tgtOut(t).cohPeakIdx)
            % ---- 独立点: 相参累积图读数(与角反相参锚定同法 → K0Coh 相消) ----
            [pPkDb, snrC] = cohReadAt(accumCfg, fftCfg, rcsCfg, tgtOut(t), dims);
            % 俯仰 z 行失配补偿: 读数增益 = gDb(θ̂el) ≤ gDb(0)=10lg64, 补回 ele=0 口径
            % (解析同式: snrAnaCalc 的 cohGainDb 实参, 一维方位 FFT 不含 z 行导向)
            gEl = cohGainDbf(zRowV, lambdaV, el) - gCoh0;
            snrVal = snrC - gEl;                             % gEl ≤ 0 → 补回 ≥0 dB
            ampVal = pPkDb;                                  % 相参功率读数 dB(噪底前)
            K0v    = rcsCfg.K0CohSeg;
            domFlg(t) = 1;
        else
            % ---- 配对成员/回退: α 域联合权重或单原子 LS 投影 ----
            if tgtOut(t).isOMP2nd == 1, aDbl = tgtOut(t).ompW2Db;
            elseif tgtOut(t).isDML2nd == 1, aDbl = tgtOut(t).dmlW2Db;
            elseif tgtOut(t).ompHas2nd && ~isnan(tgtOut(t).ompW1Db)
                aDbl = tgtOut(t).ompW1Db;              % OMP 对的母点: 联合第一分量(去邻目标泄漏)
            elseif tgtOut(t).dmlPass && ~isnan(tgtOut(t).dmlW1Db)
                aDbl = tgtOut(t).dmlW1Db;              % DML 对的母点: 联合第一分量
            else                                       % 回退: 当前方位的 LS 投影幅度
                x = double(detCfg.doaChanIn(t, chanAzi)).';
                atom = steerFitAt(A_dbf, aziGrid, az);
                aDbl = 20*log10(abs((atom' * x) / (atom' * atom)) + eps);
            end
            aDbl = aDbl + leafRDb + leafDDb;            % 精确叶损恢复(旧 winScal/compDop 已弃)
            snrVal = aDbl - (rcsCfg.flatDb + rcsCfg.shapeAt(tgtOut(t).rangIdx));
            ampVal = aDbl;
            K0v    = rcsCfg.K0Seg;
        end
        % σ_dB = SNR + 40log10(R) + 段常数 + PLL线损补偿 + 增益VS角度修正(暂 0)
        rcsDb(t)  = snrVal + 40*log10(R) + K0v(seg) + pllCompDb(R) + rcsCfg.angGainDb(az, el);
        rcsSeg(t) = seg;
        ampDb(t)  = ampVal;
        snrDb(t)  = snrVal;
        if domFlg(t) == 1
            % α 域对照(同点旧口径: 单原子 LS + 窗补 + K0Seg) —— 仅打印对比, 不进生产值;
            % 弱目标两域偏差差 = 相参图收益的直接证据
            x = double(detCfg.doaChanIn(t, chanAzi)).';
            atom = steerFitAt(A_dbf, aziGrid, az);
            aAl  = 20*log10(abs((atom' * x) / (atom' * atom)) + eps) + leafRDb + leafDDb;
            rcsAlp(t) = (aAl - (rcsCfg.flatDb + rcsCfg.shapeAt(tgtOut(t).rangIdx))) + ...
                        40*log10(R) + pllCompDb(R) + rcsCfg.K0Seg(seg) + rcsCfg.angGainDb(az, el);
        else
            rcsAlp(t) = rcsDb(t);                            % α 域点: 对照=生产值
        end
    end

    for t = 1:nTgt
        tgtOut(t).rcsDb    = rcsDb(t);
        tgtOut(t).rcsSeg   = rcsSeg(t);
        tgtOut(t).rcsAmpDb = ampDb(t);
        tgtOut(t).rcsSnrDb = snrDb(t);
        tgtOut(t).rcsDom   = domFlg(t);            % 1=相参图口径, 0=α域口径(下游/打印用)
        tgtOut(t).rcsDbAlpha = rcsAlp(t);          % α 域对照(独立点=旧口径估计, 配对点=生产值)
        % RCS 置信标志: 配对成员角间隔 < 阈值 → 低置信(0), 独立点/宽间隔 → 高置信(1)
        sepP = inf;
        if isfield(tgtOut, 'ompSepDeg') && ~isempty(tgtOut(t).ompSepDeg)
            sepP = tgtOut(t).ompSepDeg;                    % OMP 对间隔(配对时双点同存)
        elseif isfield(tgtOut, 'dmlSepDeg') && ~isempty(tgtOut(t).dmlSepDeg)
            sepP = tgtOut(t).dmlSepDeg;                    % DML 对间隔
        end
        tgtOut(t).rcsConf = double(sepP >= rcsLowConfSepDeg);
    end

    % ---------- 逐真值目标对比(场景 targets vs 点云峰值估计; 全部点含第二目标) ----------
    % 用簇内【最强点】代表目标幅度(簇边缘/副瓣检测点弱 10~20 dB, 中位数会拉偏低);
    % 同 RD 的多目标(如 3°/4.5° 对)用方位约束 ±4° 区分 —— 匹配范围含 OMP/DML
    % 第二目标点(它们才是分离出来的那个目标的代表点)
    fprintf('RCS 逐目标对比  : 真值 vs 点云峰值估计(含OMP/DML第二目标, RD+方位邻域匹配)\n');
    claimedRcs = false(1, nTgt);   % 一对一分配: 已被某真值认领的点不再给下一个(角距<窗口时两真值会抢同一点)
    if ~isfield(prm, 'stepDir') || isempty(prm.stepDir), sdR = 1; else, sdR = prm.stepDir; end
    for k = 1:numel(targets)
        dopC = mod(targets(k).gate_D + sdR*prm.cr*prm.B/(dims.dopBinNum-1) * ...
                   (2*targets(k).gate_R*3e8/(2*prm.B))/3e8 * dims.dopBinNum, dims.dopBinNum);
        idx = find(abs([tgtOut.rangIdx] - (targets(k).gate_R+1)) <= 2 & ...
                   abs([tgtOut.dopplerIdx] - (round(dopC)+1)) <= 3 & ...
                   abs([tgtOut.azimuth] - targets(k).azi) <= 4);
        if isempty(idx), continue; end           % 该目标未被检测到(如超量程演示目标)
        % 先把方位门限收紧到真值 ±1° 再取最强(窄门能分开 1.5° 对; 强者优先排除
        % 恰好压在真值方位上的弱鬼影), 窄门无点才回退 ±4° 宽门
        idxN = idx(abs([tgtOut(idx).azimuth] - targets(k).azi) <= 1.0 & ~claimedRcs(idx));
        if isempty(idxN)                          % 无未认领点: 回退全部窄门点(含已认领)
            idxN = idx(abs([tgtOut(idx).azimuth] - targets(k).azi) <= 1.0);
        end
        if ~isempty(idxN), idx = idxN; end
        % 一对一: 窗内未认领者取最强, 认领后不再给下一个真值(配对两点不抢同一点)
        [~, im] = max(rcsDb(idx));
        claimedRcs(idx(im)) = true;
        pkDb = rcsDb(idx(im));
        confTxt = '低置信';  if tgtOut(idx(im)).rcsConf == 1, confTxt = '高置信'; end
        if tgtOut(idx(im)).rcsDom == 1, domTxt = '相参图'; else, domTxt = 'α域  '; end
        alpTxt = '';                                        % 独立点并列旧口径偏差(相参收益证据)
        if tgtOut(idx(im)).rcsDom == 1
            alpTxt = sprintf(', α对照偏差 %+5.1f dB', rcsAlp(idx(im)) - targets(k).rcsDb);
        end
        fprintf('   目标%d: R=%6.1f m 真值 %5.1f dBsm → 峰值估计 %6.1f dBsm (n=%d, 偏差 %+5.1f dB, %s, %s%s)\n', ...
            k, targets(k).gate_R*3e8/(2*prm.B), targets(k).rcsDb, pkDb, numel(idx), ...
            pkDb - targets(k).rcsDb, confTxt, domTxt, alpTxt);
    end

    % ---------- 绘图(图35, 2×2) ----------
    rngM  = cell2mat(arrayfun(@(q) double(q.range), tgtOut, 'UniformOutput', false));
    aziE  = cell2mat(arrayfun(@(q) double(q.azimuth), tgtOut, 'UniformOutput', false));
    eleE  = cell2mat(arrayfun(@(q) double(q.elevation), tgtOut, 'UniformOutput', false));
    xPos  = rngM .* cosd(eleE) .* sind(aziE);
    yPos  = rngM .* cosd(eleE) .* cosd(aziE);
    segTxt = {'近场', '远场'};
    rMax10 = ceil(max(rngM)/10)*10;

    screenSize = get(0, 'ScreenSize');
    figW = screenSize(3) * 0.98;  figH = screenSize(4) * 0.80;
    figPos = [(screenSize(3)-figW)/2, (screenSize(4)-figH)/2, figW, figH];
    figRcs = figure('Position', figPos, 'Name', 'fig35_RCS');
    tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile;   % RCS-距离(颜色=方位角; 竖线=分段边界; 星=标定角反真值位置)
    scatter(rngM, rcsDb, 24, aziE, 'filled'); hold on;
    plot(rcsCfg.segR(2)*[1 1], ylim, 'k--');
    plot(rcsCfg.corner.R, rcsCfg.corner.sigDb, 'r*', 'MarkerSize', 12);  % 角反不在目标遍场景, 星=标定真值
    grid on; xlim([0, rMax10]);
    xlabel('距离 (m)'); ylabel('RCS 估计 (dBsm)');
    legend({'点云', sprintf('分段边界 %.1f m', rcsCfg.segR(2)), '标定角反'}, 'Location', 'best');
    title('RCS-距离 (颜色=方位角)');

    nexttile;   % RCS-方位
    scatter(aziE, rcsDb, 24, rngM, 'filled'); grid on;
    xlabel('方位角 (deg)'); ylabel('RCS 估计 (dBsm)');
    title('RCS-方位角 (颜色=距离)');

    nexttile;   % RCS 直方图
    histogram(rcsDb, 30); grid on;
    xlabel('RCS 估计 (dBsm)'); ylabel('点数');
    title(sprintf('RCS 分布 (中位 %.1f dBsm)', median(rcsDb)));

    nexttile;   % x-y 位置(颜色=RCS)
    scatter(xPos, yPos, 24, rcsDb, 'filled'); hold on;
    plot(0, rcsCfg.corner.R, 'r*', 'MarkerSize', 12);   % 角反 0° 布设位(x=0, y=R), 标定参考
    axis equal; grid on; colorbar;
    xlabel('x (m)'); ylabel('y (m)');
    title('点云 x-y 位置 (颜色=RCS, 星=角反)');
    sgtitle(sprintf('RCS 估计 (角反 %.0f dBsm 标定, K0=%.1f/K0Coh=%.1f dB): 独立点相参图 %d / α域 %d 点', ...
        rcsCfg.corner.sigDb, rcsCfg.K0Seg(1), rcsCfg.K0CohSeg(1), sum(domFlg == 1), sum(domFlg == 0)));
end

function [pPkDb, snrDbV] = cohReadAt(accumCfg, fftCfg, rcsCfg, pt, dims)
%COHREADAT 相参累积图读数(与 rcsCalibrate 角反相参锚定同法, 常数在 K0Coh 相消):
%   方位谱(检测峰 bin±1 邻域) × 距离维 5 点峰-插值 × 多普勒恢复比 → 功率dB / 相参SNR
%   pPkDb = 10lg(峰插值功率 × 多普勒恢复比);  snrDbV = pPkDb − (cohFloorDb + 拟合滤形)
k  = double(pt.rangIdx);  db = double(pt.dopplerIdx);
sb = mod(db - 1, dims.subBins) + 1;
b0 = double(pt.cohPeakIdx);                     % 检测链已存的本点方位峰 bin(俯仰压缩自洽)
nA = size(accumCfg.aziFftSpec, 1);              % 方位 FFT 长度(=dims.aziFftLen)
% 峰 bin ±1 邻域(回绕): 先逐 bin 距离维峰-插值, 再角维 3 点幅度抛物线顶点
% (纯 max 残余 FFT 栅格纹波 ±0.7dB@1.4°/bin —— 离轴目标散布来源; 角反 bin 居中不受影响)
pk3 = zeros(3, 1);  wbb = 0;
for bb = mod((b0-2):(b0), nA) + 1
    wbb = wbb + 1;
    c5 = double(abs(squeeze(accumCfg.aziFftSpec(bb, k-2:k+2, sb))).^2);
    [~, pj] = max(c5(2:4));  pj = pj + 1;       % 峰-插值两步(同 peakFit: 先内三点最强)
    [~, pTk, okP] = hannPk3(c5(pj-1:pj+1), dims.rngPts);
    if okP, pk3(wbb) = pTk; else, pk3(wbb) = c5(3); end           % 肩点回退原格
end
a3a = sqrt(max(pk3, 0));                        % 角维幅度抛物线: 顶点 = A2 − (A1−A3)·d/4
dA = 2*(a3a(3) - a3a(1)) / (a3a(1) + 2*a3a(2) + a3a(3) + eps);
Cpk = (a3a(2) - (a3a(1) - a3a(3)) * max(min(dA, 1), -1) / 4).^2;
dRec = 1;                                       % 多普勒恢复比: 8 根码线逐线求(肩点回退 1)
if isfield(accumCfg, 'chirpPosMap')
    binsG = double(accumCfg.chirpPosMap(k, sb, :));    % 每 TX 各自被抓的 chirpPos bin
else
    binsG = db;                                 % 兜底: 单 bin(旧口径)
end
recs = ones(numel(binsG), 1);
for tt = 1:numel(binsG)                         % 各线叶损不同 → 逐线恢复等权平均(与角反
    dbG = double(binsG(tt));                    % 相参锚定同法, 残差相消)
    dJ = mod((dbG-3):(dbG+1), dims.dopBinNum) + 1;   % 5 点全 bin(1 起, 循环回绕)
    q5 = squeeze(sum(abs(fftCfg.data_sva(k, :, dJ)).^2, 2));   % RX 功率和 [5×1]
    [~, qj] = max(q5(2:4));  qj = qj + 1;
    [~, Qpk, okD] = hannPk3(q5(qj-1:qj+1), dims.dopBinNum);
    if okD, recs(tt) = Qpk / q5(3); end
end
if ~isempty(recs), dRec = mean(recs); end
pPkDb  = 10*log10(Cpk * dRec + eps);
snrDbV = pPkDb - (rcsCfg.cohFloorDb + rcsCfg.shapeAt(k));
end

function [d, pk, ok] = hannPk3(p, N)
%HANNPK3 3 点功率峰的汉宁精确恢复(2026-11-04, 替代纯抛物线):
%   位置 = 幅度域 Grandke 插值 δ=2(a₃−a₁)/(a₁+2a₂+a₃)(抛物线位置在 δ→0.5 偏短
%   ~0.25 bin); 峰功 = p0×10^(winScalDb(δ,N)/10)(抛物线峰功在 δ=0.5 低估 1.18dB
%   —— 分数 bin 目标系统性偏低, 默认场景 285.5 bin 即中招; δ=0.25/+0.02dB,
%   δ=0.5/+0.00dB 精确)。肩点判据同 para3p(非凹/越界 → 回退原格)
a   = sqrt(max(p, 0));
den = p(1) - 2*p(2) + p(3);
ok  = den < 0;                                % 非凹 = 肩点, 回退原格
d   = 2*(a(3) - a(1)) / (a(1) + 2*a(2) + a(3) + eps);
if ~ok
    d = 0;  pk = p(2);  return
end
d = max(min(d, 0.5), -0.5);                   % 半bin目标 d≈0.5 恰在边界: 截断而非拒绝
                                              %   (噪声越过 0.5 就回退会吃下整 1.4dB 叶损)
w  = hanning(N);  nn = (0:N-1).';
pk = p(2) * 10^(-20*log10(abs(w.'*exp(-1i*2*pi*d*nn/N)) / sum(w)) / 10);
end
