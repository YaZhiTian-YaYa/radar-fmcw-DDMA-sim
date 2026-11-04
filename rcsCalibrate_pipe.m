function rcsCfg = rcsCalibrate_pipe(tgtOut, detCfg, AntCaliData, rowLut, dims, prm, targets, calibInfo, noiseCfg, accumCfg, fftCfg, filtChain, winType)
%RCSCALIBRATE_PIPE RCS 标定: 用角反(targets 中 isCorner=1)标定链路常数(SNR 比值法)
%   rcsCfg = rcsCalibrate_pipe(tgtOut, detCfg, AntCaliData, rowLut, dims, prm, targets, calibInfo, noiseCfg, filtChain, winType)
%   filtChain(可选): 滤波链身份结构(主函数自 rfInfo 透传, 命中缓存时随缓存带) ——
%   含 hpf1Fc/hpf2Fc 时滤形直接用设计档位(两级一阶解析式), 曲线拟合仅作核查; 见下方根因注释
%   winType(可选): 2D-FFT 窗类型('Hanning' 时做亚bin窗损补偿, 与 rcsCompute 同口径)
%   原理(SNR 比值标定 —— 算法工程师口径, 不需要知道滤波器设计):
%     雷达方程 SNR 形式: SNR = |α|²/N = σ_dB − 40log10(R) + C → σ_dB = SNR_dB + 40log10(R) + K0
%     SNR_dB = α_dB(20log|α| ≡ 10log|α|², 即信号功率dB) − 噪底dB(该距离门)。
%     通带内目标: 信号与噪底同 bin 同滤形 → SNR 与 IF 滤波器形状无关, 零补偿。
%     【角反特例】角反在近距, 差拍落高通阻带: 信号被压 ~37dB, 但它那个门的
%     噪底被泄漏旁瓣抬高(不随滤形走) → 直接用本门噪底会错。修法(噪声曲线拟合法):
%     噪底曲线在干净过渡段(热噪主导)本身就是滤形 → 用两级高通模型拟合该段,
%     外推到角反门, 角反的噪底参考 = 平坦段噪底 + 拟合滤形(与信号同口径被压)。
%   K0 = σ_cal_dB − SNR_cal_dB − 40log10(R_cal), 角反布设距离按 bin×ΔR 真值取。
%   分段(近场/远场): segR = [0, RfarField, Inf], 每段一个 K0Seg(结构留好, 暂共用)。
%   天线增益VS角度: 暂无数据, 预留接口(rcsCfg.angGainDb 恒 0)。
%   输出:
%     rcsCfg  .K0Seg [seg×1] 各段标定常数 (dB)     .segR 分段边界 (m)
%             .flatDb 平坦段噪底 (dB)              .shapeAt(g) 拟合滤形句柄(g=距离bin, 1起)
%             .corner 角反标定记录  .angGainDb 增益VS角度句柄(暂恒0)

    % ---------- 方位字典(最长行 36 通道, 共轭域; 与 DOA/OMP/DML 一致) ----------
    aziRange = AntCaliData.aziSteerVecRange;
    aziGrid  = aziRange(1) : aziRange(3) : aziRange(2);
    [~, aziRowBest] = max(sum(rowLut, 2));
    chanAzi  = find(rowLut(aziRowBest, :));
    A_dbf    = conj(AntCaliData.aziSteerVec(chanAzi, :));    % 信号模型域字典

    % ---------- 找角反(第一优先 isCorner=1)与其实际谱线位置(含步进频耦合) ----------
    isC    = [targets.isCorner] == 1;
    corner = targets(find(isC, 1));
    if ~isfield(prm, 'stepDir') || isempty(prm.stepDir), sd = 1; else, sd = prm.stepDir; end
    dopC   = mod(corner.gate_D + sd*prm.cr*prm.B/(dims.dopBinNum-1) * ...
                 (2*corner.gate_R*3e8/(2*prm.B))/3e8 * dims.dopBinNum, dims.dopBinNum);
    % 候选: 原始点(非 OMP/DML 第二目标), RD 邻域(±1/±2 bin), 取幅度最强
    % (校准遍检测链止于 DBF, 点云尚无 isOMP2nd/isDML2nd 字段 —— 无字段即全为原始点)
    if isfield(tgtOut, 'isOMP2nd'), isO2 = [tgtOut.isOMP2nd] == 1; else, isO2 = false(1, numel(tgtOut)); end
    if isfield(tgtOut, 'isDML2nd'), isD2 = [tgtOut.isDML2nd] == 1; else, isD2 = false(1, numel(tgtOut)); end
    isMother = ~isO2 & ~isD2;                                      % 行向量(勿转置, 防隐性扩展)
    cand = find(isMother & abs([tgtOut.rangIdx] - (corner.gate_R+1)) <= 1 & ...
                              abs([tgtOut.dopplerIdx] - (round(dopC)+1)) <= 2);
    if isempty(cand)
        error('RCS 标定: RD(%d, %d) 附近未找到角反点 —— 检查场景/检测', corner.gate_R+1, round(dopC)+1);
    end
    % 逐候选测幅度(LS 投影, 与 rcsCompute_pipe 同式), 取最强
    aBest = -inf;  kBest = cand(1);
    for k = cand
        x = double(detCfg.doaChanIn(k, chanAzi)).';
        atom = steerFitAt(A_dbf, aziGrid, double(tgtOut(k).azimuth));
        aDbl = 20*log10(abs((atom' * x) / (atom' * atom)) + eps);
        if aDbl > aBest, aBest = aDbl;  kBest = k; end
    end

    % ---------- 噪底曲线 → 平坦段噪底 + 拟合滤形(信号/噪底同口径的已知形状) ----------
    % 平坦段: 距离维 30%~80% 段的噪底中位(通带内, 滤形 ≈ 0)
    gTot   = numel(noiseCfg.floorPow);               % 距离门总数
    floorDb = 10*log10(noiseCfg.floorPow(:) + eps);  % 逐门噪底 (dB)
    flatDb  = median(floorDb(round(0.3*gTot) : round(0.8*gTot)));
    % 拟合段: 噪底明显低于平坦段(滤形可见)且热噪主导(避开近直流泄漏抬底段)的过渡带
    shapeMeas = floorDb - flatDb;                    % 仿真测量滤形 (dB, ≤0)
    gBin = (1:gTot).' - 1;                           % 0 起距离 bin(∝ 差拍频率)
    fitMask = shapeMeas < -0.5 & shapeMeas > -35 & gBin >= 8;   % 干净过渡段(门≥9 且形状可见)
    shapeAt = @(gq) zeros(size(gq));  fitOk = false; fitRms = NaN;
    hasDesign = nargin >= 10 && isfield(filtChain, 'hpf1Fc') && ~isempty(filtChain.hpf1Fc);
    binHz = prm.fs/dims.rngPts;                       % 距离 bin 频宽 (Hz)
    fitFcFreeKHz = [NaN NaN];                         % 曲线自由拟合 fc(核查用)
    if nnz(fitMask) >= 8
        % 两级【一阶】高通级联(芯片 IF 高通物理结构: 每节一阶)。
        % 【2026-11-03 根因修正】近距泄漏闪烁裙边把角反门附近噪底抬高(实测 bin4~20
        %   比平坦段高 +9 dB), 高通过渡带形状在噪底曲线里【不存在】——任何曲线拟合
        %   (含自由阶数/档位锚定)外推深区都失真(实测外推 −2.7 dB 对真实 −11.4 dB,
        %   差值全进 K0 → RCS 偏 +5.6 dB)。唯一可靠形状来源 = 设计/寄存器档位
        %   (芯片配置口径, 非滤波器设计细节): 有 filtChain.hpf1Fc/2Fc 直接用;
        %   曲线自由拟合降级为一致性核查(偏差超 100 kHz 且落在物理档位范围内才警告)。
        mdl = @(th, g) -( 10*log10(1 + (exp(th(1))./g).^2) ...
                        + 10*log10(1 + (exp(th(2))./g).^2) );
        th0 = log([175e3; 500e3]/binHz);              % 初值: 典型档 175k/500k 折 bin
        th  = th0;
        useM = fitMask(:);
        for it = 1:2                                  % 两轮: 拟合 → 剔离群(3×MAD) → 再拟合
            obj = @(t) sum((mdl(t, gBin(useM)) - shapeMeas(useM)).^2);
            opt = optimset('MaxFunEvals', 4000, 'Display', 'off');
            th  = fminsearch(obj, th, opt);
            res = mdl(th, gBin(useM)) - shapeMeas(useM);
            keep = abs(res) < max(3*1.4826*median(abs(res)), 0.5);
            idxM = find(useM);  useM(:) = false;  useM(idxM(keep)) = true;
        end
        fitFcFreeKHz = exp(th)*binHz/1e3;             % 曲线自由拟合结果(核查用)
        if hasDesign
            th = log([filtChain.hpf1Fc; filtChain.hpf2Fc]/binHz);   % 锚定设计档位
        else
            % 无设计档位(旧调用兜底): 自由拟合落到芯片档位范围(87.5~1000 kHz)之外
            % = 退化 → 从典型档位表 {175/235/350/500} kHz 选残差最优一对(精度有限)
            fitFcHz = exp(th)*binHz;
            if any(fitFcHz < 87.5e3) || any(fitFcHz > 1e6)
                gearsKHz = [175 235 350 500];
                bestS = inf;  bestG = [gearsKHz(1); gearsKHz(end)];
                for c1 = gearsKHz
                    for c2 = gearsKHz
                        thG = log([c1; c2]*1e3/binHz);
                        sG  = sum((mdl(thG, gBin(useM)) - shapeMeas(useM)).^2);
                        if sG < bestS, bestS = sG; bestG = [c1; c2]; end
                    end
                end
                th = log(bestG*1e3/binHz);
            end
        end
        shapeAt = @(gq) mdl(th, gq - 1);              % 入参 1 起 bin → 0 起
        fitFcKHz = exp(th)*binHz/1e3;                 % 最终截止频率 (kHz)
        fitRms  = sqrt(mean((mdl(th, gBin(useM)) - shapeMeas(useM)).^2));
        fitOk   = true;
    else
        warning('RCS 标定: 噪底曲线干净过渡段不足(%d 点), 滤形拟合跳过 —— 近距角反标定可能偏低', nnz(fitMask));
    end

    kBestSb = mod(tgtOut(kBest).dopplerIdx - 1, dims.subBins) + 1;   % 角反子带列(nonCohMap/chanSnapshot)
    % ---------- 角反 SNR(信号功率dB − 角反门"平坦底+拟合滤形") ----------
    rGate   = tgtOut(kBest).rangIdx;                       % 角反距离门(1起)
    floorCalDb = flatDb + shapeAt(rGate);                  % 角反门热噪参考(与信号同口径被压)
    % 亚bin窗损补偿(汉宁核同 win2dfft; 与 rcsCompute 目标侧同口径): α 取自单个
    % (距离,多普勒) bin, 峰 bin 幅度带窗叶损 —— 距离维偏移=布设整数 bin(≈0),
    % 多普勒维偏移=步进频耦合公式的小数部分, 二者按归一化窗核精确补到"bin 对中"口径
    compWinDb = [0 0];  fitSrcStr = '理论窗核';
    if nargin >= 13 && strcmpi(winType, 'Hanning') && accumCfg.nonCohMap(tgtOut(kBest).rangIdx, kBestSb) > 0
        % 拟合增量(与目标侧 peakFit 同法, 用角反自身 RD/积累数据): 距离维 5 点峰-插值 ×
        % 多普勒维恢复比; 理论 winScal 对分裂峰欠补(分数角反实测 -0.78dB 系统偏差)
        kc = tgtOut(kBest).rangIdx;  dbc = tgtOut(kBest).dopplerIdx;
        p5 = double(accumCfg.nonCohMap([kc-2:kc+2], kBestSb));
        [~, pj] = max(p5(2:4));  pj = pj + 1;
        [dpkC, PpkC, okR] = hannPk3(p5(pj-1:pj+1), dims.rngPts);  dpkC = dpkC + (pj - 3);  %#ok<NASGU>
        dJ = mod((dbc-3):(dbc+1), dims.dopBinNum) + 1;
        q5 = squeeze(sum(abs(fftCfg.data_sva(kc, :, dJ)).^2, 2));
        [~, qj] = max(q5(2:4));  qj = qj + 1;
        [~, QpkC, okD] = hannPk3(q5(qj-1:qj+1), dims.dopBinNum);
        if okR && okD
            dopRecC = QpkC / q5(3);
            incC = 10*log10(PpkC * dopRecC / double(accumCfg.nonCohMap(kc, kBestSb)));
            snrCalDb = aBest + incC - floorCalDb;        % 角反 SNR (拟合增量, 两频率维)
            fitSrcStr = sprintf('拟合增量 %+.1f dB(两频率维)', incC);
        else
            snrCalDb = aBest - floorCalDb;               % 肩点回退: 不补偿
            fitSrcStr = '无(肩点回退)';
        end
    else
        % 回退: 理论窗核补偿(旧路径)
        if nargin >= 13 && strcmpi(winType, 'Hanning')
            dRc = abs(corner.gate_R - (tgtOut(kBest).rangIdx - 1));
            dDc = abs(mod(dopC + 0.5, 1) - 0.5);
            compWinDb = [winScalDb(dRc, dims.rngPts), winScalDb(dDc, dims.dopBinNum)];
        end
        snrCalDb = aBest + sum(compWinDb) - floorCalDb;
    end

    % ---------- 标定常数(距离按角反布设真值 bin×ΔR) ----------
    Rcal = corner.gate_R * dims.DeltaR;                    % 布设距离
    sigCalDb = corner.rcsDb;                               % 角反真值 RCS (dBsm)
    K0 = sigCalDb - snrCalDb - 40*log10(Rcal);             % 链路常数 (dB)

    % ---------- 相参域角反锚定(2026-11-04): TXRX 相参累积图口径, 独立点 RCS 用 ----------
    % 动因: 非相参图无积累增益(8TX 对齐+8RX 功率平均均不改均值 SNR), 弱目标峰单元读数
    % 噪声主导、RCS 起伏大; 相参图(64 元方位 FFT)+~18dB, 弱目标可读。角反 ele=0 无
    % z 行失配损失, 读数与目标侧(rcsCompute 独立点)同法: 方位峰 bin±1 × 距离 5 点
    % 峰-插值 × 多普勒恢复比 → FFT 归一化/角度维取峰统计抬底/噪底口径全部进 K0Coh 相消
    K0Coh = NaN;  snrCohCal = NaN;  cohFloorDb = NaN;
    if nargin >= 11 && ~isempty(fftCfg) && rGate >= 3 && rGate <= dims.rngUse - 2
        % 多普勒恢复比: 相参快照 = 每 TX 从自己的 chirpPos bin 抓 1 个整数 bin, 8 根码线
        % 各有叶损 —— 逐线 5 点峰-插值求恢复比等权平均(单线拟合只代表一根线: 角反/目标
        % δ 图样不同 → 残差进 K0Coh, 实测 −1.9 dB, 2026-11-04 修)
        if isfield(accumCfg, 'chirpPosMap')
            binsC = double(accumCfg.chirpPosMap(rGate, kBestSb, :));
        else
            binsC = double(tgtOut(kBest).dopplerIdx);        % 兜底: 单 bin(旧口径)
        end
        dopRecC = dopRecMulti(fftCfg, binsC(:), rGate, dims.dopBinNum);
        if isfield(tgtOut, 'cohPeakIdx') && ~isempty(tgtOut(kBest).cohPeakIdx)
            b0c = double(tgtOut(kBest).cohPeakIdx);           % 检测链已存的本点方位峰 bin
        else
            b0c = double(accumCfg.cohPeakPos(rGate, kBestSb));% 兜底: 该单元谱峰位置
        end
        pk3c = zeros(3, 1);  wbc = 0;
        for bb = mod((b0c-2):(b0c), dims.aziFftLen) + 1       % 峰 bin ±1 邻域(回绕)
            wbc = wbc + 1;
            c5 = double(abs(squeeze(accumCfg.aziFftSpec(bb, rGate-2:rGate+2, kBestSb))).^2);
            [~, pjc] = max(c5(2:4));  pjc = pjc + 1;
            [~, pTk, okP] = hannPk3(c5(pjc-1:pjc+1), dims.rngPts);
            if okP, pk3c(wbc) = pTk; else, pk3c(wbc) = c5(3); end
        end
        a3c = sqrt(max(pk3c, 0));                             % 角维幅度抛物线(同 cohReadAt)
        dAc = 2*(a3c(3) - a3c(1)) / (a3c(1) + 2*a3c(2) + a3c(3) + eps);
        CpkC = (a3c(2) - (a3c(1) - a3c(3)) * max(min(dAc, 1), -1) / 4).^2;
        cohFloorDb = 10*log10(double(accumCfg.cohFloorPow));  % 相参域噪底(标定遍, 含取峰统计)
        snrCohCal  = 10*log10(CpkC * dopRecC + eps) - (cohFloorDb + shapeAt(rGate));
        K0Coh = sigCalDb - snrCohCal - 40*log10(Rcal);
        fprintf('RCS 相参锚定   : 角反相参读数 SNR %.1f dB(方位峰bin±1×距离峰插值×多普勒恢复 %.2f) | K0Coh = %.1f dB(独立点相参图口径)\n', ...
            snrCohCal, dopRecC, K0Coh);
    else
        warning('RCS 相参锚定: 角反门数据不足, K0Coh 缺失 —— 独立点 RCS 回退 α 域');
    end

    % ---------- 分段结构(近场/远场; 当前共用角反 K0, 结构已留好) ----------
    segR   = [0, calibInfo.RfarField, Inf];
    K0Seg  = [K0; K0];                                     % 段1(近场)=段2(远场)
    K0CohSeg = [K0Coh; K0Coh];                             % 相参域同分段结构

    % ---------- 天线增益VS角度(暂无数据, 预留接口) ----------
    angGainDb = @(az, el) zeros(size(az));   % [暂缺] 有增益方向图数据后替换此句柄

    rcsCfg = struct('K0Seg', K0Seg, 'segR', segR, 'angGainDb', angGainDb, ...
        'flatDb', flatDb, 'shapeAt', shapeAt, 'K0CohSeg', K0CohSeg, 'cohFloorDb', cohFloorDb, ...
        'corner', struct('idx', kBest, 'RD', [tgtOut(kBest).rangIdx, tgtOut(kBest).dopplerIdx], ...
                         'R', Rcal, 'ampDb', aBest, 'snrDb', snrCalDb, 'snrCohDb', snrCohCal, ...
                         'K0Coh', K0Coh, ...
                         'floorDb', floorCalDb, 'floorRawDb', floorDb(rGate), 'compWinDb', compWinDb, 'sigDb', sigCalDb));

    fprintf('RCS 角反标定   : 角反 RD(%d,%d) 布设 R=%.2f m | SNR %.1f dB (幅度 %.1f − 窗损补 %+.1f/%+.1f dB 噪底参考 %.1f) | 真值 %.0f dBsm | K0 = %.1f dB | %s\n', ...
        rcsCfg.corner.RD(1), rcsCfg.corner.RD(2), Rcal, snrCalDb, aBest, rcsCfg.corner.compWinDb, floorCalDb, sigCalDb, K0, fitSrcStr);
    if fitOk
        if hasDesign, fcSrc = sprintf('设计档位 %.0f/%.0f kHz', fitFcKHz);
        else,         fcSrc = sprintf('曲线拟合 %.0f/%.0f kHz', fitFcKHz); end
        fprintf('RCS 噪底滤形   : 平坦段 %.1f dB, 形状来源=%s, 曲线拟合 fc=%.0f/%.0f kHz(核查), 拟合段 %d 点(RMS %.2f dB), 外推角反门(%d) 滤形 %+.1f dB(本门仿真噪底 %.1f 含泄漏抬底, 不用)\n', ...
            flatDb, fcSrc, fitFcFreeKHz, nnz(useM), fitRms, rGate, shapeAt(rGate), rcsCfg.corner.floorRawDb);
        if hasDesign && all(fitFcFreeKHz >= 87.5 & fitFcFreeKHz <= 1000) ...
                && max(abs(fitFcFreeKHz - fitFcKHz)) > 100
            fprintf('   !!! 滤形核查: 曲线拟合 fc 与设计档位偏差 >100 kHz —— 检查档位配置/滤波链是否变更\n');
        end
    end
    fprintf('RCS 口径       : SNR 比值法(通带目标零补偿; 角反门噪底按噪声曲线拟合滤形重建, 不用滤波器设计)\n');
    fprintf('RCS 分段       : 近场 R<%.2f m / 远场 R≥%.2f m (两段暂用同一 K0)\n', ...
        calibInfo.RfarField, calibInfo.RfarField);

    % ---------- 标定结果图(图19: 噪底曲线 + 滤形拟合 + 角反标定点, 1×2) ----------
    % 左: 逐距离门噪底(仿真测量) + 平坦段电平 + 角反门两个参考点(本门仿真值 vs 重建);
    % 右: 拟合段滤形(仿真测量点 vs 两级高通模型) —— RCS 标定方法学的证据图。
    screenSize = get(0, 'ScreenSize');
    figW = screenSize(3) * 0.98;  figH = screenSize(4) * 0.62;
    figPos = [(screenSize(3)-figW)/2, (screenSize(4)-figH)/2, figW, figH];
    figRcsCal = figure('Position', figPos, 'Name', 'fig19_RCSCalib');
    tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile;   % 左: 噪底曲线
    plot(gBin, floorDb, 'b-', 'LineWidth', 0.8); hold on;
    plot([gBin(1), gBin(end)], flatDb*[1 1], 'k--');
    plot(rGate-1, rcsCfg.corner.floorRawDb, 'rx', 'MarkerSize', 11, 'LineWidth', 2);
    plot(rGate-1, floorCalDb, 'r*', 'MarkerSize', 11);
    grid on; xlabel('距离 bin (0 起)'); ylabel('噪底功率 (dB)');
    legend({'逐门噪底(仿真测量)', sprintf('平坦段 %.1f dB', flatDb), ...
            '角反门仿真值(含泄漏抬底)', '角反门重建参考'}, 'Location', 'best');
    title('逐距离门噪底与角反门参考');

    nexttile;   % 右: 拟合滤形
    plot(gBin(fitMask), shapeMeas(fitMask), 'b.', 'MarkerSize', 10); hold on;
    if fitOk
        plot(gBin, mdl(th, gBin), 'r-', 'LineWidth', 1.2);
        plot(rGate-1, shapeAt(rGate), 'r*', 'MarkerSize', 11);
        legTxt = {'拟合段仿真滤形', '两级高通模型(拟合)', '角反门外推点'};
    else
        legTxt = {'拟合段仿真滤形(点数不足, 未拟合)'};
    end
    xline(rGate-1, 'k--');
    grid on; xlabel('距离 bin (0 起)'); ylabel('相对平坦段 (dB)');
    legend(legTxt, 'Location', 'south');
    title(sprintf('噪底滤形拟合(%d 点, RMS %.2f dB)', nnz(fitMask), fitRms));
    sgtitle(sprintf('RCS 标定(角反 %.0f dBsm, R = %.1f m): SNR = %.1f dB → K0 = %.1f dB', ...
        sigCalDb, Rcal, snrCalDb, K0));

    rcsCfg.figRcsCal = figRcsCal;               % 句柄随结构体返回(main 图序 fig18)
end

function cDb = winScalDb(delta, N)
%WINSCALDB 汉宁窗(对称式, 同 win2dfft 的 hanning())归一化峰响应补偿量(dB >= 0):
%   单音偏移 |delta| bin 时峰 bin 幅度相对 bin 对中的损失; delta 可为向量
w  = hanning(N);
nn = (0:N-1).';
cDb = zeros(size(delta));
for k = 1:numel(delta)
    cDb(k) = -20*log10(abs(w.'*exp(-1i*2*pi*delta(k)*nn/N)) / sum(w));
end
end

function dR = dopRecMulti(fftCfg, binsGrab, k, dopBinNum)
%DOPRECMULTI 8 根 DDMA 码线各自被抓 bin 的多普勒叶损恢复比(等权平均):
%   相参快照 = 每 TX 从自己的 chirpPos bin 抓 1 个整数 bin, 8 根线各有叶损;
%   单线拟合只代表一根线, 角反/目标 δ 图样不同 → 残差进 K0Coh。
%   逐线 5 点峰-插值 → rec_t = 峰功率/被抓bin功率, 肩点该线回退 1
dR = 1;
if isempty(binsGrab), return; end
recs = ones(numel(binsGrab), 1);  okAny = false;
for tt = 1:numel(binsGrab)
    db = double(binsGrab(tt));
    dJ = mod((db-3):(db+1), dopBinNum) + 1;            % 5 点全 bin(1 起, 循环回绕)
    q5 = squeeze(sum(abs(fftCfg.data_sva(k, :, dJ)).^2, 2));
    [~, qj] = max(q5(2:4));  qj = qj + 1;
    [~, Qpk, okD] = hannPk3(q5(qj-1:qj+1), dopBinNum);
    if okD, recs(tt) = Qpk / q5(3);  okAny = true; end
end
if okAny, dR = mean(recs); end
end

function [d, pk, ok] = hannPk3(p, N)
%HANNPK3 3 点功率峰的汉宁精确恢复(同 rcsCompute_pipe, 2026-11-04 替代纯抛物线):
%   位置 = 幅度域 Grandke 插值(抛物线位置在 δ→0.5 偏短 ~0.25 bin);
%   峰功 = p0×10^(winScalDb(δ,N)/10)(抛物线峰功在 δ=0.5 低估 1.18dB)
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

