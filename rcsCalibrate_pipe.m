function rcsCfg = rcsCalibrate_pipe(tgtOut, detCfg, AntCaliData, rowLut, dims, prm, targets, calibInfo, noiseCfg)
%RCSCALIBRATE_PIPE RCS 标定: 用角反(targets 中 isCorner=1)标定链路常数(SNR 比值法)
%   rcsCfg = rcsCalibrate_pipe(tgtOut, detCfg, AntCaliData, rowLut, dims, prm, targets, calibInfo, noiseCfg)
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
    if nnz(fitMask) >= 8
        % 两级高通级联模型(工程通用先验): shape(g) = −10·[log10(1+(gc1/g)^(2n1)) + log10(1+(gc2/g)^(2n2))]
        mdl = @(th, g) -( 10*log10(1 + (exp(th(1))./g).^(2*exp(th(2)))) ...
                        + 10*log10(1 + (exp(th(3))./g).^(2*exp(th(4)))) );
        th0 = log([8 1; 31 2]);  th = th0(:);        % 初值: 两级高通的典型量级(拟合定标)
        useM = fitMask(:);
        for it = 1:2                                  % 两轮: 拟合 → 剔离群(3×MAD) → 再拟合
            obj = @(t) sum((mdl(t, gBin(useM)) - shapeMeas(useM)).^2);
            % 4 参数默认预算 800 次常不够收敛(退出时打"超过最大数目"提示),
            % 放宽到 4000 并静默 —— 目标函数只有百来个点, 多算不耗时
            opt = optimset('MaxFunEvals', 4000, 'Display', 'off');
            th  = fminsearch(obj, th, opt);
            res = mdl(th, gBin(useM)) - shapeMeas(useM);
            keep = abs(res) < max(3*1.4826*median(abs(res)), 0.5);
            idxM = find(useM);  useM(:) = false;  useM(idxM(keep)) = true;
        end
        shapeAt = @(gq) mdl(th, gq - 1);              % 入参 1 起 bin → 0 起
        fitRms  = sqrt(mean((mdl(th, gBin(useM)) - shapeMeas(useM)).^2));
        fitOk   = true;
    else
        warning('RCS 标定: 噪底曲线干净过渡段不足(%d 点), 滤形拟合跳过 —— 近距角反标定可能偏低', nnz(fitMask));
    end

    % ---------- 角反 SNR(信号功率dB − 角反门"平坦底+拟合滤形") ----------
    rGate   = tgtOut(kBest).rangIdx;                       % 角反距离门(1起)
    floorCalDb = flatDb + shapeAt(rGate);                  % 角反门热噪参考(与信号同口径被压)
    snrCalDb = aBest - floorCalDb;                         % 角反 SNR (dB)

    % ---------- 标定常数(距离按角反布设真值 bin×ΔR) ----------
    Rcal = corner.gate_R * dims.DeltaR;                    % 布设距离
    sigCalDb = corner.rcsDb;                               % 角反真值 RCS (dBsm)
    K0 = sigCalDb - snrCalDb - 40*log10(Rcal);             % 链路常数 (dB)

    % ---------- 分段结构(近场/远场; 当前共用角反 K0, 结构已留好) ----------
    segR   = [0, calibInfo.RfarField, Inf];
    K0Seg  = [K0; K0];                                     % 段1(近场)=段2(远场)

    % ---------- 天线增益VS角度(暂无数据, 预留接口) ----------
    angGainDb = @(az, el) zeros(size(az));   % [暂缺] 有增益方向图数据后替换此句柄

    rcsCfg = struct('K0Seg', K0Seg, 'segR', segR, 'angGainDb', angGainDb, ...
        'flatDb', flatDb, 'shapeAt', shapeAt, ...
        'corner', struct('idx', kBest, 'RD', [tgtOut(kBest).rangIdx, tgtOut(kBest).dopplerIdx], ...
                         'R', Rcal, 'ampDb', aBest, 'snrDb', snrCalDb, ...
                         'floorDb', floorCalDb, 'floorRawDb', floorDb(rGate), 'sigDb', sigCalDb));

    fprintf('RCS 角反标定   : 角反 RD(%d,%d) 布设 R=%.2f m | SNR %.1f dB (幅度 %.1f − 噪底参考 %.1f) | 真值 %.0f dBsm | K0 = %.1f dB\n', ...
        rcsCfg.corner.RD(1), rcsCfg.corner.RD(2), Rcal, snrCalDb, aBest, floorCalDb, sigCalDb, K0);
    if fitOk
        fprintf('RCS 噪底滤形   : 平坦段 %.1f dB, 过渡段拟合 %d 点(RMS %.2f dB), 外推角反门(%d) 滤形 %+.1f dB(本门仿真噪底 %.1f 含泄漏抬底, 不用)\n', ...
            flatDb, nnz(useM), fitRms, rGate, shapeAt(rGate), rcsCfg.corner.floorRawDb);
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
