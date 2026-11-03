function fitRep = peakFit_pipe(tgtOut, accumCfg, noiseCfg, fftCfg, rcsCfg, AntCaliData, rowLut, dims, prm, winType, targets, anaFun)
%PEAKFIT_PIPE 峰值拟合对比模块(独立解耦: 只读现成数据, 不改任何上游模块/tgtOut)
%   fitRep = peakFit_pipe(tgtOut, accumCfg, noiseCfg, fftCfg, rcsCfg, AntCaliData, ...
%                         rowLut, dims, prm, winType, targets, anaFun)
%   分流: tgtOut.peak_type ∈ {3,4} = 峰值点(做拟合); 0 = 延伸点(不拟合, 保留原始读数)。
%   拟合对象 = 距离 bin(中频频率) + 多普勒 bin(多普勒频率)两个频率维(不直接拟距离/速度,
%   bin→距离/速度的精确转换含步进频解耦, 另行处理; 本模块 SNR/RCS 只做比值恢复不碰转换):
%     距离 bin    nonCohMap 距离维 5 点, 峰-插值两步(先最强 cell 再 3 点抛物线)
%     多普勒 bin  fftCfg.data_sva 全分辨率 5 点, 同款两步(RX 功率和) → 恢复比 dopRec
%     ncSNR       距离维拟合峰功率 × dopRec / floorPow(两频率维都拟)
%     cohSNR      cohPeakPow 距离维拟合 × dopRec / cohFloorPow
%     RCS     复用 ncSNR 拟合增量(同方法同源: rcsDb + (ncSnrFit − nonCohSNR),
%             增量只含距离维窗损恢复, 多普勒项相消; 不再自建通道域估计)
%   打印: 逐真值目标"真值 / 原始 / 拟合"三方对比(ncSNR/cohSNR 另附 anaFun 解析值);
%   fitRep 逐点返回拟合值, 供主函数后续决定是否采纳写回(本模块不改 tgtOut)。


nT = numel(tgtOut);
fitRep = struct('isPeak', false(nT,1), 'rngFit', nan(nT,1), 'spdFit', nan(nT,1), ...
                'ncSnrFit', nan(nT,1), 'cohSnrFit', nan(nT,1), 'rcsFit', nan(nT,1));
for t = 1:nT
    fitRep.isPeak(t) = any(tgtOut(t).peak_type == [3 4]);
    if ~fitRep.isPeak(t), continue; end
    k  = tgtOut(t).rangIdx;                     % 距离 bin(1 起)
    db = tgtOut(t).dopplerIdx;                  % 解模糊后全 bin(1 起)
    sb = mod(db - 1, dims.subBins) + 1;         % 子带列(nonCohMap/cohPeakPow/chanSnapshot 第 2/3 维)
    if k < 3 || k > dims.rngUse - 2, continue; end   % 边界 bin 无邻域, 不拟合
    % 拟合对象 = 距离 bin(中频频率)与多普勒 bin(多普勒频率)两个频率维 —— 不直接拟
    % 距离/速度: bin→距离/速度的精确转换含步进频解耦, 另行处理(本模块只做比值恢复)
    % ① 距离 bin: nonCohMap 距离维 5 点, 先取内三点最强 cell 再绕它 3 点抛物线
    %    (峰-插值两步: 谱线真峰常不在检测 bin, 直接以检测 bin 为中心的抛物线会遇单调段失效)
    p5 = double(accumCfg.nonCohMap([k-2:k+2], sb));
    [~, pj] = max(p5(2:4));  pj = pj + 1;                   % 内三点(2..4)最强位置
    [dpk, Ppk, okR] = hannPk3(p5(pj-1:pj+1), dims.rngPts);
    dpk = dpk + (pj - 3);                                   % 相对检测 bin 总偏移(可超 ±0.5)
    fitRep.rngFit(t) = (k - 1 + dpk) * dims.DeltaR;         % 仅显示用粗转换
    % ② 多普勒 bin: fftCfg 全分辨率 5 点, 同款峰-插值两步(RX 功率和; 循环回绕)。
    %    累积链每个 TX 只抓 1 个整数 bin(jointAccum 的 chirpPos), 谱线落 bin 间时单 bin
    %    读数带窗叶损(本实现 δ=0.76 → −3.4dB) —— 多普勒恢复比 = 拟合峰功率/被抓 bin 功率
    dJ = mod((db-3):(db+1), dims.dopBinNum) + 1;            % 5 点全 bin(1 起)
    q5 = squeeze(sum(abs(fftCfg.data_sva(k, :, dJ)).^2, 2)); % [5×1] 功率
    [~, qj] = max(q5(2:4));  qj = qj + 1;
    [ddk, Qpk, okD] = hannPk3(q5(qj-1:qj+1), dims.dopBinNum);
    ddk   = ddk + (qj - 3);
    dopRec = Qpk / q5(3);                                   % 多普勒维恢复比(两域同窗, 比值域无关)
    fitRep.spdFit(t) = (db - 1 + ddk) * dims.Vres;          % 谱线速度粗显示; bin→速度精确解耦后做
    if t == 1   % DBG3: 首点五点分布
        fprintf(['    [DBG3] db=%d q5bins=%s q5dB=%s qj=%d ddk=%+.2f dopRec=%.2f' newline], ...
            db, mat2str(mod((db-3):(db+1), dims.dopBinNum)), ...
            mat2str(round(10*log10(q5 + eps), 1)), qj, ddk, dopRec);
    end
    % ③ cohSNR 距离维: cohPeakPow 5 点, 同款峰-插值两步
    c5 = double(accumCfg.cohPeakPow([k-2:k+2], sb));
    [~, cj] = max(c5(2:4));  cj = cj + 1;
    [~, Cpk, okC] = hannPk3(c5(cj-1:cj+1), dims.rngPts);
    if ~okR || ~okD || ~okC
        % 肩点(抛物线前提不成立): 全部回落原始读数, 不拟合
        fitRep.rngFit(t) = tgtOut(t).range;   fitRep.spdFit(t) = tgtOut(t).speed;
        fitRep.ncSnrFit(t) = tgtOut(t).nonCohSNR;  fitRep.cohSnrFit(t) = tgtOut(t).cohSNR;
        fitRep.rcsFit(t) = tgtOut(t).rcsDb;
        continue
    end
    % ④ ncSNR/cohSNR: 距离维拟合峰功率 × 多普勒维恢复比(两频率维都拟)
    fitRep.ncSnrFit(t)  = 10*log10(Ppk * dopRec / noiseCfg.floorPow(k));
    fitRep.cohSnrFit(t) = 10*log10(Cpk * dopRec / accumCfg.cohFloorPow);
    % ⑤ RCS: 生产口径已在 rcsCompute 双域化(2026-11-04: 独立点=TXRX 相参图读数+K0Coh,
    %    配对点=α域联合权重) —— 本模块为解耦对比, RCS 行直接引用生产值(旧 ncSNR 拟合
    %    增量与非相参域绑死, 对相参域生产值会跨域混叠, 已弃)
    fitRep.rcsFit(t) = tgtOut(t).rcsDb;
end

% ---- 逐真值目标三方对比打印(代表点=簇内窄方位门最强, 与 rcsCompute 同口径) ----
fprintf(['峰值拟合对比    : 峰值点 %d / 延伸点 %d(不拟合, 保留原始读数)' newline], sum(fitRep.isPeak), sum(~fitRep.isPeak));
if ~isfield(prm, 'stepDir') || isempty(prm.stepDir), sd = 1; else, sd = prm.stepDir; end
claimedFit = false(1, numel(tgtOut));   % 一对一分配(同 rcsCompute 口径): 配对两目标不抢同一点
claimedR   = false(1, numel(tgtOut));   % RCS 代表的独立认领表(按最近角+拟合值)
for g = 1:numel(targets)
    dopC = mod(targets(g).gate_D + sd*prm.cr*prm.B/(dims.dopBinNum-1) * ...
               (2*targets(g).gate_R*3e8/(2*prm.B))/3e8 * dims.dopBinNum, dims.dopBinNum);
    idx = find(abs([tgtOut.rangIdx] - (targets(g).gate_R+1)) <= 2 & ...
               abs([tgtOut.dopplerIdx] - (round(dopC)+1)) <= 3 & ...
               abs([tgtOut.azimuth] - targets(g).azi) <= 4);
    if isempty(idx), continue; end
    idxN = idx(abs([tgtOut(idx).azimuth] - targets(g).azi) <= 1);
    if ~isempty(idxN), idx = idxN; end
    % 代表点: 与链路预算审计同判据(最近期望单元, 距离差+循环多普勒差加权) —— 两处打印同点,
    % 原始/拟合数值可直接对照(按拟合值选代表会选到簇内另一点, 与审计行对不上)
    dAbsT = abs([tgtOut(idx).rangIdx] - (targets(g).gate_R + 1));
    dCirT = mod(abs([tgtOut(idx).dopplerIdx] - (round(dopC)+1)), dims.dopBinNum);
    dCirT = min(dCirT, dims.dopBinNum - dCirT);
    okC = ~claimedFit(idx);
    if ~any(okC), okC = true(size(idx)); end        % 全被认领则回退
    dTmp = dAbsT + dCirT/dims.subBins;  dTmp(~okC) = inf;
    [~, im] = min(dTmp);
    t  = idx(im);  claimedFit(t) = true;
    an = anaFun(targets(g).gate_R*dims.DeltaR, targets(g).rcsDb, targets(g).ele);
    splitNames = {'独立点','OMP点1','OMP点2','DML点1','DML点2'};   % 与 main_sim splitTag 同源
    fprintf(['  目标%d @真值 R=%.2f m, v=%+.2f m/s, RCS=%.1f dBsm | 代表点: %s, 峰值点=%d, 簇内 %d 点' newline], ...
        g, targets(g).gate_R*dims.DeltaR, targets(g).gate_D*dims.Vres, targets(g).rcsDb, ...
        splitTagStr(tgtOut(t)), fitRep.isPeak(t), numel(idx));
    fprintf(['    距离   (m) : 真值 %8.3f | 原始 %8.3f | 拟合 %8.3f  (原始→拟合偏差 %+7.3f → %+7.3f)' newline], ...
        targets(g).gate_R*dims.DeltaR, tgtOut(t).range, fitRep.rngFit(t), ...
        tgtOut(t).range - targets(g).gate_R*dims.DeltaR, fitRep.rngFit(t) - targets(g).gate_R*dims.DeltaR);
    fprintf(['    速度 (m/s) : 线真值 %+8.3f | 原始 %+8.3f | 拟合 %+8.3f   (点云速度=耦合谱线速度, 解模糊前口径)' newline], ...
        dopC*dims.Vres, tgtOut(t).speed, fitRep.spdFit(t));
    % RCS 行代表点 = 簇内最强峰点(与 rcsCompute 逐目标对比同口径; 审计匹配点可能选到
    % 分裂簇的较弱 cell, 低估 RCS) —— SNR 两行仍用审计同判据(与审计行数值对齐)
    okR = ~claimedR(idx);
    if ~any(okR), okR = true(size(idx)); end
    rTmp = fitRep.rcsFit(idx);  rTmp(~okR) = -inf;      % 窗内未认领者取最强(纯最近角会选中贴角弱肩点)
    [~, iR] = max(rTmp);
    tR = idx(iR);  claimedR(tR) = true;
    if tgtOut(tR).rcsDom == 1, domR = '相参图'; else, domR = 'α域'; end
    fprintf(['    RCS (dBsm) : 真值 %8.2f | 生产 %8.2f (%s, 代表点=簇内最强峰[%s], 非审计匹配点; 无独立解析参考——雷达方程往返=输入回显)' newline], ...
        targets(g).rcsDb, tgtOut(tR).rcsDb, domR, splitTagStr(tgtOut(tR)));
    fprintf(['    ncSNR (dB) : 原始 %6.2f | 拟合 %6.2f | 解析 %6.2f' newline], ...
        tgtOut(t).nonCohSNR, fitRep.ncSnrFit(t), an.SNR_dB);
    fprintf(['    cohSNR (dB): 原始 %6.2f | 拟合 %6.2f | 解析 %6.2f' newline], ...
        tgtOut(t).cohSNR, fitRep.cohSnrFit(t), an.SNRcoh_dB);
end
end

function [d, pk, ok] = hannPk3(p, N)
%HANNPK3 3 点功率峰的汉宁精确恢复(2026-11-04 替代纯抛物线, 同 rcsCompute_pipe):
%   位置 = 幅度域 Grandke 插值 δ=2(a₃−a₁)/(a₁+2a₂+a₃)(抛物线位置在 δ→0.5 偏短
%   ~0.25 bin); 峰功 = p0×10^(winScalDb(δ,N)/10)(抛物线峰功在 δ=0.5 低估 1.18dB
%   —— 分数 bin 目标系统性偏低, 默认场景 285.5 bin 即中招)。肩点判据同旧 para3p
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

function s = splitTagStr(p)
%SPLITTAGSTR 点角色标签(防御式): 空/缺失/非标量 → '未知(补查)' —— 若打印出未知, 说明
%   splitTag 链路有漏打点(OMP 初始化+配对, DML 补打), 需回两模块核查
n = {'独立点','OMP点1','OMP点2','DML点1','DML点2'};
if ~isfield(p, 'splitTag') || ~isscalar(p.splitTag) || isempty(p.splitTag)
    s = '未知(补查)';  return
end
s = n{p.splitTag + 1};
end
