function [tgtOut, rcsCfg, figRcs] = rcsCompute_pipe(tgtOut, detCfg, AntCaliData, rowLut, dims, prm, targets, rcsCfg, winType)
%RCSCOMPUTE_PIPE RCS 估计: 对全部点云逐点算 RCS, 写回点云并绘图(SNR 口径)
%   [tgtOut, rcsCfg, figRcs] = rcsCompute_pipe(tgtOut, detCfg, AntCaliData, rowLut, ...
%                                              dims, prm, targets, rcsCfg)
%   逐点幅度来源(均为"目标复幅度 α"同一量纲, 积累快照单位):
%     - 原始点(含 DML 换过方位的第一目标): 在其当前方位上做表格原子 LS 投影 α = aᴴx/(aᴴa);
%     - OMP 第二目标: 用 OMP 联合 LS 权重(ompW2Db);
%     - DML 第二目标: 用 DML 联合 LS 权重(dmlW2Db)。
%   RCS 公式(SNR 比值, 分段, 用角反标定常数):
%     σ_dB = SNR_dB + 40log10(R) + K0Seg(段) + [G(0°) − G(θ)]
%     SNR_dB = 2·α_dB − 10log10(floorPow(该距离门))   (信号功率 |α|² / 该门噪底)
%     信号与噪底同 bin 同滤形 → SNR 与 IF 滤波器(高通/低通/纹波)形状无关, 无需任何
%     滤形补偿(算法工程师口径: 滤形未知); 分段按 R 与 segR 比较(近场/远场);
%     增益修正项暂无方向图数据 = 0。
%   写回字段: rcsDb(估计 dBsm) / rcsSeg(1 近场, 2 远场) / rcsAmpDb(幅度 dB) / rcsSnrDb(SNR dB)。
%   绘图(图35, 2×2): RCS-距离 / RCS-方位 / RCS 直方图 / x-y 位置(颜色=RCS);
%   并打印逐真值目标对比表(场景 targets 的 rcsDb 真值 vs 点云峰值估计)。

    % ---------- 方位字典(最长行 36 通道, 共轭域) ----------
    aziRange = AntCaliData.aziSteerVecRange;
    aziGrid  = aziRange(1) : aziRange(3) : aziRange(2);
    [~, aziRowBest] = max(sum(rowLut, 2));
    chanAzi  = find(rowLut(aziRowBest, :));
    A_dbf    = conj(AntCaliData.aziSteerVec(chanAzi, :));    % 信号模型域字典

    % ---------- 亚bin窗损补偿(汉宁核同 win2dfft; SVA 自适应无固定核不补) ----------
    % 距离维: 逐点用亚bin插值偏移精确补; 多普勒维: DDMA 16 槽后亚bin相位不可逐点知,
    % 按目标速度均匀分布取窗响应均值(期望口径) —— 角反侧(rcsCalibrate)用其确定性
    % 耦合位置精确补, 两侧同折算到"bin 对中"口径, RCS 不再随 bin 位置抖动
    doWinComp = nargin >= 9 && strcmpi(winType, 'Hanning');
    if doWinComp
        % 多普勒维常数: 16 相位全链扫描标定(2026-11-03, gate_R=285 整数 bin 单目标;
        % 理论均值 +1.75 dB 下偏差均值 +2.12 / std 0.63 -> 净共模 = −0.37 dB)。
        % 差异来源: 实际幅度是"峰分裂邻点 max 选择"读数(赢家常为 0.5~0.6 bin 偏移的
        % 邻点), 与固定 bin 窗读不是同一口径, 无闭式 —— 照 winStat 口径收实测标定常数;
        % 残差散布 ±0.6 dB 为选择统计固有, 如实保留
        compDopDb = -0.37;
        fprintf(['RCS 窗损补偿   : 距离维逐点(亚bin插值偏移) + 多普勒维扫描标定 %+.2f dB(16 相位全链标定, 残差散布 ±0.6 dB)' newline], compDopDb);
    else
        compDopDb = 0;
        fprintf(['RCS 窗损补偿   : %s 无固定核不补偿(口径如实)' newline], 'SVA 自适应窗');
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

    for t = 1:nTgt
        if tgtOut(t).isOMP2nd == 1                 % OMP 第二目标: 用 OMP 联合 LS 权重
            aDbl = tgtOut(t).ompW2Db;
        elseif tgtOut(t).isDML2nd == 1             % DML 第二目标: 用 DML 联合 LS 权重
            aDbl = tgtOut(t).dmlW2Db;
        elseif tgtOut(t).ompHas2nd && ~isnan(tgtOut(t).ompW1Db)
            aDbl = tgtOut(t).ompW1Db;              % OMP 对的母点: 联合第一分量(去邻目标泄漏, 与第二目标对称)
        elseif tgtOut(t).dmlPass && ~isnan(tgtOut(t).dmlW1Db)
            aDbl = tgtOut(t).dmlW1Db;              % DML 对的母点: 联合第一分量(ampW(1), dmlRefine 已存)
        else                                       % 原始点: 当前方位的 LS 投影幅度
            x = double(detCfg.doaChanIn(t, chanAzi)).';
            atom = steerFitAt(A_dbf, aziGrid, double(tgtOut(t).azimuth));
            aDbl = 20*log10(abs((atom' * x) / (atom' * atom)) + eps);
        end
        if doWinComp                                     % 亚bin窗损: 距离逐点 + 多普勒均值
            dRt = abs((tgtOut(t).rangIdx - 1) - double(tgtOut(t).range)/dims.DeltaR);
            aDbl = aDbl + winScalDb(dRt, dims.rngPts) + compDopDb;
        end

        R = double(tgtOut(t).range);               % 插值距离 (m)
        seg = find(R >= rcsCfg.segR(1:end-1) & R < rcsCfg.segR(2:end), 1);
        if isempty(seg), seg = numel(rcsCfg.K0Seg); end
        az = double(tgtOut(t).azimuth);  el = double(tgtOut(t).elevation);
        % SNR 口径: 信号功率dB(= α_dB, 因 20log|α| ≡ 10log|α|²) − 噪底参考(平坦底+拟合滤形)。
        % 通带目标滤形 ≈ 0 → 纯 SNR 零补偿; 近距点(如角反)信号被压, 噪底参考同口径被压,
        % 滤形来自噪声曲线拟合(不用滤波器设计) → SNR 仍与滤形无关
        snrDb(t) = aDbl - (rcsCfg.flatDb + rcsCfg.shapeAt(tgtOut(t).rangIdx));
        % σ_dB = SNR + 40log10(R) + 段常数 + 增益VS角度修正(暂 0; 标定角反在 0°)
        rcsDb(t)  = snrDb(t) + 40*log10(R) + rcsCfg.K0Seg(seg) + ...
                    rcsCfg.angGainDb(az, el);      % 注: G(0°)−G(θ) 已含负号方向, 有数据后核对符号
        rcsSeg(t) = seg;
        ampDb(t)  = aDbl;
    end

    for t = 1:nTgt
        tgtOut(t).rcsDb    = rcsDb(t);
        tgtOut(t).rcsSeg   = rcsSeg(t);
        tgtOut(t).rcsAmpDb = ampDb(t);
        tgtOut(t).rcsSnrDb = snrDb(t);
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
        % 窗内可能含同 RD 的邻近角目标(如 3°/4.5° 对, 间隔仅 1.5°)与弱旁瓣/鬼影点:
        % 先把方位门限收紧到真值 ±1° 再取最强(窄门能分开 1.5° 对; 强者优先排除
        % 恰好压在真值方位上的弱鬼影), 窄门无点才回退 ±4° 宽门
        idxN = idx(abs([tgtOut(idx).azimuth] - targets(k).azi) <= 1.0 & ~claimedRcs(idx));
        if isempty(idxN)                          % 无未认领点: 回退全部窄门点(含已认领)
            idxN = idx(abs([tgtOut(idx).azimuth] - targets(k).azi) <= 1.0);
        end
        if ~isempty(idxN), idx = idxN; end
        % 一对一最近角分配(配对两目标角距可小于 ±1° 窗, 各取最强会选中同一点):
        % 先按 |方位−真值| 最近, 平手再比 rcsDb; 认领后不再给下一个真值
        [~, im] = max(rcsDb(idx));                       % 窗内未认领者取最强(认领表防配对两点抢同一点)
        claimedRcs(idx(im)) = true;
        pkDb = rcsDb(idx(im));
        confTxt = '低置信';  if tgtOut(idx(im)).rcsConf == 1, confTxt = '高置信'; end
        fprintf('   目标%d: R=%6.1f m 真值 %5.1f dBsm → 峰值估计 %6.1f dBsm (n=%d, 偏差 %+5.1f dB, %s)\n', ...
            k, targets(k).gate_R*3e8/(2*prm.B), targets(k).rcsDb, pkDb, numel(idx), pkDb - targets(k).rcsDb, confTxt);
    end

    % ---------- 绘图(图23, 2×2) ----------
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
    sgtitle(sprintf('RCS 估计 (角反 %.0f dBsm 标定, K0=%.1f dB): 近场 %d 点 / 远场 %d 点', ...
        rcsCfg.corner.sigDb, rcsCfg.K0Seg(1), sum(rcsSeg == 1), sum(rcsSeg == 2)));
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
