function [cCal, calibInfo] = calib0deg(fftCfg, encodeSeq, nSamples, nChirps, prm, R1, v1, antenna_t, antenna_r, lambda)
%CALIB0DEG 0° 校准(含近场补偿): 用已知目标1 提取 64 通道相位 → 1×64 校准表
%   [cCal, calibInfo] = calib0deg(fftCfg, encodeSeq, nSamples, nChirps, prm, ...
%                                 R1, v1, antenna_t, antenna_r, lambda)
%
%   原理:
%     目标1 在 0°(方位/俯仰 0)已知距离 R1 处。理想平面波下 64 通道相位全 0,
%     仿真测量相位 = 通道失配 + 【近场球面波相位】(目标不在远场时)。
%     近场补偿: 目标在 y 轴(正视方向)上, TX i / RX j 到目标的单程距离
%       Rt_i = sqrt(x_ti^2 + R1^2 + z_ti^2),  Rr_j = sqrt(x_rj^2 + R1^2 + z_rj^2)
%     近场引入的额外相位(相对远场平面波 2R1) = -k*(Rt_i + Rr_j - 2R1), k = 2π/λ。
%     从仿真测量相位中把它加回去(等效除掉 e^{-jk ΔR}), 剩下的才是通道失配 →
%     校准表 = 共轭相位(单位模, 只校相位不校幅度)。
%     每个TX的目标1 谱线位置: 距离 bin = 差拍/bin间隔(全 TX 共用);
%     速度 bin = 多普勒 + DDMA 码位×(nChirps/16) + 步进频耦合 f_step·τ·nChirps。
%
%   校准检查图(单图): 校准源目标1(可能带近场效应) × 校准表 → 全阵 64 元 +
%   最长行 角度谱, 与等幅导向 0° 理想阵因子对比(虚线) —— 校准正确时主峰在 0°。
%
%   输入:
%     fftCfg.data_sva, encodeSeq, nSamples, nChirps, prm(fs/B/Ts/Tc/cr/fStart),
%     R1/v1(目标1 已知距离/速度), antenna_t/antenna_r [3×8] 天线相位中心(居中后, m)
%   输出:
%     cCal [1×64] 复数校准表, |cCal|=1; 列顺序同 A_manifold(TX 内循环/RX 外循环)
%     calibInfo: .m .rBin .dopBins .phRaw .phCode .phErr(近场补偿后的通道相位)
%                .phNF [nTX×nRX] 被扣除的近场相位(rad)
%                .RfarField  可不做近场补偿的最小距离(m, 双程 ≤ π/8 判据)
%                .figCal     校准检查图句柄(两条角度谱 + 理想参考)
    slope = prm.B / prm.Ts;                    % 调频斜率
    dfBin = prm.fs / nSamples;                 % 距离 bin 间隔 (Hz)
    if ~isfield(prm, 'stepDir') || isempty(prm.stepDir), prm.stepDir = 1; end
    fStep = prm.stepDir * prm.cr * prm.B / (nChirps - 1);    % 步进频 chirp 间起始频移(带符号: 下跳频为负)
    k = 2*pi/lambda;                           % 波数(lambda 由主函数算好传入)

    nGroups = numel(encodeSeq);
    slots = find(encodeSeq) - 1;               % 各 TX 的码槽位(0 起)
    nTX = numel(slots);
    nRX = size(fftCfg.data_sva, 2);
    win = hanning(nChirps);

    % ================= 目标通道提取(通用: 已知 R/v → 各 TX bin → 邻域搜索) =================
    extract = @(R, v) extractCh(fftCfg, slots, nGroups, nChirps, win, ...
                                slope, dfBin, fStep, lambda, prm.Tc, R, v);

    [m, rBin, dopFound, phCode] = extract(R1, v1);

    % ---- 近场补偿: 扣除已知 R1 的球面波相位(目标在正视 y 轴上) ----
    phNF = nearFieldPh(antenna_t, antenna_r, R1, k);

    % ---- 通道相位: 近场补偿 → 去均值 → 共轭 → 单位模 ----
    % 【不扣解析码相位 phCode】: 仿真验证(理想链路)提取相位 = 每 TX 常量 + 每 RX
    % 常量(行列可分, 不可分离残差 0°), 即 DDMA 码/窗相位常量与 TX 通道相位【不可
    % 分辨】且同属校准表应覆盖的内容 → 一并留在表内。注意: 窗相位含 πk(N-1)/N
    % 项, 在奇偶不同的多普勒 bin 上相差 180°——用本表校准其它 bin 的通道时,
    % 若 (k'-k) 为奇数需补 180°。
    phRaw = angle(m);                            % [nTX×nRX]
    phErr = phRaw - phNF;                        % 扣近场球面波相位(距离相关项)
    phErr = phErr - mean(phErr(:));              % 去公共相位
    cCal  = reshape(exp(-1i*phErr), 1, []);      % [1×64], TX 内循环/RX 外循环(同 A_manifold)

    % ---- 可不做近场补偿的最小距离: 双程最坏 (a_t²+a_r²)/(2R)≤D²/(4R), 相位≤π/8 ----
    Dx = max([antenna_t(1,:) antenna_r(1,:)]) - min([antenna_t(1,:) antenna_r(1,:)]);
    Dz = max([antenna_t(3,:) antenna_r(3,:)]) - min([antenna_t(3,:) antenna_r(3,:)]);
    D  = hypot(Dx, Dz);                          % 整阵对角尺寸
    RfarField = 4*D^2/lambda;

    % ================= 校准检查图(单图: 校准源 × 校准表 角度谱) =================
    theta = linspace(-90, 90, 1801) * pi/180;
    xv = antenna_t(1,:).' + antenna_r(1,:);      % [nTX×nRX] 虚拟阵元 x
    zv = antenna_t(3,:).' + antenna_r(3,:);
    zRows = unique(zv(:));
    nInRow = zeros(numel(zRows), 1);
    for r = 1:numel(zRows), nInRow(r) = sum(zv(:) == zRows(r)); end
    [~, rMax] = max(nInRow);
    rowSel = (zv(:) == zRows(rMax));

    % ---- 幅度归一化(不做幅度校准, 与雷达处理一致): 通道取单位模, 只保留相位 ----
    % 说明: 波束形成若直接用通道幅度 |m| 加权, 各通道幅度差(TX_A/RX_A 设计差异等)
    % 会成为无设计的随机锥削 → 旁瓣抬升、方向图畸变。
    % 雷达的口径是【幅度归一化】= 通道单位模 × 相位校准表(表本身 |c|=1)。
    mCal1 = exp(1i*angle(m .* exp(-1i*phNF))) .* reshape(cCal, nTX, nRX);   % 目标1: 归一化+校准

    % ---- 理想参考(等幅导向 0° 的阵因子, 虚线): 全阵与最长行 ----
    BIdA = abs(exp(1i*k*xv(:)*sin(theta)).' * ones(numel(xv), 1));       % 理想 全阵
    BIdB = abs(exp(1i*k*xv(rowSel)*sin(theta)).' * ones(nnz(rowSel), 1)); % 理想 最长行

    B1a = abs(exp(1i*k*xv(:)*sin(theta)).' * mCal1(:));         % 全阵 64 元
    B1b = abs(exp(1i*k*xv(rowSel)*sin(theta)).' * mCal1(rowSel)); % 最长行

    screenSize = get(0, 'ScreenSize');
    figW = screenSize(3) * 0.98;  figH = screenSize(4) * 0.55;
    figPos = [(screenSize(3)-figW)/2, (screenSize(4)-figH)/2, figW, figH];
    figCal = figure('Position', figPos, 'Name', 'fig18_CalibCheck');
    plot(theta*180/pi, 20*log10(B1a/max(B1a)), 'b-', 'LineWidth', 1); hold on;
    plot(theta*180/pi, 20*log10(B1b/max(B1b)), 'r-', 'LineWidth', 1);
    plot(theta*180/pi, 20*log10(BIdA/max(BIdA)), 'k--', 'LineWidth', 1);
    plot(theta*180/pi, 20*log10(BIdB/max(BIdB)), 'k:', 'LineWidth', 1);
    ylim([-40 2]); grid on;
    xlabel('角度 (deg)'); ylabel('归一化幅度 (dB)');
    legend({sprintf('全阵 64 元(角反, 幅度归一化+校准)'), ...
            sprintf('最长行 %d 元(角反, 同上)', nInRow(rMax)), ...
            '理想 全阵(等幅导向 0°)', '理想 最长行(等幅导向 0°)'}, ...
           'Location', 'northeast');
    title(sprintf('校准源角反 × 校准表 角度谱(R = %.2f m, 近场已补偿)', R1));

    calibInfo = struct('m', m, 'rBin', rBin, 'dopBins', dopFound, ...
                       'phRaw', phRaw, 'phCode', phCode, 'phErr', phErr, ...
                       'phNF', phNF, 'RfarField', RfarField, 'figCal', figCal);
    fprintf('0°校准       : 目标1 @ bin%d, 各TX速度bin [', rBin);
    fprintf('%d ', dopFound); fprintf('], 近场补偿后表相位跨度 rms %.1f°, |c|=1, 免近场补偿距离 %.2f m\n', ...
        rad2deg(sqrt(mean(phErr(:).^2))), RfarField);
end

% ----------------------------------------------------------------
% 目标通道提取: 已知 R/v → 距离bin + 各TX速度bin(含步进频耦合) → 邻域搜索
% ----------------------------------------------------------------
function [m, rBin, dopFound, phCode] = extractCh(fftCfg, slots, nGroups, nChirps, win, ...
                                                  slope, dfBin, fStep, lambda, Tc, R, v)
    cc = physconst('lightspeed');
    nTX = numel(slots);
    nRX = size(fftCfg.data_sva, 2);
    tau = 2*R/cc;
    fb = slope*tau + 2*v/lambda;               % 差拍(含多普勒项)
    rBin = round(fb / dfBin);
    % 速度 bin = 多普勒(2v/λ×Tc×N) + DDMA 码位×N/16 + 步进频耦合 f_step·τ·N
    dopExp = round(mod( (2*v/lambda)*Tc*nChirps + slots*nChirps/nGroups ...
                  + fStep*tau*nChirps, nChirps ));
    m = zeros(nTX, nRX);
    dopFound = zeros(1, nTX);
    phCode = zeros(1, nTX);
    for i = 1:nTX
        best = 0;  bR = rBin;  bD = dopExp(i);  mi = zeros(nRX, 1);
        for rr = max(0, rBin-1) : min(size(fftCfg.data_sva,1)-1, rBin+1)
            for dd = max(0, dopExp(i)-2) : min(nChirps-1, dopExp(i)+2)
                h = fftCfg.data_sva(rr+1, :, dd+1).';      % [nRX×1]
                if norm(h) > best, best = norm(h);  bR = rr;  bD = dd;  mi = h; end
            end
        end
        m(i, :) = mi.';  dopFound(i) = bD;
        codeRef = exp(-1i * deg2rad(slots(i)*360/nGroups) * (0:nChirps-1));
        phCode(i) = angle( sum(codeRef(:) .* win .* exp(-2i*pi*bD*(0:nChirps-1).'/nChirps)) );
    end
end

% ----------------------------------------------------------------
% 近场相位 [nTX×nRX]: φ = -k(Rt_i + Rr_j - 2R), 目标在正视 y 轴
% ----------------------------------------------------------------
function phNF = nearFieldPh(antenna_t, antenna_r, R, k)
    Rt = sqrt(antenna_t(1,:).^2 + R^2 + antenna_t(3,:).^2);
    Rr = sqrt(antenna_r(1,:).^2 + R^2 + antenna_r(3,:).^2);
    phNF = -k * (Rt(:) + Rr(:).' - 2*R);
end
