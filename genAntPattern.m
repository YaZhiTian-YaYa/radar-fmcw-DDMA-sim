function antPat = genAntPattern(nTX, nRX)
%GENANTPATTERN 车载雷达天线方向图(逐天线复数二维矩阵): 方位/俯仰各 −90°~+90°, 1° 网
%   antPat = genAntPattern(nTX, nRX)
%   【仅定义数据 + 检查用, 未接入任何模块】—— 方向图经图1~6 人工检查确认后,
%   再决定接入哪些模块; 在此之前其他模块一律不得使用 antPat。
%   主形状工程锚点(车载雷达典型):
%     峰值     : (az, el) = (0°, 0°) 处 ≈15 dBi
%     方位面   : 不组阵(单元宽波束): ±70° 仍约 0 dBi, 无副瓣
%     俯仰面   : 组阵: ±20° 降到 0 dBi, 主瓣裙边外每侧 5 个等高副瓣(−25 dB, 近切比雪夫)
%   可分离近似: G(az, el) = G0 + G_az(az) + G_el(el), 锚点全部参数化(旋钮见下)。
%
%   通道间不一致性(模拟实测方向图的逐天线差异, 当前用随机模型代替):
%     幅度: 逐天线增益偏置 ±0.4 dB + 全角度小起伏(±0.25/±0.15 dB)
%           + 大角度波纹(±1.5 dB, |az|>35° 渐起 —— 大角度性能退化逐天线不同)
%     相位: 恒定偏置 ±10° + 方位/俯仰线性斜率(指向微偏) + 方位/俯仰正弦起伏 ±5°
%   rng 种子固定(7), 每次运行生成的 16 张"实测替代图"完全相同, 结果可复现;
%   将来拿到实测方向图后, 逐个替换 txGainC/rxGainC 对应页即可(维度不变)。
%   复数约定: 幅值 = 10^(G_dB/20)(场量纲, |·|² = 功率增益), 相位 = 远场相位(rad)。
%
%   输出 antPat:
%     .txGainC [nEl × nAz × nTX] / .rxGainC [nEl × nAz × nRX]  复数方向图
%     .azDeg [1×nAz] / .elDeg [nEl×1]  角度网格
%     .gainPeakDb / 各锚点旋钮(带出便于检查对照)

    % ---------- 主形状锚点旋钮 ----------
    gainPeakDb = 15;        % 峰值增益 (dBi) @ 0°/0°
    azEdgeDeg  = 70;        % 方位锚点: |az| = 此值处落到 0 dBi(高斯滚降, 无副瓣)
    azFallDb   = 15;        % 方位锚点落差 (dB)
    elEdgeDeg  = 20;        % 俯仰锚点: |el| = 此值处主瓣降到 0 dBi
    elFallDb   = 15;        % 俯仰主瓣落差 (dB)
    slLevelDb  = 25;        % 俯仰副瓣电平 (dB, 相对峰值; 切比雪夫式等电平)
    slNlobe    = 5;         % 每侧副瓣个数
    slNullDepthDb = 20;     % 副瓣峰到峰间零深 (dB)
    slP        = 0.5;       % 瓣形指数(|sin|^p)

    % ---------- 通道间不一致性旋钮(随机模型幅度) ----------
    ampOffMaxDb = 0.4;      % 逐天线增益偏置上限 (dB)
    ripAzMaxDb  = 0.25;     % 方位向形状起伏上限 (dB, 全角度)
    ripElMaxDb  = 0.15;     % 俯仰向形状起伏上限 (dB)
    ripFarMaxDb = 1.5;      % 方位大角度波纹上限 (dB): |az|>35° 渐起, 边缘满幅 ——
                            %   模拟大角度性能退化, 且逐天线退化规律不同(实测常见)
    phOffMaxDeg = 10;       % 相位恒定偏置上限 (deg)
    phSlpMaxDeg = 0.15;     % 相位线性斜率上限 (deg/deg, 等效指向微偏)
    phRipMaxDeg = 5;        % 相位起伏上限 (deg, 方位/俯仰各一组)

    % ---------- 网格与主形状 ----------
    azDeg = (-90 : 1 : 90);                     % 方位网格 (deg)
    elDeg = (-90 : 1 : 90).';                   % 俯仰网格 (deg)
    a = abs(azDeg);  e = abs(elDeg);
    G_az = -azFallDb * (a / azEdgeDeg).^2;      % 方位: 高斯滚降到 ±70° 处 0 dBi
    mainDb = -elFallDb * (e / elEdgeDeg).^2;    % 俯仰主瓣裙边
    eX = elEdgeDeg * sqrt(slLevelDb / elFallDb);            % 裙边降到副瓣电平的角度 ≈25.8°
    slLobeDeg = (90 - eX) / (slNlobe + 0.5);                % 瓣距: 恰容 5 峰, 末峰后压零
    inSL = e > eX;
    lobeDb = -slLevelDb - slNullDepthDb * abs(sin(pi * (e - eX) / slLobeDeg)).^slP;
    G_el = mainDb;
    G_el(inSL) = max(mainDb(inSL), lobeDb(inSL));
    gainDb0 = gainPeakDb + G_el + G_az;         % 主形状 (dBi) [nEl × nAz]

    % ---------- 逐天线差异(固定种子, 可复现) ----------
    rng(7);
    nAnt = nTX + nRX;
    ampOff = ampOffMaxDb  * (2*rand(1,nAnt) - 1);      % 增益偏置
    ripAzA = ripAzMaxDb  * (2*rand(1,nAnt) - 1);       % 方位起伏幅度(全角度, 小)
    ripAzF = 2 + 2*rand(1,nAnt);
    ripAzP = 2*pi*rand(1,nAnt);
    ripElA = ripElMaxDb  * (2*rand(1,nAnt) - 1);
    ripElF = 2 + 2*rand(1,nAnt);
    ripElP = 2*pi*rand(1,nAnt);
    ripFarA = ripFarMaxDb * (2*rand(1,nAnt) - 1);      % 大角度波纹幅度(逐天线不同)
    ripFarF = 3 + 4*rand(1,nAnt);                      % 波纹频率(周期/180°)
    ripFarP = 2*pi*rand(1,nAnt);
    envFar = max(0, (a - 35)/55).^1.5;                 % 波纹包络: 35° 外渐起, 90° 满幅
    phOff  = phOffMaxDeg * (2*rand(1,nAnt) - 1);       % 相位偏置
    phSlpA = phSlpMaxDeg * (2*rand(1,nAnt) - 1);       % 方位相位斜率(指向微偏)
    phSlpE = phSlpMaxDeg * (2*rand(1,nAnt) - 1);
    phRipA = phRipMaxDeg * (2*rand(1,nAnt) - 1);       % 方位相位起伏
    phRipF = 2 + 2*rand(1,nAnt);
    phRipP = 2*pi*rand(1,nAnt);
    phRipEA = phRipMaxDeg * (2*rand(1,nAnt) - 1);      % 俯仰相位起伏(独立一组)
    phRipEF = 2 + 2*rand(1,nAnt);
    phRipEP = 2*pi*rand(1,nAnt);

    txGainC = zeros(size(gainDb0,1), size(gainDb0,2), nTX);
    rxGainC = zeros(size(gainDb0,1), size(gainDb0,2), nRX);
    for k = 1:nAnt
        ripDb = ripAzA(k)*sin(2*pi*ripAzF(k)*azDeg/180 + ripAzP(k));       % [1×nAz] 全角度小起伏
        ripFar = ripFarA(k)*sin(2*pi*ripFarF(k)*azDeg/180 + ripFarP(k)) .* envFar;  % 大角度波纹
        ripDbE = ripElA(k)*sin(2*pi*ripElF(k)*elDeg/180 + ripElP(k));      % [nEl×1]
        gDb = gainDb0 + ampOff(k) + ripDb + ripFar + ripDbE;               % 该天线幅度 (dBi)
        phDeg = phOff(k) + phSlpA(k)*azDeg + phSlpE(k)*elDeg ...
              + phRipA(k)*sin(2*pi*phRipF(k)*azDeg/180 + phRipP(k)) ...
              + phRipEA(k)*sin(2*pi*phRipEF(k)*elDeg/180 + phRipEP(k));    % 该天线远场相位
        gC = 10.^(gDb/20) .* exp(1i*deg2rad(phDeg));
        if k <= nTX, txGainC(:, :, k) = gC; else, rxGainC(:, :, k - nTX) = gC; end
    end

    antPat = struct('azDeg', azDeg, 'elDeg', elDeg, ...
        'txGainC', txGainC, 'rxGainC', rxGainC, ...
        'gainPeakDb', gainPeakDb, 'azEdgeDeg', azEdgeDeg, 'azFallDb', azFallDb, ...
        'elEdgeDeg', elEdgeDeg, 'elFallDb', elFallDb, ...
        'slLevelDb', slLevelDb, 'slNlobe', slNlobe, 'slNullDepthDb', slNullDepthDb);
end
