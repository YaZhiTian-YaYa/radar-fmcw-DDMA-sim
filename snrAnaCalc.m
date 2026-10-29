function an = snrAnaCalc(rfInfo, prm, Rm, rcsDb, nCoh, nRXavg, isRealSamp, winStat, filtG, cohStat, cohGainDb)
%SNRANACALC 谱线 SNR 解析计算(工程近似, 统一物理口径——95/96/97/97_IQ 各工程同一子函数)
%   SNR = Pr − 10lg(kT·NF) + 10lg(Ts) + 10lg(nCoh) + winStat − 3.0×(实采)
%         + 滤形 − PLL跟踪线损失 − 自扩散抬底
%   各项: 雷达方程接收功率(TXRX 非相参积累口径, Pt = 每 TX 功率——8TX 对齐 + 8RX 功率
%   平均均不改均值 SNR, 解析与单线同式) / 噪声 PSD(仿真同温同 NF) / 快时间相干积分(=
%   距离 bin 宽倒数, 采样窗 Ts) / 慢时间 nCoh chirp 相干(DDMA 码移不损失) / 窗与统计
%   净项(随度量口径标定, 见下) / 实采样镜像劈半罚 / 滤形(缺省 = LPF×HPF×数字半带×
%   带边均衡 @目标差拍) / PLL 两项(见下)。
%   【PLL 跟踪两项(2026-10-29 显式化, 勿再藏进 winStat)】逐 chirp iid 随机相位 φ_n
%   (扫频跟踪 σδf = 系数×B/Ts 与跳频误差 σf, 残余相位 = 2π·σ·τ, τ=2R/c)把谱线功率
%   按 e^{−σφ²} 劈分, 两项后果都随 R² 增长:
%     ① 线相干损失 = 10lg(e^{−σφ²}) = −4.343·σφ² (300m/系数7.6e-9 → −2.8 dB);
%     ② 自扩散抬底: 劈出的 (1−e^{−σφ²}) 份额被 DDM 8 根码线平铺到整个多普勒行
%        (逐 chirp iid → 平铺), 抬高行中位噪底 = 10lg(1 + 8·(1−e^{−σφ²})·S0/nCoh),
%        S0 = 无 PLL 项线 SNR(线性) —— 弱目标(S0·8/nCoh≪1)时→0, 强目标显著
%        (20dBsm@210m 仿真观测抬 +5.4 dB), 即强目标动态范围天花板: SNR → nCoh/(8(1−e^{−σφ²}))。
%   两项合成仿真测量 SNR = S0·e^{−σφ²}/(1 + 8(1−e^{−σφ²})·S0/nCoh), 2026-10-29 对
%   210m/20dBsm 与 300m/10dBsm 两目标验证一致(±0.5 dB)。
%   【winStat 构成(非物理增益, 勿按积累解读)】双维窗峰损/噪底带宽收窄之净 + 噪底
%   中位统计 + 跨 bin 跨瓣 + PN/非线性肩部等小量; 逐天线损耗差的均值功率已显式列为
%   chTerm(2026-10-29 从 winStat 拆出, 勿再藏回)。【随度量口径与窗类型标定】:
%   97 三芯片(1T2R/局部行中位) = −6.0; 96(4RX/全帧/局部行中位) = −3.0; 95 DDMA(点云
%   nonCohSNR 口径) 按窗分档: Hanning = −2.6 / SVA = −1.1(2026-10-29 同探针标定,
%   SVA 收回矩形 ENBW 的大半, 见主流程 winStatTbl)。两缺省分别取 −6.0 与
%   rfInfo.filtGainAtR(Rm)。
%   【通道相参积累 SNR(2026-10-29 新增, 同日改 DBF 口径)】nTX×nRX 通道经幅相校准后同相
%   叠加: 公式口径 = 单线 SNR + 10lg(nTX·nRX)(64 通道无幅度差/俯仰 0° 时 +18.1 dB)。
%   【DBF 口径(优先, 用户 2026-10-29, 同日改归一化)】传入 cohGainDb 时用实际阵列
%   计算的增益(cohGainDbf.m): 10lg(|Σ_c e^{jk·z_c·sinθel}|²/N), 通道幅度归一化
%   (算法方案口径: DBF 权重归一化, 不含 TX_A/RX_A 锥削 —— 幅度差均值已由 chTerm
%   显式列项, 勿重复计入), z_c = 阵元俯仰坐标 —— 俯仰角≠0 时一维方位 FFT 对不齐
%   的 z 行相位自动降低增益(此时"强行 N 元相参"低于公式值, 属预期),
%   θel=0 时恒等于公式值 10lg(N)。
%   但仿真口径在【角度维取峰】(cohPeakPow = 方位谱取 max): 噪声参考(全图中位)含
%   max-of-角度bin 统计抬升 —— 理论值由 cohStatCalc.m 蒙特卡洛计算(随雷达结构
%   nChan/nFft/nRXavg 自动变, 勿写死常数), cohStat 传【理论统计值, 勿凑差为零】;
%   残差(理论略优于仿真)如实显示。
%   【rxTerm = 0】RX 功率平均不增加(平均)SNR —— 收益是噪声更平均、峰值更低 → 门限
%   可降; "非相参增益"是门限收益的近似等效(饱和型), 非真实 SNR 提升(用户结论
%   2026-10-28; 曾误写 +10log(nRX) 高估 10lg(nRX) dB, 已修)。
    if nargin < 8 || isempty(winStat), winStat = -6.0; end
    if nargin < 9, filtG = []; end
    if nargin < 10 || isempty(cohStat), cohStat = 0; end
    if nargin < 11, cohGainDb = []; end
    cL    = physconst('lightspeed');
    tauR  = 2*Rm/cL;                                             % 双程时延 (s)
    sig2  = 0;                                                   % 逐 chirp iid 相位方差累积
    if isfield(rfInfo, 'pllTrackOn') && rfInfo.pllTrackOn
        sig2 = sig2 + (2*pi*rfInfo.pnTrackCoef*(prm.B/prm.Ts)*tauR).^2;
    end
    if isfield(rfInfo, 'hopErrOn') && rfInfo.hopErrOn
        sig2 = sig2 + (2*pi*rfInfo.hopErrHz*tauR).^2;
    end
    rcs   = 10^(rcsDb/10);
    Pr_W  = 10^(rfInfo.Pt_dBm/10)*1e-3 * 10^(rfInfo.Gt_dBi/10) * 10^(rfInfo.Gr_dBi/10) ...
            * rfInfo.lambda^2 * rcs / ((4*pi)^3 * Rm^4);        % 远场单站雷达方程(单 TX)
    an.Pr_dBm    = 10*log10(Pr_W*1e3);                           % 接收功率 (dBm)
    an.psd_dBmHz = 10*log10(rfInfo.Pn0_W_Hz*1e3);                % 噪声 PSD (dBm/Hz) = 10lg(kT·NF)
    an.tTerm     = 10*log10(prm.Ts);                             % 快时间相干积分(= bin 宽倒数)
    an.dTerm     = 10*log10(nCoh);                               % 慢时间相干(nCoh chirp)
    an.rxTerm    = 0;                                            % RX 功率平均不计增益(见上)
    an.winStat   = winStat;                                      % 窗与统计净项(口径标定值)
    an.sampP     = 3.0*double(isRealSamp);                       % 实采镜像劈半罚(复采 0)
    if isempty(filtG)
        an.filtG = rfInfo.filtGainAtR(Rm);                       % LPF×HPF×数字半带×均衡 @目标 fb
    else
        an.filtG = filtG;                                        % 外部口径(如 ADC 片噪份额)
    end
    nTxL  = 1;  if isfield(rfInfo, 'nTX'), nTxL = rfInfo.nTX; end    % 平铺行内的码线数
    an.chTerm    = 0;                                            % 通道锥削均值损耗(dB, 2026-10-29
    if isfield(rfInfo, 'chLossDb')                               %   显式化, 原藏在 winStat):
        an.chTerm = rfInfo.chLossDb;                             %   逐天线损耗差(TX_A/RX_A)的
    end                                                          %   64 通道平均功率 ≈ −0.8 dB
    an.pllTerm   = -4.343*sig2;                                  % ① PLL 线相干损失(dB)
    s0Lin        = 10^((an.Pr_dBm - an.psd_dBmHz + an.tTerm + an.dTerm + an.rxTerm ...
                       + an.chTerm + an.winStat - an.sampP + an.filtG)/10);  % 无 PLL 项线 SNR
    an.selfFloor = 10*log10(1 + nTxL*(1 - exp(-sig2))*s0Lin/nCoh);   % ② 自扩散抬底(dB)
    an.SNR_dB    = 10*log10(s0Lin*exp(-sig2)) - an.selfFloor;
    if isempty(cohGainDb)
        an.cohTerm = 10*log10(nTxL * nRXavg);                    % 公式口径(无阵列信息时)
    else
        an.cohTerm = cohGainDb;                                  % DBF 口径(实际阵列, 含锥削/俯仰失配)
    end
    an.cohStat   = cohStat;                                      % 相参口径统计项(cohStatCalc.m 理论值)
    an.SNRcoh_dB = an.SNR_dB + an.cohTerm + an.cohStat;          % 通道相参 SNR(对照点云 cohSNR)
end
