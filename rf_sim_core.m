function [adcRawData, rfInfo] = rf_sim_core(nSamples, nChirps, TX_A, RX_A, antenna_t, antenna_r, ddmShift, prm, targets, swOv, printCfg)
%RF_SIM_CORE 8T8R MIMO 毫米波汽车雷达 FMCW 射频链路仿真子函数(DDMA, 实采样 ADC)
%   adcRawData = rf_sim_core(nSamples, nChirps, TX_A, RX_A, antenna_t, antenna_r)
%   输入(仅 6 项; 全部射频开关/旋钮/目标设置都在本函数内部, 不经主函数):
%     nSamples    每 chirp 实采样点数
%     nChirps     每帧 chirp 数(多普勒 FFT 长度)
%     TX_A, RX_A  [1×8] TX/RX 通道幅度加权
%     antenna_t, antenna_r  [3×8] TX/RX 天线相位中心坐标 [x;y;z] (m)
%     ddmShift    [1×nTX]  各 TX 的 DDMA 码槽位(0 起; 主函数由 encodeSeq 算好传入)
%     swOv        可选: 结构体覆盖上方开关/旋钮(单因素对比实验用)
%     printCfg    可选: 是否打印射频配置块(默认 true)。同一套配置跑两遍时
%                 (校准遍+目标遍)后跑的一遍传 false, 配置只打印一次;
%                 目标表/ADC量程/溢出计数等数据相关行仍每遍都打。
%   输出:
%     adcRawData  [nSamples×nChirps×nRX] ADC 原始数据, int16 counts 口径
%     (含全部射频损伤; 16bit 溢出环绕等数字异常也在本函数内闭合)
%% ========================================================================
%  rf_sim_core.m —— RF 链路仿真子函数(功能最全版, 含全部射频损伤模型)
%  8T8R MIMO 毫米波汽车雷达 FMCW 仿真（DDMA 波形，实采样 ADC 版本）
%
%  -------------------------- 实采样仿真原理 ------------------------------
%  复采样(IQ)链路：混频后 I/Q 两路 ADC，得到复数差拍信号，可观测带宽 ±fs/2；
%  实采样链路：  混频 + 低通滤波后只有一路 ADC，信号为实数，可观测带宽 [0, fs/2)。
%
%  仿真做法：射频/基带等效信号仍按复数生成（物理正确），在 ADC 环节取 real()。
%            real{ s_tx * conj(s_rx) } 恰好就是实混频器输出的差频项
%            cos(φ_tx − φ_rx)，因此对复差拍信号取实部即为实采样的正确模型。
%
%  与复采样相比的处理差异（实采样核心）：
%   1) 实信号 FFT 频谱共轭对称：每个目标在 +fb 和 −fb（即第 N−k 个 bin）
%      各有一个镜像峰。距离维 FFT 只保留前 nSamples/2+1 个 bin（0 ~ fs/2），
%      镜像峰直接丢弃 —— 若把全部 bin 都画出来，会看到"目标翻倍"，这就是
%      直接改参数+取 real 后图形不对的根本原因；
%   2) 同幅度目标的峰值比复采样低 6 dB（一半能量分给了镜像峰）；
%   3) 最大不模糊差拍频率由 fs/2 决定 → 最大不模糊距离 Rmax = c*fs/(4*slope)；
%   4) 速度信息保存在距离 bin 复数相位随 chirp 的演化中，实采样不损失
%      速度维信息，距离 FFT 之后的 2D-FFT 处理流程与复采样完全一致；
%
%  链路仿真项（编号与顺序同下方「仿真控制开关区」）：
%   1) 步进频：chirp 间起始频移（占比 = 跳频带宽/扫频带宽，端点锚定 →
%      帧射频范围 [fStart, fStart+(1+cr)·B]）；目标与泄漏码族随距离产生多普勒耦合偏移；
%      跳频误差（enableHopErr, 本文件新增）：每次跳频的起始频率相对理想值
%      N(0, hopErr_Hz²) 正态随机（σ=hopErr_Hz）、逐 chirp 独立【不累积】→ 差拍频率不变, 残余为
%      逐 chirp 随机相位 2π·δf·τ ∝ 目标时延 → 远目标多普勒码线去相关、
%      展宽成噪底（近目标几乎不受影响）；
%   2) 接收机白噪声：空暗室底噪 Pn = k·T0·B_NEB·NF（与目标无关），
%      经链路增益 G_rx 折算到 ADC 域与信号同一绝对刻度后叠加；
%   3) IF 两级模拟高通（enableHpf 可关）：压 0 频泄漏载波与近距离目标幅度
%      （防 ADC 饱和），巴特沃斯、截止/阶数可配置；近距底噪被压低后随频率
%      按滤波器形状回升（对数频率下呈抛物线）；hpfFloor_dB = 绕过高通的
%      泄漏旁路耦合（近距码梳高度的直接控制量之一）；
%   4) IF 抗混叠低通：射频链路建模幅度按雷达方程加权（∝ sqrt(σ)/(Rt·Rr)），
%      巴特沃斯模拟低通（平滑滚降，非数字硬截断）；超量程目标先被 LPF 压制，
%      残余分量采样后混叠成"鬼影"（视在距离 ≠ 真实距离、多普勒符号翻转）；
%      滤波在抽取之前 → 与硬件"先模拟滤波、后采样折叠"次序一致；
%   5) TX 移相器实际实现：DDMA 相位量化到 N-bit 移相器可达值 + 逐 chirp
%      正态随机误差（两个独立开关）；码步进非移相器 LSB 整数倍时，不可实现
%      相位按就近值替代 → 空子带出现泄漏（部分子带，位置取决于替代值）；
%   6) 带内起伏：天线+射频链随回波扫频射频的确定性响应，在 chirp 内表现
%      为调幅（AM）→ 距离维 FFT 旁瓣展宽（约 -20dB 量级）；峰峰值可调，
%      8TX/8RX 可独立设置形状相位；
%   7) 扫频非线性：共锁相环的确定性非线性 ε(t)，四种形状可选（单凸/单凹/
%      先凸再凹/先凹再凸，端点锚定 → 带宽不变；峰值 = 非线性度×B），回波与
%      本振相减后残余调频 ≈ τ·dε/dt ∝ 距离 → R-FFT 展宽随距离增大；形状
%      不同的残余调频轨迹（单向扫/往返扫）会产生不同的 R-FFT 谱形；
%   8) 相位噪声 + PLL 环路滤波器（两项独立开关）：共锁相环 PN（dBc/Hz 拟合
%      点可编辑）随机相位法合成时域实现；dechirp 共模对消后残余 PSD ∝
%      |1−e^{−j2πfτ}|² → 目标谱两侧出现肩部，肩部电平随距离（τ）抬高；
%      PLL 环路低通整形使输出 PN 在环路带宽处陡峭截断（所有目标肩部截断
%      偏移相同，与距离无关）；
%   9) 零频泄漏（enableLeakDC 独立开关）：片内直通 τ=0 → 差拍=0 → bin0 出
%      码梳；bin0 行慢时间序列为实数 → 实采样多普勒共轭镜像（默认关）；
%  10) 近距泄漏多径族（enableLeakNear 统一开关）：leakPaths 每行一条非零延迟
%      路径 [差拍频率Hz, 电平]（地弹+屏蔽腔多次反射 / 环境残余反射 / 天线互耦），
%      差拍 = 斜率×路径延迟 → 落在各自距离 bin → 近距码梳（复数行→无实采样镜像）；
%      leakRes_dB = 全部泄漏路径的全局电平修正（近距码梳高度直接控制量）；
%  11) TX-RX 收发隔离泄漏【总开关】（enableLeakage）：泄漏载波携带 DDMA 与
%      每路径固定相位；隔离度矩阵 isoMat 决定各路径幅度（零频/近距多径块
%      的细分开关需与总开关同时开）；泄漏 PN 裙边（leakPnTau）+ 每通道
%      1/f^α 闪烁噪声（PSD ∝ 泄漏功率）→ 近距底噪抬高且随距离缓降；
%  12) TX-TX 发发隔离泄漏 + PA 三阶交调（两项独立开关，交调依赖耦合存在）：
%      TXm 信号串入 TXn 链路经 PA 放大后从 TXn 天线再辐射（携带 TXm 的
%      DDMA 码、固定随机相位、几何路径按 TXn）→ DDMA 各槽相干污染（虚拟
%      阵列误差，幅度 ~耦合度）；多码复合信号过 |x|²x 非线性 → 2k_n−k_m
%      组合码落入空子带 → 空子带相同距离 bin 出现微弱目标信号；
%  13) RX-RX 收收隔离泄漏（enableRxCoupling 可关）：RXm 信号串入 RXn 前端
%      （固定随机相位、线性叠加于本链噪声/滤波之前）→ 通道间串扰（虚拟
%      阵列误差）。RX 交调不模拟：µV 级小信号的 LNA 交调 ≪ 热噪不可见；
%  14) ADC 时域量程（int16 counts 口径）：ADC 原始数据以 ±32768 counts 为单位
%      输出（不模拟整数量化）；时域样点不得超出 ±32768（超限=削顶, 告警提示），
%      FFT 后超过满量程属积累增益、正常。固定增益 G_rx_dB 标定：无强目标时
%      （弱目标+泄漏+底噪）时域约 ±200 counts，强目标最多用到接近满量程；
%  15) 馈线长度色散（enableFeedDisp, 本文件新增）：8TX+8RX 共 16 条馈线不等长，
%      各通道引入与频率成正比的相位差 φ(f)=2π(f−f_cal)·ΔL/c；假设已在全带
%      （扫频 B+步进频）中心频率 f_cal 校准 → 中心频率处 16 通道相位相等，
%      偏离中心按线长差计算 → 阵列通道相位失配随频率变化（测角误差/虚拟
%      阵列畸变，随偏离中心频率增大）；
%  ========================================================================

%% ========================== 仿真控制开关区 ==========================
%  组织约定：每个功能一块 = 一行 [开关] enableXxx + 该功能的全部 [旋钮] 参数，
%  块间用横线隔开；新增功能的开关与参数也按此放在对应块内。
%  所有开关均为逻辑量 true/false；[旋钮] 只出现在所属开关的下方。

% ---------- 1. 步进频（chirp 间起始频移） ----------
coupling_ratio = prm.cr;   % 步进频占比 = 总跳频带宽 / 扫频带宽 B（0 = 无步进频, 由主函数 prm 注入）
enableHopErr  = false;     % [开关] 步进频跳频误差【本文件新增功能】
hopErr_Hz     = 1e5;       % [旋钮] 跳频误差标准差 σ (Hz, 正态分布)：第 n 个 chirp 的
                           %   起始频率 = 理想值 + N(0, hopErr_Hz²)。逐 chirp 独立抽取、
                           %   【不累积】——始终以理想频率为基准, 与前面 chirp 的
                           %   误差无关（不是在上一次实际值上继续累加）。
                           %   物理含义：频率综合器每次跳频不准, 整个 chirp 随起始
                           %   频率整体偏移；本振与发射共源 → 差拍频率不变
                           %  （= slope·τ）, 残余为逐 chirp 随机相位 2π·δf·τ,
                           %   随目标时延/距离增大 → 远目标多普勒码线展宽成噪底

% ---------- 2. 接收机白噪声（空暗室底噪, 与目标无关） ----------
enableNoise = true;        % [开关] 接收机白噪声
NF_dB       = 11;          % [旋钮] 噪声系数 (dB)，控制底噪高低（典型 10~20）

% ---------- 3. IF 模拟高通 ×2（片上连续时间有源滤波器: Gm-C / active-RC） ----------
% 真实 MMIC/SoC 的 IF 高通是两节【一阶】串联(TI mmWaveStudio 的 HPF1/HPF2 各自
% 选档正是这个结构, 典型档位 175/235/350/500 kHz 一类)。只挡直流偏置、收发直耦
% 和极近距杂波 —— 阶数高了群延迟大、冲激响应振铃, 会把近距小目标弄坏, 所以
% 故意做温和: 每节一阶(无振铃), 过渡带缓, 两节串联等效二阶但阻尼最大。
enableHpf  = true;         % [开关] 两级模拟高通
gpuMode    = 'auto';       % [旋钮] 目标回波张量算力: 'auto'(有GPU用GPU) | 'on' | 'off'(纯CPU张量)
hpf1Fc     = 175e3;        % [旋钮] 高通 1(混频输出级)截止 (Hz) —— TI 典型档 175/235/350/500k
hpf1Order  = 1;            % [旋钮] 高通 1 阶数(固定 1: 一阶节)
hpf2Fc     = 500e3;        % [旋钮] 高通 2(ADC 前级)截止 (Hz) —— 同上档位表
hpf2Order  = 1;            % [旋钮] 高通 2 阶数(固定 1: 一阶节; 两节串联=等效二阶无振铃)
hpfFloor_dB = -60;         % [旋钮] 泄漏旁路耦合 (dB)：部分泄漏绕过高通直达后级的
                           %   等效耦合度 → 决定近距码梳高度（±5 dB → 码梳 ±5 dB）。
                           %   实现方式：滤波输出后叠加 泄漏×10^(hpfFloor_dB/20)
%（不用"响应下限 max()"实现：那会把冲激响应拉长到数千
%   样点，破坏分段线性卷积 → 段尾环绕伪影）

% ---------- 4. IF 抗混叠低通（片上连续时间有源滤波器: Gm-C / active-RC, 巴特沃斯类） ----------
% 【设计约束链, 改参数前必读】nSamples=1024/nChirps=256 是硬件内存上限 →
%   fb(350m)/奈奎斯特 = 350/(Ns·c/4B) = 97.5%(与 fs 无关, 不可绕过)——通带边与
%   混叠边只隔 2.5%, 单级模拟滤波器"带内缓降 + 带边强抑制"不可兼得, 三选一:
%   A 牺牲距离分辨率: B 降到 ~120 MHz(ΔR 1.25 m) → 350m 只占奈奎斯特 55%,
%     4 阶即可"缓降+抑制"(350m −0.3 dB / 奈奎斯特 −12 dB / 45M −25 dB);
%   B 保 ΔR=0.70 m(B=214 MHz 不动): fc 上移让带内缓降, 代价是带边抑制弱
%     (370~400 m 强目标鬼影折回 300~350 m 段, 多普勒符号翻转可甄别);
%   C 对齐实测: 若目的就是复现真实芯片滤波器, 不设计——拿实测噪底曲线的
%     droop 形状反推 lpfOrder/lpfFc(与 ADC 噪声标定同思路)。
% 【两段式已解耦该矛盾】(见下方两段式块): 模拟级宽通带只管镜像区、带内零形状,
%   数字半带悬崖定义信号带边 —— 缓降与强抑制兼得。当前响应(2026-10-28 复算, 含
%   带边噪声均衡 enableHbComp): 信号 @16.7M(300m) −0.0 | @19.5M(350m) −0.2(残余
%   droop) | @20M(奈奎斯特) −3.4 | @24M −95 | @30M −104 dB(悬崖不变); 噪底(折叠
%   对合计×均衡)全带 ±0.2dB 平(对齐芯片实测)。输出 NEB ≈ 20.2M(均衡后数值积分)。
enableLpf = true;          % [开关] 抗混叠低通(模拟域)
lpfOrder  = 3;             % [旋钮] 阶数(巴特沃斯, SOS 实现; 模拟域只管宽处, 可缓)
% ---- 【两段式链路 2026-10 定稿, 对齐真实芯片"宽模拟 + 数字悬崖"架构】----
% chipOsr = 降采样倍数(芯片 ADC 实采样率 = chipOsr×fs):
%   1: 模拟LPF → ADC@fs → 直出(旧行为, 带边由模拟滤波器定义)
%   2: 模拟LPF(宽) → ADC@2fs → 数字半带LPF → ÷2
%   4: 模拟LPF(宽) → ADC@4fs → 半带LPF1 → ÷2 → 半带LPF2 → ÷2(真实芯片 160M 链)
% 【照 96 修法上提主函数】降采样倍数主函数与射频函数都要用 —— 主函数用它算真正
% 的数据率与雷达各项指标(fs = fsAdc/chipOsr), 射频核用它设计数字半带链并抽取 →
% 唯一定义处在主函数、经 prm.chipOsr 注入此处(本文件不再自带旋钮)。
% 模拟截止 = 0.45×ADC率(远离信号带, 带内零形状, 只管 >ADC奈奎斯特 的镜像区);
% 信号带边(fs/2)由末级半带 FIR 的陡截止定义(数字墙, 线性相位)。
% 仿真内部表示率 = 2×ADC率(编程手段: 用离散信号模仿模拟滤波作用于真实频率,
% 与硬件倍率无关; chipOsr=1 时内部 2× 即旧行为 osf=2)。
chipOsr   = prm.chipOsr;   % 降采样倍数(1 | 2 | 4; 唯一定义处在主函数, 经 prm 注入)
assert(any(chipOsr == [1 2 4]), 'chipOsr 只支持 1/2/4');
hbTaps    = [63 127];      % [旋钮] 各级半带 FIR 抽头数(奇数, 4k±1; 第 s 级取 hbTaps(min(s,end)))
                           %   末级是带边悬崖, 要长: 127 抽头@80M 域 → 过渡带 ≈±1.8 MHz
hbBeta    = 8;             % [旋钮] 半带 Kaiser 窗 β(8 ≈ −60 dB 阻带)
enableHbComp = true;       % [开关] 带边噪声均衡(照 96 修法 enableCicComp; 2026-10-28
                           %   对齐芯片实测: 实测远端噪底【无】滤波形状 —— 不均衡时末级半带
                           %   过渡带把白噪压出 ~2dB 下垂裙边, 实测平)
compFreezeV  = 1.0;        % [旋钮] 均衡冻结点(v = f/(fs/2), 补偿到此为止; 1.0 = 全带)
hbCompTaps   = 127;        % [旋钮] 均衡 FIR 抽头数(奇数; 63 → ±0.4dB 纹波, 127 → ±0.2dB)
fsAdc     = prm.fsAdc;     % 芯片 ADC 实采样率 (Hz, 主函数唯一定义 = 160M)
assert(abs(fsAdc - chipOsr*prm.fs) < 1, ...
    'prm 口径不自洽: fsAdc(%gM) 必须 = chipOsr×fs(%d×%gM) —— 改降采样倍率请改主函数 chipOsr', ...
    fsAdc/1e6, chipOsr, prm.fs/1e6);
lpfFc     = 0.45 * fsAdc;  % 模拟截止 = 0.45×ADC率(chipOsr=1 时 = 0.45×fs, 同旧口径)
assert(lpfFc < 0.48*fsAdc, ...
    'lpfFc(%g MHz) 过分贴近/超过 ADC 奈奎斯特 %g MHz —— 抗混叠失效', lpfFc/1e6, fsAdc/2/1e6);
% 【2026-10 事故留档】单级 LPF 时代曾写死 24 MHz(当时以为 fs=60M), fs 改 40M 后
%   24M>奈奎斯特 → 20~24.6M 噪声抽取折叠回远端, 底噪反抬 +2.5 dB(实测复现)。
%   两段式后模拟截止挂 ADC 率(结构上永远远离信号带), 带边交给数字半带 —— 该类
%   配置错误从结构上消除。
% ---------- 5. TX 移相器实际实现（DDMA 相位量化 + 随机误差） ----------
enablePhaseQuant  = true;  % [开关] 相位量化到 N-bit 可实现值
phaseBits         = 7;     % [旋钮] 移相器位数（8-bit → 360/256 = 1.40625°/LSB）
enablePhaseErr    = true;  % [开关] 附加逐 chirp 正态随机相位误差
phaseErrSigma_deg = 1;     % [旋钮] 随机误差标准差 (deg rms，典型 1~5）

% ---------- 6. 天线/射频链带内幅度起伏（确定性 AM → R-FFT 展宽） ----------
enableRipple = true;       % [开关] 带内幅度起伏
ripplePP_dB  = 3;          % [旋钮] 起伏峰峰值 (dB/通道；逐通道形状相位见 W 权重区）

% ---------- 7. 扫频非线性（共锁相环 → 8TX 同非线性；残余调频∝距离） ----------
enableNonlin = true;       % [开关] 扫频非线性
nonlinPct    = 0.017;        % [旋钮] 非线性度 (%) = 最大瞬时频偏 / 扫频带宽 B
nonlinType   = 'hump';     % [旋钮] 形状: 'hump'单凸 | 'dip'单凹 | 'hump_dip'先凸再凹 | 'dip_hump'先凹再凸

% ---------- 8. 相位噪声 + PLL 环路滤波器（两项独立开关；环路滤波只整形 PN） ----------
enablePhaseNoise = true;   % [开关] 锁相环相位噪声（共锁相环 → 8TX 一致）
pnProfile        = [1e3 -49.25; 10e3 -65.84; 100e3 -70.80;...
                    1e6 -95.42; 10e6 -121.76; 100e6 -130.03];  % [旋钮] PN 拟合点 [Hz, dBc/Hz]
                           %   （近端 = 环路内特性；远端截断行为由本块环路滤波器决定）
enableLoopFilter = true;    % [开关] PLL 环路低通整形（仅 enablePhaseNoise=true 时有作用；
                           %   关闭 = 输出 PN 按拟合点外推、肩部不截断）
loopBW_Hz        = 800e3;  % [旋钮] 环路带宽 (Hz)：输出 PN 在此偏移处被陡峭截断，
                           %   所有目标肩部截断在同一偏移（与距离无关）
loopOrder        = 4;      % [旋钮] 环路滤波器阶数（截断陡度, PSD ~12n dB/oct）

% ---------- 8b. PLL 扫频跟踪误差（快斜坡下环路滞后 ∝ B/Ts → 肩膀随带宽抬高） ----------
% 芯片公司口径：chirp 带宽增大 → 相位噪声对 RD 图的恶化增大。机理 = 环路有限
% 带宽跟踪快扫频时, 每个 chirp 的实际频率轨迹存在随机滞后 δf, 其 σ ∝ 调频
% 斜率 S = B/Ts。与步进频跳频误差（第 1 块）同注入点同机理: δf 同进发射与本振
% → 差拍不变, 残余 = 逐 chirp 随机相位 2π·δf·τ（随距离增大）→ RD 图目标两侧
% 肩膀。与 PN 机理（φ(t)−φ(t−τ) 连续去相关）可叠加、互相独立。
% 与 pnProfile 的关系：pnProfile 是【当前标定带宽 B=214MHz 下】的拟合, 已把该
% 带宽下的跟踪误差吸收进去 —— 所以默认系数 = 0, 改 B 后的【额外】恶化量由
% 本项补足（一个实测点只能标定一个 B 的总和, 拆不开两项, 这是诚实口径）。
enablePllTrack = true;     % [开关] PLL 扫频跟踪误差（开关开 + 系数 0 = 行为不变）
pnTrackCoef    = 7.6e-9;   % [旋钮] 跟踪误差系数 (s)：σ_δf = 系数×B/Ts (Hz), 逐 chirp
                           %   独立正态抽取。换芯片/换带宽后实测肩膀比仿真高时, 由
                           %   超出量反推: 多普勒域非相参能量占比 ≈ σ_φ², 其中
                           %   σ_φ = 2π·σ_δf·τ (rad, τ=目标双程时延) →
                           %   系数 = σ_φ / (2π·τ·B/Ts), 再微调至肩膀对齐。
                           %   量级参考: σ_φ=0.05 rad @300m/B=214M → ≈ 7e-10 s。
                           %   【2026-10-29 两轮重标定(用户实测口径)】原 1e-8 时 300m 线损
                           %   −4.8dB, 实测 SNR 损耗随距离变化没这么大 → 8.9e-9(−3.8dB)
                           %   仍偏高 → 7.6e-9, 300m 线损 −2.8dB(σφ≈0.80 rad)。

% ---------- 9. 零频泄漏（片内直通 τ=0 → 差拍=0 → bin0；独立开关） ----------
enableLeakDC   = true;    % [开关] 恰直流泄漏路径。开启后 bin0 出码梳，且 bin0 行
                           %   慢时间序列为实数 → 实采样多普勒共轭镜像（实测 bin0
                           %   完全看不到 → 默认关；做镜像机理演示时再开）。
                           %   幅度基准 = 隔离度矩阵 isoMat（下方 TX-RX 收发隔离
                           %   泄漏块）× leakRes_dB + leakDC_dB；需总开关同时开
leakDC_dB      = 0;        % [旋钮] 直流路径相对电平 (dB, 与 leakRes_dB 相加后作用于
                           %   隔离度基准幅度；片内直通通常最强, 0 = 满值）

% ---------- 10. 近距离非 0 bin 泄漏多径族（地弹+屏蔽腔多次反射 / 环境残余反射） ----------
enableLeakNear = true;     % [开关] 近距非零 bin 泄漏多径族【统一开关】——
                           %   一次性启停 leakPaths 全部近距路径；false → 近距码梳
                           %   全部消失，只剩闪烁噪声抬底（对比实验用）。
                           %   幅度基准同零频块：isoMat × (leakRes_dB + 逐路径电平)
leakPaths     = [ 58.6e3,  -20;  % [旋钮] 近距泄漏多径表：每行一条路径 [差拍频率Hz, 相对电平dB]。
                  117.2e3, -48;  %   差拍频率 fb = 路径延迟 τ×slope → 视在距离 = c·fb/(2·slope)，
                  5.86e3,  -50 ];%   物理固定、与采样配置解耦（bin 口径会随 fs/采样点数漂移，
                           %   Hz 不会 —— 照 96 修法 Hz 化, 表值 = 原 1.5/3.0/0.15 bin
                           %   @fs40M/1024 点等值换算）。
                           %   相对电平 = 逐路径修正，与全局 leakRes_dB 相加后
                           %   作用于隔离度基准幅度（0 = 满值）。
                           %   默认三行 = 三类物理机制（每条独立固定随机相位）：
                           %     58.6 kHz : 地弹 + 屏蔽腔内多次反射（视在 ~1.05 m, 主残余, 已标定）
                           %     117.2 kHz: 环境残余反射（视在 ~2.1 m, 暗室吸收材料残余 0.5~2 m）
                           %     5.86 kHz : 天线互耦（视在 ~0.1 m → 亚 bin, 多数落 bin0 被丢弃）
                           %   自由增删行 = 加/减机制；恰直流（fb=0）用上方『零频泄漏』块独立控制
leakRes_dB    = -5;        % [旋钮] 泄漏残余电平修正 (dB)：【近距码梳高度直接控制量】，
                           %   ±X dB → 全部码梳峰值同步 ±X dB（作用于 leakPaths 所有
                           %   路径与零频路径的泄漏幅度，与滤波器形状无关；逐路径
                           %   差异用表中第二列相对电平）。物理含义：泄漏路径的净
                           %   损耗/绕射损耗；0 = 按隔离度满值。标定：-5 与实测对齐

% ---------- 11. TX-RX 收发隔离泄漏【总开关】（隔离度矩阵 + PN裙边 + 闪烁噪声） ----------
enableLeakage = true;      % [开关] TX-RX 收发隔离泄漏总开关 = 零频路径 + 近距多径族 +
                           %   闪烁噪声三者随它整体启停（上面两块的细分开关需与
                           %   本开关同时开才有效）
isoMat        = -70*ones(8,8);  % [旋钮] 8×8 TX→RX 隔离度 (dB)，可任意改
% isoMat(5,7)   = -35;       %   TX1→RX4 隔离度差（示例）
% isoMat(6,7)   = -35;       %   TX2→RX4 隔离度差（示例）
leakPnTau     = 0.5e-9;    % [旋钮] 泄漏 PN 等效去相关延迟 (s)，板级 ps~ns 级；
flickerK      = 2e-11;     % [旋钮] 闪烁噪声系数: S(f) = P_leak·K/f^flickerAlpha
flickerAlpha  = 0.5;       % [旋钮] 闪烁频谱指数（0.5 匹配实测：RX4 近距 ≈+10 dB 缓降，
                           %   近最大测距才与热噪合并；其余 RX ≪ 白噪声不可见）

% ---------- 12. TX-TX 发发隔离泄漏 + PA 三阶交调（两项独立开关，交调依赖耦合） ----------
enableTxCoupling = true;     % [开关] 发发隔离泄漏
ttCorr_dB        = 0;      % [旋钮] 修正值 (dB)：串入信号经后级其它电路（如 PA）的净修正
ttIsoMat = [ NaN  -47 -39 -52 -44 -49 -41 -53;   % [旋钮] 8×8 发发隔离度 (dB)，m 行→n 列，
             -45  NaN -50 -42 -48 -37 -53 -46;   %   m≠n 取值于 -35 ~ -55 dB，对角 NaN 不用
             -38 -51  NaN -46 -54 -43 -40 -49;
             -50 -43 -47  NaN -39 -52 -45 -41;
             -46 -54 -40 -48  NaN -42 -51 -38;
             -42 -48 -53 -39 -50  NaN -46 -44;
             -53 -41 -46 -55 -43 -49  NaN -40;
             -48 -39 -51 -44 -37 -46 -42  NaN ];
enablePaIM3 = true;        % [开关] PA 三阶交调（与发发耦合同块：无耦合时多码复合
                           %   不存在 → 无空子带交调信号；交调作用于耦合后的信号）
paIM3_dBc   = -4;        % [旋钮] 双音 IM3 电平 (dBc, 负值 = IM3 低于载波 X dB, 典型 −20 ~ −5)。
                           %   空子带泄露 ≈ 耦合度 × 10^(paIM3_dBc/20), 线性: ±1 dB → 泄露 ±1 dB

% ---------- 13. RX-RX 收收隔离泄漏（RXm 串入 RXn 前端, 线性叠加） ----------
enableRxCoupling = true;   % [开关] 收收隔离泄漏
rrCorr_dB        = 0;      % [旋钮] 修正值 (dB)（RX 间影响 ≈ 隔离度，默认 0）
rrIsoMat = [ NaN  -44 -52 -38 -49 -43 -54 -46;   % [旋钮] 8×8 收收隔离度 (dB)，m 行→n 列，
             -51  NaN -46 -53 -40 -48 -44 -39;   %   m≠n 取值于 -35 ~ -55 dB，对角 NaN 不用
             -39 -47  NaN -44 -52 -41 -49 -53;
             -53 -41 -48  NaN -45 -55 -42 -47;
             -46 -54 -43 -50  NaN -39 -52 -44;
             -42 -49 -55 -43 -47  NaN -46 -51;
             -50 -38 -45 -54 -41 -53  NaN -43;
             -47 -52 -40 -46 -49 -42 -55  NaN ];

% ---------- 14. ADC 量程与接收增益（时域限 ±32768 counts; 不模拟整数化） ----------
adcBits        = 16;        % [旋钮] ADC 位数(12/14/16...); 数据恒以 2 字节 int16 输出,
                            %   不足 16bit 时左移 (16−bits) 位补齐(行业惯例)。
                            %   只在本模块指定, 主函数不经 prm 传入
adcFullScale_V = 0.6;       % [旋钮] ADC 满量程电压 (V, 正半量程口径: ±FS ↔ ±32768 counts,
                            %   即差分峰峰 = 2×FS = 0.5 Vpp, 车载雷达 ADC 常见范围低端)。
                            %   数据缩放关系: counts = 电压 / (FS/2^(bits-1)), 1 LSB = FS/2^15
                            %   = 7.6 μV。信号与热噪在链路里是【电压】, 不随 FS 变 → 它们的
                            %   counts ∝ 1/FS(FS 减半 counts 加倍); 但 ADC 芯片噪声的 counts
                            %   与 FS 无关(手册 SNR 是满幅正弦口径, 折算时 FS 约掉)→ 改 FS 会
                            %   微调热噪/芯片噪配比, SNR 随之微变, 不是纯数值缩放。用途: 仿真
                            %   rawdata 与实测 counts 标尺对齐(SNR 已吻合的场合); 若只想同比
                            %   放大信号+全部噪声(SNR 完全不动), 用数字增益 digGainDb。
% 接收链为【固定增益】（无自动增益），增益 G_rx_dB 在"雷达参数"区。标定口径：
% 无强目标时（弱目标+收发隔离泄漏+底噪）时域约 ±200 counts（与实测一致）；
% 强目标最多允许逼近满量程，时域超出 ±32768 = 硬件削顶（运行时告警提示）。
% 2D-FFT 之后的数值可远超 ±32768（FFT 相干积累增益），不属超限。
enableInt16Wrap = true;   % [开关] 16bit 溢出异常模拟【本文件新增】—— 超出
                           %   ±32768 的样点按【补码环绕】处理（+32768→−32768,
                           %   +40000→−25536）, 不是削顶：模拟硬件 ADC 输出 16 位
                           %   有符号整型、后级（FPGA 寄存器/驱动）按补码直接解析
                           %   而未做饱和 → 高位溢出数值翻到对面, 时域"翻折"、
                           %   频域谐波+宽带杂散 → RD 底噪抬升。
                           %   false = 理想（不环绕, 超限仅告警）
enableInt16Quant = true;  % [开关] ADC 整数量化【本文件新增功能】—— 真实雷达
                           %   ADC 输出 2 字节整数(int16): 数据四舍五入到 ADC 的
                           %   ±1 LSB 网格(2^(bits-1) 层), 再按行业惯例【左移
                           %   (16−bits) 位】装入 int16(占满容器, 12bit→×16)。
                           %   false = 理想(counts 连续值, 不量化, ADC-count 口径)
adcChipSnrDb = 60.2;        % [旋钮] 芯片手册 SNR (dB, 满幅正弦口径; 供应商给的值)。
                            %   实测核对结论: 按 60.2 建模时近直流底噪只比平坦段低 ~3 dB,
                            %   实测低 ~28 dB → 实际 ADC 噪声比手册口径低 ~28 dB(可能是
                            %   最差规格/其它带宽口径, 或片内 ADC 后还有数字直流消除)。
                            %   因此默认不用本值, 改用下方实测标定电压; 把标定值置空
                            %   (adcChipSigmaV_cal = [])即退回本手册口径
adcChipSigmaV_cal = 25.9e-6; % [旋钮] ADC 噪声实测标定 (V rms, 输入参考, 优先于手册 SNR)
                            %   来源: 实测近直流坑深 −28 dB —— 坑底 = 热噪被两级高通压
                            %   死后露出的 ADC 噪声 → σ_chip = σ_热噪(V)×10^(−28/20)
                            %   = 650.5 μV × 0.0398 ≈ 26 μV(按 G_rx=54 dB 的热噪电压)。
                            %   电压口径与满量程无关: 改 adcFullScale_V 坑深自动保持;
                            %   改 G_rx 坑深随热噪同向变化(物理正确)。等效 SNR ≈ 88.7 dB
digGainDb    = 0;           % [旋钮] 数字增益(dB) —— 纯数值缩放: 在芯片噪声之后、
                            %   量化/溢出检查之前把整个数据乘 10^(dB/20), 信号与
                            %   全部噪声(热噪+芯片噪)同比例, SNR 与噪声占比都不变。
                            %   用途: 仿真 rawdata 与实测对齐幅度标尺(实测 counts
                            %   偏大且 SNR 一致时, 用它对齐; 若实测 SNR 也更高,
                            %   应改 G_rx_dB 而不是它)。对应真实芯片 DFE/位宽
                            %   打包一类的数字增益级, 默认 0 dB 不生效
enableAdcChipNoise = true;  % [开关] ADC 芯片自身噪声—— 采样保持/
                           %   缓冲器输入参考白噪, 电压域 σ=(满量程/√2)/10^(SNR/20),
                           %   加在量化前。与 adcBits 无关(物理噪声不随位数变);
                           %   勿用"减小 adcBits"代替 —— 那是换一颗理想低位 ADC

% ---------- 15. 馈线长度色散（TX/RX 16 条馈线不等长 → 通道相位随频率失配） ----------
enableFeedDisp = true;    % [开关] 馈线长度色散【本文件新增功能】
feedLenScale   = 0.015;     % [旋钮] 最大长度差距 (m)：馈线长度差矩阵中 |最大值|=1
feedLenDev = [ 0.30 -0.50  0.80 -0.20  1.00 -0.70  0.10  0.45 ; ...
              -0.85  0.25 -0.40  0.60 -0.15  0.95 -0.55  0.05 ];  % [旋钮] 16 条馈线的
                           %   相对长度差（无量纲形状, |最大值|=1, 可任意改）：
                           %   前 8 个 = TX1~TX8, 后 8 个 = RX1~RX8；
                           %   实际长度差 = feedLenScale × 本值 (m)
                           %   物理含义：馈线不等长 → 各通道相位差随频率线性变化
                           %   φ(f) = 2π·(f − f_cal)·ΔL/c, f_cal = 全带中心频率
                           %  【已假设校准】：中心频率处 16 通道相位相等（虽然
                           %   不完全符合实际, 但等效于出厂做了中心频点校准），
                           %   其它频率的相位差按线长差计算

%% 常数
c = physconst('lightspeed');

%% ============================ 雷达参数 ============================
% 共享雷达常数全部由主函数经 prm 结构体注入(唯一定义处在主函数, 此处不再重复):
startFreq     = prm.fStart;        % 扫频起始频率
adcSampleRate = prm.fs;        % 输出数据率 fs（降采样后, 单路实数; ADC 物理率 = prm.fsAdc）
% nSamples / nChirps 由主函数注入(函数入参)

nSamplesValid = nSamples/2 + 1;% 实采样距离维有效 bin 数（0 ~ fs/2，含 fs/2 点）

FFT_N0 = 128;                  % 角度二维 CFAR FFT 点数（预留）
% FFT_d0 阵列实际间距由主函数注入(antenna_t/antenna_r 已按它缩放), 此处不重复定义

B      = prm.B;                % Chirp 有效带宽（= 采样窗 Ts 内扫过的带宽）
Ts     = nSamples/adcSampleRate;   % ADC 采样窗 = 25.6 us
T_idle = prm.Tc - Ts;          % chirp 间保护时间
Tc     = prm.Tc;               % chirp 周期 = 37 us（必须大于 ADC 采样窗 Ts！）
adcStartDelay = 5e-6;          % ADC 相对 chirp 起始的采样延迟 (s)：真实链路不在 chirp
                               % 起点采样、首部留几 µs —— 泄漏单音在 chirp 边界的相位
                               % 跳变（DDMA 码步进+步进频跳变）会激发高通建立瞬态
                               %（幅度可达原始泄漏电平），ADC 窗避开它。≤ T_idle；
                               % 设 0 恢复从 chirp 起点采样（会看到首/末端瞬态伪值）
PRF    = 1/Tc;                 % chirp 重复频率

FreqSlope = B/Ts;              % 调频斜率 ≈ 8.36 MHz/us

DeltaR = c/(2*B);              % 距离分辨率 = 0.701 m（只由带宽 B 决定，与采样方式无关）
Rmax   = c*adcSampleRate/(4*FreqSlope);   % 实采样最大不模糊距离 ≈ 359 m

% IF 抗混叠低通滤波器（位于 dechirp 之后、ADC 之前）：参数与使能已移至
% 文件前部"仿真控制开关"区（enableLpf / lpfOrder / lpfFc）。只负责抗混叠
% 与带外压制，与 PN 肩部截断无关（肩部截断由 PLL 环路滤波器产生，
% 见开关区 enableLoopFilter / loopBW_Hz / loopOrder）。通带上缘 8.5 MHz → 有效不模
% 糊距离 ≈ c·8.5e6/(2·slope) ≈ 254 m。

% 发射功率与天线增益（雷达方程用，8 个通道均相同）
Pt_dBm        = 13;            % 每个 TX 发射功率 (dBm)
radomeLossDb  = 0;             % [旋钮] 天线罩损耗 (dB, 单程): TX 辐射与 RX 接收各过一次罩 → 双程共 2 dB
rfChainLossDb = 0;             % [旋钮] 射频链路损耗 (dB): RX 天线口到接收机前的插入损耗(位于低噪放之前, 劣化 SNR)
Gt_dBi = 15 - radomeLossDb;    % 每个 TX 天线增益 (dBi) = 15 − 天线罩损耗(发射路径过罩)
Gr_dBi = 15 - radomeLossDb - rfChainLossDb;  % 每个 RX 天线增益 (dBi) = 15 − 天线罩(接收路径过罩) − 射频链路
Pt = 10^(Pt_dBm/10)*1e-3;      % 发射功率 (W)
Gt = 10^(Gt_dBi/10);           % TX 天线增益（线性）
Gr = 10^(Gr_dBi/10);           % RX 天线增益（线性）
% 注: 两个损耗折进 Gt/Gr → 雷达方程自动生效(逐目标 Pr 打印与回波幅度 ampEq 都用 Pt·Gt·Gr):
%     双程合计压低回波 2×radomeLossDb + rfChainLossDb = 4 dB, 绝对幅度与输入SNR 同降;
%     G_rx 只是 ADC 电压刻度(信号/噪声同乘), 不承载这两个损耗 —— 损耗必须加在方程侧才会劣化 SNR。

% RX 链路增益：把 RX 输入端功率折算成 ADC 域电压幅度 v = sqrt(P·G_rx)。
% 同时作用于信号与噪声，因此不影响信噪比，只设定 ADC 电压刻度（counts 数）。
% 固定物理增益（无自动增益）。标定：54 dB 时无强目标场景（仅底噪/弱目标/
% 泄漏残余）时域峰值 ≈ ±200 counts（与实测对齐）；每 +6 dB 时域 counts ×2。
G_rx_dB = 54;                  % RX 链路总增益 (dB)
G_rx    = 10^(G_rx_dB/10);     % 线性

if ~isfield(prm, 'stepDir') || isempty(prm.stepDir), prm.stepDir = 1; end
stepDir = prm.stepDir;          % 步进频方向: +1 上跳频(每 chirp 起始频率更高) / −1 下跳频(更低)
f_step = stepDir * c/2/(DeltaR*(nChirps-1)) * coupling_ratio;  % chirp 间起始频移(带符号) = ±coupling_ratio·B/(nChirps-1)

f_center = startFreq + B/2;              % 中心载频
f_center = f_center + f_step*(nChirps-1)/2;
lambda   = c/f_center;                   % 雷达工作波长

Vres = c/(2*f_center*(Tc*nChirps));      % 速度分辨率
Vmax = c/f_center/4/(Tc);                % 最大不模糊速度 ±Vmax

nTX       = size(antenna_t, 2);    % TX 通道数(由主函数注入的阵列读出, 不硬编码)
nRX       = size(antenna_r, 2);    % RX 通道数(同上)
nGroups   = prm.nGroups;           % DDMA 分组数(主函数 = numel(encodeSeq), 随编码序列联动)
groupSize = nChirps/nGroups;       % 每组 chirp 数

%% ============================ 目标参数 ============================
% 目标由主函数以【结构体数组】注入(targets): 每个目标一个元素, 字段:
%   .gate_R  距离 gate(bin)   .gate_D  法向速度 gate(bin, 产生多普勒; 可超 512 折叠同 bin)
%   .azi     方位角(deg)       .ele     俯仰角(deg)
%   .rcsDb   RCS(dBsm)         .isCorner  是否角反(校准/检查模块用)
%   .vTan    切向速度 (m/s, 水平面内垂直视线方向; 雷达测不到, 跟踪算法负责推测)
%   【切向速度的物理】v⃗ = vRad·û + vTan·t̂ (û=视线单位矢量, t̂=水平面内垂直 û),
%   位置 p(t) = p0 + v⃗·t → 精确距离 |p(t)| = sqrt(R0² + 2R0·vRad·t + (vRad²+vTan²)t²)
%   (p0·t̂=0 → 切向只进二阶项: 不产生多普勒, 但帧内横向走动改变视角 → 阵列相位缓变/角度微移)
Obj_amount = numel(targets);
gate_R = [targets.gate_R];
gate_D = [targets.gate_D];
Obj_Angle_azi = [targets.azi];
Obj_Angle_ele = [targets.ele];
RCS_dB  = [targets.rcsDb];
RCS     = 10.^(RCS_dB/10);      % 线性 σ (m²), 功率计算用(始终实数)
% 复数 RCS(目标回波初相, 工程含义: 反射系数 σ∠φ 随目标形状/材质/姿态变化,
% 两个同距同速目标即使 RCS 幅度相同, 初相不同 → 阵列上的干涉图样完全不同)
if isfield(targets, 'rcsPhDeg')
    RCS_c = RCS .* exp(1i*deg2rad([targets.rcsPhDeg]));   % 复数 σ(幅度方程用)
else
    RCS_c = RCS;                                          % 无该字段 → 全部同相(原有行为)
end
Obj_R   = DeltaR*gate_R;                 % 各目标距离 (m)
Obj_v   = Vres*gate_D;                   % 各目标法向速度 (m/s, v>0=远离)
if isfield(targets, 'vTan'), Obj_vTan = [targets.vTan]; else, Obj_vTan = zeros(1, Obj_amount); end
Obj_v2  = Obj_v.^2 + Obj_vTan.^2;                % 全速度模平方(精确距离式的二阶项系数)

% 由距离和角度计算目标笛卡尔坐标 + 法向/切向单位矢量(速度分解用)
Obj_uHat = [cosd(Obj_Angle_ele) .* sind(Obj_Angle_azi); ...
            cosd(Obj_Angle_ele) .* cosd(Obj_Angle_azi); ...
            sind(Obj_Angle_ele)];                 % û 视线单位矢量 [3×K]
Obj_tHat = [cosd(Obj_Angle_azi); ...
            -sind(Obj_Angle_azi); ...
            zeros(1, Obj_amount)];                % t̂ 水平切向单位矢量(⟂û) [3×K]
Obj_xy = Obj_R .* Obj_uHat;                       % 初始位置 [3×K]

% ---- 接收机白噪声（与目标无关的空暗室底噪） ----
% 现实中底噪在无目标的暗室中标定：Pn = k·T0·B_NEB·NF，折算到 RX 输入端，
% 与场景中有什么目标完全无关。B_NEB 为 IF 低通等效噪声带宽（巴特沃斯）。
% 链路增益 G_rx 把功率折算为 ADC 域电压幅度 v = sqrt(P·G_rx)（信号/噪声
% 同一刻度，G_rx 不影响信噪比）。白噪声被 LPF 限带后经 ADC 采样，
% 各样点为独立实高斯，标准差 noiseSigma = sqrt(Pn·G_rx)。
k_B     = physconst('Boltzmann');
T0      = 273.15 + 80*1;                                        % 实际工作噪声温度 (K)
% ---- 数字抽取链设计(chipOsr≥2): nDStg 级半带 FIR, 每级截止 = 该域 fs/4 ----
% 半带构造: h[n] = 0.5·sinc(n/2)·kaiser(L,β) —— 偶数点(除中心)自动为 0(真半带
% 结构: 一半系数为零、对称线性相位), 归一化 DC 增益 = 1。末级截止恰为 fs/2 = 信号带边。
% hbMagOf(s,f) = 第 s 级在物理频率 f (Hz) 处的幅度响应(供目标表/曲线图/filtGain 用)
nDStg = round(log2(chipOsr));                        % 数字抽取级数(0/1/2)
hbH = cell(1, max(nDStg,1));  hbFc = zeros(1, max(nDStg,1));  hbR = zeros(1, max(nDStg,1));
for s = 1:max(nDStg,1)
    R_s  = fsAdc / 2^(s-1);                          % 该级工作域采样率 (Hz)
    L    = hbTaps(min(s, numel(hbTaps)));            % 该级抽头数
    nnS  = -(L-1)/2 : (L-1)/2;
    h    = 0.5*sinc(nnS/2) .* kaiser(L, hbBeta).';   % sinc 偶点自动为零 → 真半带
    hbH{s} = h / sum(h);                             % DC 增益归一
    hbFc(s) = R_s/4;  hbR(s) = R_s;
end
nnOf = @(h) -(numel(h)-1)/2 : (numel(h)-1)/2;   % FIR 对称中心轴(行向量)
hbMagOf = @(sI, fHz) arrayfun(@(fq) abs(sum(hbH{sI} .* exp(-1j*2*pi*(fq/hbR(sI))*nnOf(hbH{sI})))), fHz);

% ---- 带边噪声均衡设计(照 96 修法 enableCicComp; 在 B_noise 之前 —— NEB 需计入均衡) ----
% 输出域白噪响应 = 两级 ÷2 折叠后的【折叠对合计】H_pair(f) = sqrt(Σ_k [H²(k·fs+f) + H²((k+1)·fs−f)]),
% k = 0..chipOsr−1(HB 链各级阻带已把远端 alias 压死, 主要贡献 = 末级过渡带对 f 与 fs−f)。
% 不均衡时 H_pair 在末级过渡带(≈0.93~1×fs/2)下垂 → 噪底出现滤波形状(实测没有)。
% 均衡目标 = H_pair 本身(不是带内幅度!): G = 1/max(H_pair, cap) —— 噪底按构造拉平;
% 96 的"全额均衡抬折叠噪声 +2dB"教训是均衡带内幅度的做法, 对折叠对合计均衡不发生。
% 信号得 G×H(f): 350m(v=0.975) 残余 droop ≈ −0.2dB, 奈奎斯特点 ≈ −3.4dB(未用带尾);
% 带外强目标折叠残余同乘 G(峰值 +2.7dB, 只影响最后 ~10 bin, "鬼影 −4~−7dB 折回"口径不变)。
if nDStg > 0
    gAna = @(fb) (1-enableLpf) + enableLpf*(1./sqrt(1 + (fb/lpfFc).^(2*lpfOrder)));   % 模拟 LPF 幅度
    gDig = @(fb) arrayfun(@(fq) prod(arrayfun(@(sI) hbMagOf(sI, fq), 1:nDStg)), fb);  % 半带链幅度
    gChainA = @(fb) gAna(fb) .* gDig(fb);                    % ADC 域原生频率的链路幅度
    vG = linspace(0, 1, 513).';  fG = vG * prm.fs/2;         % 输出频率栅格(0 ~ fs/2)
    Hpair = zeros(size(fG));
    for kA = 0:chipOsr-1                                     % 折叠到输出 bin 的全部原生频点
        Hpair = Hpair + gChainA(kA*prm.fs + fG).^2 + gChainA((kA+1)*prm.fs - fG).^2;
    end
    Hpair = sqrt(Hpair);
    if enableHbComp
        capV  = interp1(vG, Hpair, compFreezeV);             % 冻结点处的合计幅度(截断增益上限)
        bComp = fir2(hbCompTaps-1, vG, 1./max(Hpair, capV)); % 逆 droop(阶数 = 抽头−1, 奇抽头
        bComp = bComp/sum(bComp);                            %   Type I → 奈奎斯特点增益可非零; 直流增益 1)
        gC    = @(fb) arrayfun(@(fq) abs(sum(bComp .* exp(-1j*2*pi*(fq/prm.fs)*(0:numel(bComp)-1)))), fb);  % 均衡幅度(DTFT 直算, 标量/向量都吃)
    else
        bComp = 1;  gC = @(fb) ones(size(fb));               % 关 → 全通(=旧行为, 噪底带下垂)
    end
    % 输出噪声带宽 = (折叠对合计 × 均衡)²的数值积分 —— 抽取只搬谱不改总功率, 折叠
    % 贡献天然含在 H_pair 里; 均衡后 ≈ fs/2(拉平的噪底), 不均衡 ≈ 末级半带 NEB
    B_noise = trapz(fG, (Hpair .* gC(fG)).^2);
else
    bComp = 1;  gC = @(fb) ones(size(fb));                   % 无数字链 → 无 droop, 全通
    B_noise = lpfFc * (pi/(2*lpfOrder)) / sin(pi/(2*lpfOrder));  % 巴特沃斯 NEB ≈ 1.02*fc
end
Pn_W    = k_B * T0 * B_noise * 10^(NF_dB/10);         % 单通道噪声功率 (W, RX输入端)
Pn_dBm  = 10*log10(Pn_W*1e3);                         % (dBm)
noiseSigma = sqrt(Pn_W * G_rx);                       % ADC 域噪声标准差 (链路后, 等效电压)
% 链路【前】的白噪声标准差：白噪声在滤波链前加入、经完整 IF 链后
% 总功率收敛到 noiseSigma²（NEB 标定），近距 bin 底噪被高通压低几十 dB。
noiseSigmaW = noiseSigma * sqrt(adcSampleRate/(2*B_noise));

% 【探针钩子】可选第 10 参 swOv: 结构体覆盖任意【定义在上方】的开关/旋钮
%   (单因素对比实验用, 见 quantProbe.m; 日常不传 → 行为与原版完全一致)。
%   位置在参数区末尾、数据计算(W 加权)之前; 覆盖后重算受牵连的派生量
%   (G_rx 线性值、噪声 σ) —— 否则只改 G_rx_dB 不会传导到噪声。
%   注意: chipOsr 已上提主函数(经 prm 注入)、hbTaps/enableHbComp 派生链太长
%   (fsAdc/lpfFc/半带/均衡/NEB/osf 全要重算), 均不经 swOv 覆盖 —— 两段式对比
%   实验请直接改主函数的 chipOsr(射频核自动跟随)或本文件 hbTaps/均衡开关。
% 可选第 11 参 printCfg(默认 true): false 时只静默"雷达参数"配置块与配置类
%   状态行(算力/芯片噪声σ/量化), 目标表与量程检查等数据行仍打印 ——
%   两遍射频同配置时, 配置只从第一遍(或目标遍)打印一次, 不重复刷屏。
if nargin < 11 || isempty(printCfg), printCfg = true; end
if nargin >= 10 && ~isempty(swOv)
    fnOv = fieldnames(swOv);
    for qOv = 1:numel(fnOv)
        eval([fnOv{qOv} ' = swOv.(fnOv{qOv});']);    %#ok<EV2I> 覆盖已定义的开关/旋钮
    end
    G_rx = 10^(G_rx_dB/10);
    noiseSigma = sqrt(Pn_W * G_rx);
    noiseSigmaW = noiseSigma * sqrt(adcSampleRate/(2*B_noise));
end

%% ======================= 参数合法性检查与打印 =======================
assert(Tc > Ts, 'chirp 周期 Tc 必须大于 ADC 采样窗 Ts，否则chirp相互重叠！');
assert(adcStartDelay <= T_idle, 'adcStartDelay 必须小于等于 chirp 间保护时间 T_idle！');
% 注意：不再禁止目标超过 Rmax —— 超量程目标由模拟 LPF 压制，残余混叠成鬼影

onoffStr = @(b) strrep(strrep(char(string(b)), 'true', '开'), 'false', '关');  % 开关状态转文字

if printCfg
% ---- 配置块只随第一遍射频打印(目标表/量程检查等数据行在块外, 每遍都打) ----
fprintf('================ 雷达参数 ================\n');
fprintf('调频斜率       : %.3f MHz/us\n', FreqSlope/1e12);
dirTxt = '上跳频(+1, 每 chirp 起始频率更高)';  if stepDir < 0, dirTxt = '下跳频(−1, 每 chirp 起始频率更低)'; end
% 帧射频范围随方向: 上跳 [f0, f0+(1+cr)B]; 下跳 [f0−crB, f0+B] —— 总带宽同为 (1+cr)·B
fFrameLo = startFreq + min(0, f_step*(nChirps-1));            % 全帧最低射频
fFrameHi = startFreq + B + max(0, f_step*(nChirps-1));        % 全帧最高射频
fprintf('步进频         : %s, 占比 %.2f → 帧射频范围 %.4f ~ %.4f GHz (总带宽 %.1f MHz), 末chirp起点 %.6f GHz\n', ...
        dirTxt, coupling_ratio, fFrameLo/1e9, fFrameHi/1e9, ...
        (1+coupling_ratio)*B/1e6, (startFreq + f_step*(nChirps-1))/1e9);
% 波形频率三要素(含扫频+跳频全跨度): 全帧最低 / 全帧中点 / 全帧最高
fprintf('波形频率        : 起始 %.4f GHz → 中心 %.4f GHz → 结束 %.4f GHz (含扫频+跳频全跨度)\n', ...
        fFrameLo/1e9, f_center/1e9, fFrameHi/1e9);
% 俯仰无模糊范围: 俯仰虚拟行间距(数据实测) × 中心频率波长 → 栅瓣-free 俯仰角范围
% (坐标按 1 nm 量化后再去重: 居中平移使 z 为半整数栅格, 浮点求和尾数会破坏 unique)
% 俯仰/方位基本间距 = 相邻间距的最大公约数(spacingGCD, 非 min —— 不规则布阵下
% min 只是最近的一对, gcd 才是栅格基线, 例如 0/3/5 → 间距 3,2 → 基线 1)
dzRow   = spacingGCD(antenna_t(3,:).' + antenna_r(3,:));   % 俯仰行基线间距 (m)
lamC    = c/f_center;                                  % 中心频率波长 (m)
eleUa   = asind(min(1, lamC/(2*dzRow)));               % 俯仰无模糊范围 (≥1 时为 ±90)
dzAzi   = spacingGCD(antenna_t(1,:).' + antenna_r(1,:));   % 方位基线间距 (m)
aziUa   = asind(min(1, lamC/(2*dzAzi)));
fprintf('俯仰无模糊范围  : ±%.1f° (俯仰行间距 %.3f mm = %.2fλ, 中心频率 %.4f GHz; 超出将出现俯仰栅瓣)\n', ...
        eleUa, dzRow*1e3, dzRow/lamC, f_center/1e9);
fprintf('方位无模糊范围  : ±%.1f° (方位最小间距 %.3f mm = %.2fλ)\n', ...
        aziUa, dzAzi*1e3, dzAzi/lamC);
fprintf('跳频误差       : %s, σ=%.0f Hz/chirp 正态分布（逐 chirp 独立、不累积; 残余随机相位 ∝ 时延 → 远目标多普勒去相关）\n', ...
        onoffStr(enableHopErr), hopErr_Hz);
fprintf('馈线色散       : %s, 最大长度差 %.1f cm（16 条 TX/RX 馈线形状 × 缩放因子）, 全带中心 %.4f GHz 处校准（中心相位相等, 偏离中心按线长差计算）\n', ...
        onoffStr(enableFeedDisp), 100*feedLenScale, f_center/1e9);
fprintf('输出数据率     : %.1f MHz (ADC 实采 %.0fM ÷ 降采样倍率 chipOsr=%d, 倍率已上提主函数), 奈奎斯特带宽 [0, %.1f) MHz\n', ...
        adcSampleRate/1e6, fsAdc/1e6, chipOsr, adcSampleRate/2/1e6);
fprintf('模拟IF低通     : %s, 巴特沃斯 %d 阶, 截止 %.1f MHz(=0.45×ADC率%gM) | 响应: ', ...
        onoffStr(enableLpf), lpfOrder, lpfFc/1e6, fsAdc/1e6);
for fq = (prm.fs/2)/1e6*[0.5 0.7 0.85 0.95 1 1.25 1.5]   % 打印频点随 fs 联动(×奈奎斯特, MHz)
    hq = 1/sqrt(1+(fq*1e6/lpfFc).^(2*lpfOrder));
    fprintf('%.4gM:%.1f ', fq, 20*log10(hq));
end
fprintf('dB (1.0×=输出奈奎斯特 %g MHz; 模拟只管宽处, 带边见数字链)\n', prm.fs/2/1e6);
fprintf('数字抽取链     : chipOsr=%d → ADC@%gM, %d 级半带 FIR', chipOsr, fsAdc/1e6, nDStg);
for s = 1:nDStg, fprintf(' | HB%d %d抽头@%gM域/截止%gM', s, numel(hbH{s}), hbR(s)/1e6, hbFc(s)/1e6); end
fprintf(' | 带边噪声均衡 %s(%d抽头@输出域, 冻结 v=%.2f → 噪底全带拉平对齐芯片实测) | 输出NEB=%.2fM\n', ...
        onoffStr(enableHbComp && nDStg > 0), numel(bComp), compFreezeV, B_noise/1e6);
fprintf('模拟高通(两级) : %s, HPF1 %.0f kHz/%d阶 + HPF2 %.0f kHz/%d阶, 阻带残余 %.0f dB（=近距泄漏码梳高度, ±5dB→码梳±5dB）\n', ...
        onoffStr(enableHpf), hpf1Fc/1e3, hpf1Order, hpf2Fc/1e3, hpf2Order, hpfFloor_dB);
noiseState = strrep(char(string(enableNoise)), 'true', '开启');
noiseState = strrep(noiseState, 'false', '关闭（RD 图为无噪理想情况）');
fprintf('接收机白噪声   : %s, NF = %.0f dB, 等效噪声带宽 %.2f MHz, 空暗室底噪 Pn = %.1f dBm (RX输入端), ADC域σ = %.3f mV ≈ %.0f counts\n', ...
        noiseState, NF_dB, B_noise/1e6, Pn_dBm, noiseSigma*1e3, ...
        noiseSigma/(adcFullScale_V/2^(adcBits-1)));
fprintf('TX 移相器      : %d-bit (%.4f°/LSB), 量化模拟: %s, 随机误差模拟: %s (σ = %.1f°)\n', ...
        phaseBits, 360/2^phaseBits, onoffStr(enablePhaseQuant), onoffStr(enablePhaseErr), phaseErrSigma_deg);
fprintf('带内起伏       : %s, 峰峰值 %.1f dB/通道, 回波扫频RF域确定响应 → chirp内调幅 → R-FFT展宽\n', ...
        onoffStr(enableRipple), ripplePP_dB);
fprintf('扫频非线性     : %s, %.1f%%, 形状 %s (最大频偏 %.2f MHz, 共锁相环 → 残余调频∝距离)\n', ...
        onoffStr(enableNonlin), nonlinPct, nonlinType, nonlinPct/100*B/1e6);
fprintf('相位噪声       : %s, 拟合点: ', onoffStr(enablePhaseNoise));
fprintf('%.0fHz→%.0f ', pnProfile.');
fprintf('dBc/Hz\n');
fprintf('PLL 环路滤波器 : %s, 低通 %d 阶, 带宽 %.0f kHz → 输出 PN 在此偏移处陡峭截断, 所有目标肩部截断偏移相同\n', ...
        onoffStr(enableLoopFilter), loopOrder, loopBW_Hz/1e3);
fprintf('PLL 扫频跟踪   : %s, 系数 %.1e s → σ_δf = %.2f kHz (= 系数×B/Ts), 残余随机相位 σ_φ(300m) = %.3f rad\n', ...
        onoffStr(enablePllTrack), pnTrackCoef, pnTrackCoef*FreqSlope/1e3, ...
        2*pi*pnTrackCoef*FreqSlope*2*300/c);
if pnTrackCoef > 0 && pnTrackCoef*FreqSlope > 1e6   % σδf > 1 MHz: 远超环路 ps~ns 抖动量级
    fprintf(['  ★警告: σ_δf = %.2g Hz 超出物理合理范围(环路时延抖动 ps~ns → σδf 应在 kHz 级),\n' ...
             '   仅可用于机理演示 —— 正常标定请把系数压回 1e-10~1e-9 量级\n'], pnTrackCoef*FreqSlope);
end
plkDb = 10*log10(Pt*sum(10.^(isoMat/10),1)*1e3);   % 各 RX 总泄漏功率 (dBm)
fprintf('TX-RX 泄漏     : %s, 各RX总泄漏 %.1f ~ %.1f dBm, 近距多径族 %s（全局修正 %+.0f dB）: ', ...
        onoffStr(enableLeakage), min(plkDb), max(plkDb), onoffStr(enableLeakNear), leakRes_dB);
for pp = 1:size(leakPaths, 1)
    fprintf('fb%.1fk(%+.0fdB) ', leakPaths(pp,1)/1e3, leakRes_dB + leakPaths(pp,2));
end
fprintf('非直流→无镜像; 恰直流 %s(%+.0fdB→bin0+镜像); 闪烁噪声抬高近距底噪\n', ...
        onoffStr(enableLeakDC), leakRes_dB + leakDC_dB);
ttDb = ttIsoMat + ttCorr_dB;  ttDb = ttDb(~isnan(ttDb));
fprintf('TX-TX 发发泄漏 : %s, 有效耦合 %.0f ~ %.0f dB (隔离度+修正%d dB), 固定随机相位, 串入分量携带原 TX 的 DDMA 码\n', ...
        onoffStr(enableTxCoupling), min(ttDb), max(ttDb), ttCorr_dB);
fprintf('PA 三阶交调    : %s, 双音 IM3 = %.0f dBc → 空子带组合码(2k_n−k_m)幅度 ≈ IM3×耦合度\n', ...
        onoffStr(enablePaIM3), paIM3_dBc);
rrDb = rrIsoMat + rrCorr_dB;  rrDb = rrDb(~isnan(rrDb));
fprintf('RX-RX 收收泄漏 : %s, 有效耦合 %.0f ~ %.0f dB (隔离度+修正%d dB), 固定随机相位, 线性叠加(无RX交调)\n', ...
        onoffStr(enableRxCoupling), min(rrDb), max(rrDb), rrCorr_dB);
fprintf('发射功率/增益  : Pt = %d dBm, Gt = %g dBi(含天线罩-%g), Gr = %g dBi(含天线罩-%g+链路-%g)（8 通道相同, 双程损耗共 -%g dB）\n', ...
        Pt_dBm, Gt_dBi, radomeLossDb, Gr_dBi, radomeLossDb, rfChainLossDb, 2*radomeLossDb + rfChainLossDb);
fprintf('ADC 量程       : %d bit 有符号整型（±%d counts = %.2f V, 1 LSB = %.1f μV）, 固定增益 G_rx = %d dB（无强目标时域约 ±200 counts）\n', ...
        adcBits, 2^(adcBits-1), adcFullScale_V, adcFullScale_V/2^(adcBits-1)*1e6, G_rx_dB);
fprintf('距离分辨率     : %.3f m,  最大不模糊距离 : %.1f m\n', DeltaR, Rmax);
fprintf('速度分辨率     : %.3f m/s, ±Vmax : %.2f m/s\n', Vres, Vmax);
fprintf('chirp 周期     : %.1f us (ADC 窗 %.1f us + 保护 %.1f us), 帧长 %.2f ms\n', ...
        Tc*1e6, Ts*1e6, T_idle*1e6, Tc*nChirps*1e3);
fprintf('---- bin 标定（RD 图横/纵轴即为 bin 序号）----\n');
fprintf('距离 bin: 第0个 = 0 m, 最后一个(第%d个) = %.2f m, 每bin宽 %.3f m\n', ...
        nSamplesValid-1, (nSamplesValid-1)*DeltaR, DeltaR);
fprintf('速度 bin: 第0个 = 0 m/s, 最后一个(第%d个) = %.2f m/s, 每bin宽 %.3f m/s（速度模糊, 等效±Vmax）\n', ...
        nChirps-1, Vres*(nChirps-1), Vres);
end                                             % if printCfg(配置块)
fprintf('================ 目标设置 ================\n');
for k = 1:Obj_amount
    % 远场近似的差拍频率与 LPF 增益（精确权重 W 按各通道真实几何距离计算）
    f_bk = 2*FreqSlope*Obj_R(k)/c + 2*Obj_v(k)/lambda;
    Hk   = 1/sqrt(1 + (f_bk/lpfFc)^(2*lpfOrder));
    % 两级模拟高通在差拍频点的解析增益(巴特沃斯高通 |H| = 1/√(1+(fc/f)^2n));
    % 近距目标 fb 落阻带会被显著压低(如角反 fb≈0.1MHz 处级联约 -37dB),
    % 逐目标存表供下游(RCS 标定)按已知滤形修正
    if enableHpf
        Hh1 = 1/sqrt(1 + (hpf1Fc/max(f_bk, 1))^(2*hpf1Order));
        Hh2 = 1/sqrt(1 + (hpf2Fc/max(f_bk, 1))^(2*hpf2Order));
    else
        Hh1 = 1;  Hh2 = 1;
    end
    % 两段式链路解析增益: 模拟LPF@真实fb → ADC 折叠(>fsAdc/2) → 逐级半带
    % (每级先滤波后 ÷2: 超过该级奈奎斯特的残余折到带内, 由后续级继续压) → 带边均衡
    fk = f_bk;  gk = Hk;
    if fk > fsAdc/2, fk = fsAdc - fk; end           % ADC 采样折叠
    for s = 1:nDStg
        gk = gk * hbMagOf(s, fk);
        if fk > hbR(s)/2, fk = hbR(s) - fk; end     % 该级 ÷2 后的折叠残余频率
    end
    filtGainDb(k) = 20*log10(gk * Hh1 * Hh2 * gC(fk));  % IF 链路总增益(含数字链+均衡) @ 差拍频点 (dB)
    filtFbHz(k)   = f_bk;
    % 接收功率（远场单站近似 Rt=Rr=R）: Pr = Pt·Gt·Gr·λ²·σ / ((4π)³·R⁴)
    Pr_k = Pt*Gt*Gr*lambda^2*RCS(k) / ((4*pi)^3 * Obj_R(k)^4);
    ampK = sqrt(RCS(k))/Obj_R(k)^2 * gk;            % 幅度 ∝ sqrt(Pr) ∝ sqrt(σ)/R² × 链路
    if k == 1, ampRef = ampK; end
    relDb = 20*log10(ampK/ampRef);                  % 相对目标1 的幅度 (dB)
    if f_bk < adcSampleRate/2
        fprintf(['目标%d: R = %7.2f m (bin %3d), v法 = %+6.2f m/s (bin %3d), v切 = %+5.1f m/s, ', ...
                 'σ = %6.1f m², fb = %5.2f MHz, 链路 %+6.2f HPF %+6.2f dB, Pr = %7.1f dBm, ', ...
                 '输入SNR %+5.1f dB, 幅度 %+6.1f dB\n'], ...
            k, Obj_R(k), round(f_bk/(adcSampleRate/nSamples)), ...
            Obj_v(k), mod(gate_D(k), nChirps), Obj_vTan(k), RCS(k), f_bk/1e6, ...
            20*log10(gk), 20*log10(Hh1*Hh2), 10*log10(Pr_k*1e3), 10*log10(Pr_k/Pn_W), relDb);
    else
        binAlias = round((adcSampleRate - f_bk)/(adcSampleRate/nSamples));
        fprintf(['目标%d: R = %7.2f m (超量程!), v = %+7.2f m/s, σ = %6.1f m², ', ...
                 'fb = %5.2f MHz > fs/2, Pr = %7.1f dBm, 输入SNR %+5.1f dB\n'], ...
            k, Obj_R(k), Obj_v(k), RCS(k), f_bk/1e6, 10*log10(Pr_k*1e3), 10*log10(Pr_k/Pn_W));
        fprintf(['       链路压制 %+6.2f dB 后残余混叠: 鬼影出现在视在距离 bin %d (%.1f m), ', ...
                 '幅度 ≈ %+6.1f dB (相对目标1), 且鬼影多普勒符号翻转(%+d→%+d)\n'], ...
            20*log10(gk), binAlias, binAlias*DeltaR, relDb, ...
            mod(gate_D(k), nChirps), mod(-gate_D(k), nChirps));
    end
end
fprintf('注: DDMA 使每个目标在多普勒维出现 8 个拷贝, 间隔 %d 个 bin。\n', nChirps/nGroups);

%% ============================ 阵列几何 ============================
%% ============================ 阵列几何 ============================
% TX_A / RX_A / antenna_t / antenna_r 由主函数注入(函数入参):
%   TX_A, RX_A    [1×8] TX/RX 通道幅度加权
%   antenna_t/r   [3×8] 天线相位中心坐标 [x;y;z] (m), y = 0 平面内

% 虚拟阵列坐标（TX 与 RX 相位中心矢量叠加, 顺序同 arrayManifold.m:
% TX 内循环 / RX 外循环 → TX1RX1..TX8RX1, TX1RX2..）
antenna_tr = arrayManifold(antenna_t, antenna_r);

%% ============================ 回波仿真 ============================
t = (0:nSamples-1)/adcSampleRate;             % 快时间轴 (s)

% DDMA 码槽位 ddmShift 由主函数注入(= find(encodeSeq==1) − 1, 各 TX 的槽位从 0 起)
DDMA_P = ddmShift * 360/nGroups;     % 各 TX 每 chirp 相位步进 (deg)

% ---- DDMA 移相器实际值模拟 ----
% 理想所需相位: φ_ideal(jj,n) = -DDMA_P(jj)·(n-1)  (deg)
% 实际移相器只能取 360/2^phaseBits 的整数倍 → 就近量化到最接近的可实现值；
% 每个 chirp 的相位设置不同，误差建模为逐 chirp 独立的正态随机误差。
% 注：当前 DDMA 步进为 22.5°=16×1.40625°，恰是 8-bit 的整数倍 LSB，
%     因此 8-bit 量化不引入误差（与理想一致）；改用更细的 DDMA 网格或
%     更低位数（≤3-bit）时量化效应才会显现。随机误差则会造成 DDMA
%     多普勒槽间泄漏（泄漏比 ≈ 20·log10(σ_rad/√nChirps)）。
phaseRes  = 360/2^phaseBits;                          % 移相器分辨率 (deg)
phiIdeal  = -DDMA_P(:)*(0:nChirps-1);                 % 理想相位 (nTX × nChirps, deg)
quantOn   = double(enablePhaseQuant);                 % 开关（false=0 不量化）
phiQuant  = (1-quantOn)*phiIdeal + quantOn*round(phiIdeal/phaseRes)*phaseRes;
DDMA_phaseQ = deg2rad(phiQuant);                      % 量化后的确定性相位 (rad)
% 注意：随机误差不在这里预生成 —— 在下方 chirp 循环内逐 chirp 独立抽取，
% 确保每个 TX × 每个 chirp 都是一次新的 N(0, σ²) 抽样。

% ---- TX-TX 发发隔离泄漏（TXm 信号串入 TXn 射频链路, 经 TXn 的 PA 放大后
%      从 TXn 天线再辐射）----
% TX n 天线实际辐射 = 自身 DDMA 调制 chirp + Σ_{m≠n} c_mn·(TX m 的 DDMA 调制 chirp)。
%   c_mn = 10^((ttIsoMat+ttCorr_dB)/20)·e^{jθmn}，θmn 为每对通道的固定相位
%   （黄金比例伪随机表 —— 只与通道对有关, 不随时间/频率/chirp/运行变化）。
% 固定随机相位防止各通道泄漏相位关系过好（避免相干相消/增强），
% 保证后续用 DDMA 通道泄露的能量也能测到目标。
% 串入分量从 TX n 天线辐射 → 几何路径按 TX n 位置（代码沿用天线 jj 的
% tau/几何, 仅把 DDMA 相位标签换成等效向量 DDMA_eff）。
% 对 DDMA 解码的影响：串入分量携带 TX m 的 DDMA 码 → 落在 TX m 的多普勒槽,
% 与直接路径同槽相干叠加（随机相位 → 不会相消）→ 污染虚拟阵列元素/
% 测角（幅度 ~c_mn 的阵列误差）。
txcOn = double(enableTxCoupling);                     % 开关（false=0 → cTx=I 无耦合）
ttLin = 10.^((ttIsoMat + ttCorr_dB)/20);              % 线性耦合幅度
ttLin(isnan(ttLin)) = 0;                              % 对角置 0（自发不耦合）
ttPh  = 2*pi*mod(((1:nTX).')*(1:nTX)*0.61803398875, 1);  % 固定伪随机相位表 θmn
cTx   = eye(nTX) + txcOn*(ttLin .* exp(1j*ttPh));     % cTx(m,n): m 串入 n 输出的系数

% ---- PA 三阶交调（发发耦合的必然结果：多码复合信号过非线性 PA）----
% 线性耦合只把能量搬进已占用的码槽；PA 的三阶非线性 y = x·(1+im3·|x|²)
% 作用于同载波多码复合包络 D_n = Σ c_mn·e^{jθm} 时，产生 θ_i−θ_j+θ_k 组合码
% —— 特别是 2k_n−k_m 分量（经典 2f1−f2 交调的相位码版本），幅度 ≈ im3×c_mn，
% 落入未占用的空子带 → 空子带相同距离 bin 出现微弱目标信号（与实测一致）。
% 双音等幅标定（精确）：两等幅单音过 y=x·(1+im3·(|x|²−1)) 时
% IM3/载波 = im3/(1+2·im3) → im3 = r/(1−2r), r = 10^(paIM3_dBc/20)。
% paIM3_dBc 为负值（IM3 低于载波 X dB, 典型 −20 ~ −5）。
% 空子带泄露电平 ≈ 耦合度 × im3（线性、不饱和）: paIM3_dBc ±1 dB → 泄露 ±1 dB。
% 注：无耦合时 |D|=1 → y=x 精确恒等 —— 交调空子带信号必须与发发耦合同时存在。
rIM3 = double(enablePaIM3) * 10^(paIM3_dBc/20);
rIM3 = min(rIM3, 0.45);                       % 防 r→0.5 发散（极限 PA）
im3Lin = rIM3 / max(1 - 2*rIM3, 0.1);         % 归一化三阶系数

% ---- RX-RX 收收隔离泄漏（RXm 信号串入 RXn 的前端, 线性叠加）----
% RX n 链路实际输入 = 自身信号 + Σ_{m≠n} r_mn·(RX m 的信号)。r_mn 为固定
% 复系数（隔离度+修正, 相位为固定伪随机表 —— 只与通道对有关, 不随时间/
% 频率/chirp/运行变化）。共锁相环/本振 → 串入信号经同一混频 → beat 域
% 直接线性叠加即可（在滤波链与实采样之前）。
% 不加 RX 三阶交调：串入与本链回波同为 µV 级小信号, LNA 交调产物比回波
% 低 ~75 dB、比热噪低 100+ dB, 不可见 —— 与 TX 情况不同（那里串入信号
% 要与主信号一起满功率过 PA, 交调不可避免）。
rrxOn = double(enableRxCoupling);                     % 开关（false=0 → rMat=I 无耦合）
rrLin = 10.^((rrIsoMat + rrCorr_dB)/20);              % 线性耦合幅度
rrLin(isnan(rrLin)) = 0;                              % 对角置 0
rrPh  = 2*pi*mod(((1:nRX).')*(1:nRX)*0.41421356237, 1);  % 固定伪随机相位表（与 TX 表不同）
rMat  = eye(nRX) + rrxOn*(rrLin .* exp(1j*rrPh));     % rMat(m,n): m 串入 n 输入的系数

% ---- 扫频非线性（所有 TX 共用一个锁相环 → 同频、同非线性；本振同源） ----
% 瞬时频率 f(t) = f_m + slope·t + ε(t)，四种可选形状（v = t/tLast ∈ [0,1]）：
%   'hump'     单凸:     ε = +A·sin(πv)
%   'dip'      单凹:     ε = −A·sin(πv)
%   'hump_dip' 先凸再凹: ε = +A·sin(2πv)   （前半凸、后半凹）
%   'dip_hump' 先凹再凸: ε = −A·sin(2πv)
% 统一写成 ε = s·A·sin(m·πv)，共同约束：ε(0)=ε(tLast)=0（首末采样点频率
% 不变、带宽保持 B）、缓变连续、ε' 两端为 0 与线性段平滑衔接。
% 步进频端点保证：锚定是【逐 chirp】的（每个 chirp 起止频点精确），因此
% 帧级起止频率恒为 startFreq 与 startFreq + (1+coupling_ratio)·B，
% 对任意 coupling_ratio 成立 —— 非线性只弯曲 chirp 内部轨迹，不动端点。
% 非线性相位 = 2π·∫ε·dt（解析闭式）：∫εdt = s·A·(tLast/(mπ))·(1−cos(mπv))，
% 同时加到发射相位与回波相位（回波按其发射时刻 delta_t = t−τ 求值）。
% dechirp 残余调频 δf ≈ τ·dε/dt = s·A·τ·(mπ/tLast)·cos(mπv)，对 R-FFT 的影响：
%   m=1: 差拍频率在 chirp 内单向扫一遍 → 展宽为连续团状、左右对称；
%        凸/凹只翻转扫频方向（时间反演），|FFT| 幅度谱形状相同；
%   m=2: 差拍频率往返扫两遍且在 chirp 中央出现驻留点（窗函数加权最大处）
%        → 能量堆积到一侧形成单边强瓣，hump_dip 与 dip_hump 互为镜像，
%        且残余调频幅度是 m=1 的 2 倍 → 展宽更宽。
switch nonlinType
    case 'hump',     nlSign = +1; nlHarm = 1;
    case 'dip',      nlSign = -1; nlHarm = 1;
    case 'hump_dip', nlSign = +1; nlHarm = 2;
    case 'dip_hump', nlSign = -1; nlHarm = 2;
    otherwise, error('nonlinType 必须为 hump / dip / hump_dip / dip_hump 之一');
end
tLast     = (nSamples-1)/adcSampleRate;       % 首末采样点时间间隔 (s)
nonlinAmp = (nonlinPct/100)*B;                % 最大瞬时频偏 = 非线性度×B (Hz)
nonlinPh  = @(tq) nlSign*nonlinAmp*tLast/(nlHarm*pi) ...
            .* (1 - cos(nlHarm*pi*tq/tLast)); % 2π·∫εdt/(2π)
nlOn      = double(enableNonlin);             % 开关（false=0 → 理想线性扫频）

% 各天线到目标的单程距离（决定发射/接收双程延迟）
Obj_t_R_all = zeros(nTX, Obj_amount);
for jj = 1:nTX
    Obj_t_R_all(jj, :) = sqrt(sum((Obj_xy - antenna_t(:, jj)).^2));
end
Obj_r_R_all = zeros(nRX, Obj_amount);
for ii = 1:nRX
    Obj_r_R_all(ii, :) = sqrt(sum((Obj_xy - antenna_r(:, ii)).^2));
end

% ---- 雷达方程幅度加权 + 模拟抗混叠低通（ADC 之前） ----
% MIMO 双站雷达方程: Pr = Pt·Gt·Gr·λ²·σ / ((4π)³·Rt²·Rr²)
%   （发射单程 Rt、接收单程 Rr 各平方一次；远场单站近似退化为 R⁴）
% IF 幅度 ∝ sqrt(Pr) ∝ sqrt(Pt·Gt·Gr·G_rx)·λ·sqrt(σ) / ((4π)^1.5·Rt·Rr)
% Pt/Gt/Gr/G_rx 为实际链路参数（8 通道相同）；TX_A/RX_A 为通道幅度增益。
% W 保持绝对物理刻度（不按目标归一化），与噪声的 ADC 域标准差同一刻度，
% 信噪比与 G_rx 取值无关。
%
% 模拟低通: 巴特沃斯 |H(f)| = 1/sqrt(1+(f/fc)^(2n))
%   每条 TX-目标-RX 通路是一个差拍单频分量，乘以 H(fb) 即为 LTI 系统
%   对单频的稳态响应 —— 数学上精确等效模拟滤波器（非数字硬截断）；
%   滤波后残余的带外分量在采样后自然混叠到带内形成鬼影。
%
% ---- 带内起伏（天线+射频链路的确定性幅度-频率响应，非随机） ----
% 形状 = 两个固定频率正弦的叠加（缓变、非单调、完全确定），峰峰值归一化
% 到 ripplePP_dB。评估域为【回波扫频射频】u = (f_RF − startFreq)/B：
% 步进频开启时，整帧实际覆盖的射频范围为 startFreq ~ startFreq+(1+cr)·B，
% u 全程 span [0, 1+coupling_ratio] —— 峰峰归一化也在此【全带】上进行
%（即 ripplePP_dB 是整个 (1+cr)·B 射频范围上的峰峰起伏，与硬件口径一致）。
% 链路增益随 f 缓变 → 等效于对发射/接收波形的调幅(AM)，dechirp 后差拍
% 单音被 AM 包络调制 → 距离维 FFT 旁瓣展宽（量级 ~ripplePP_dB/2）。
%（对比：若起伏在 IF 中频链 vs 差拍频率，则对单音只是常数增益、无展宽。）
% 8 个 TX / 8 个 RX 通过 txRipplePhase / rxRipplePhase 可独立设置形状相位，
% 当前全为 0（即 16 个通道响应完全一致，TX+RX 同相叠加峰峰最大 2×ripplePP_dB）。
rippleShape = @(u, ph) sin(2*pi*1.3*u + ph) + 0.4*sin(2*pi*2.7*u + 2.1 + ph);
txRipplePhase = zeros(1, nTX);         % TX 各通道起伏形状相位 (rad)，可独立设置
rxRipplePhase = zeros(1, nRX);         % RX 各通道起伏形状相位 (rad)，可独立设置
uSample = linspace(0, 1+coupling_ratio, 2001);  % 全带 [0, 1+cr] 归一化栅格
ppTx = zeros(1, nTX);                  % 各 TX 通道形状峰峰值（归一化前）
ppRx = zeros(1, nRX);                  % 各 RX 通道形状峰峰值（归一化前）
for jj = 1:nTX
    s = rippleShape(uSample, txRipplePhase(jj));
    ppTx(jj) = max(s) - min(s);
end
for ii = 1:nRX
    s = rippleShape(uSample, rxRipplePhase(ii));
    ppRx(ii) = max(s) - min(s);
end
rippleOn = double(enableRipple);       % 开关（false=0 → 起伏增益恒为 1）

% ---- 快时间过采样生成（消除"生成时即折叠"的建模缺陷）----
% 差拍信号以 osf·fs 生成 → IF 滤波链在过采样域对【真实频率】作用（fb>fs/2
% 的超量程目标真实存在、吃满阻带衰减），抽取回 fs 时才发生折叠 —— 与硬件
% "先模拟滤波、后采样折叠"的次序完全一致，无需任何混叠补偿。
osf = 2*chipOsr;    % 仿真过采样倍率 = 2×芯片 ADC 倍率(编程手段: 模仿模拟滤波对
                    % 真实频率作用需 2×ADC 率; chipOsr=1 时 =2 即旧行为)
fs2 = osf*adcSampleRate;                % 过采样率 (Hz)
i0  = round(adcStartDelay*fs2);         % ADC 首采样在过采样栅格上的偏移（0-based）
nS2 = i0 + osf*nSamples;                % 过采样快时间点数：覆盖 [0, Ts+adcStartDelay]
                                        %（扫频窗相应加长；ADC 只采 [adcStartDelay, +Ts)）
t2c = ((0:nS2-1)/fs2).';                % 过采样快时间列向量 (s)

W = zeros(nTX, nRX, Obj_amount);       % 每 TX-RX-目标 通路的最终权重
for jj = 1:nTX
    for ii = 1:nRX
        Rt = Obj_t_R_all(jj, :);       % 发射天线 → 目标 单程距离 (1×K)
        Rr = Obj_r_R_all(ii, :);       % 目标 → 接收天线 单程距离 (1×K)
        % 雷达方程幅度（含实际 Pt/Gt/Gr/G_rx 与通道增益，绝对刻度; 复数 σ 含初相）
        ampEq = sqrt(Pt*Gt*Gr*G_rx) * lambda .* sqrt(RCS_c) ./ ((4*pi)^1.5 * Rt .* Rr) ...
                * (TX_A(jj) * RX_A(ii));
        W(jj, ii, :) = ampEq;
    end
end

% IF 模拟滤波链（位于 dechirp 后、ADC 前）：两级高通 + 抗混叠低通，在【过
% 采样域 fs2】作用于全部谱成分 —— 目标载波、PN/带内起伏边带、泄漏单音、
% 闪烁噪声近端，以及 10~20 MHz 的带外真实频率分量（滤波后才抽取折叠）。
% 实现方式：巴特沃斯三节【因果 IIR】（butter 系数 + filter 带状态跨 chirp），
% 与真实模拟滤波器完全同构：因果（只依赖过去）、跨 chirp 连续。
% 教训（此前频域实现的三类伪影，实测均无）：
%   • ifft(fft(x).*H) 循环卷积：chirp 末端绕回头部 → 首末端爆出超满量程伪值；
%   • 零相位/线性相位响应：冲激非因果、需要"未来"输入 → 段尾反因果爬升；
%   • 倒谱最小相位：DC 深陷点(响应=0)使冲激带长数值尾 → 稳态抑制变差 7 倍。
%   • 两级高通: 压泄漏单音 + 压近距离目标幅度（防 ADC 饱和）；低通只负责
%     抗混叠与带外压制；PN 肩部"瞬间截断"不在这里（在 PLL 环路滤波器）。
if enableHpf
    [bH1, aH1] = butter(hpf1Order, hpf1Fc/(fs2/2), 'high');   % 高通1（混频输出级）
    [bH2, aH2] = butter(hpf2Order, hpf2Fc/(fs2/2), 'high');   % 高通2（ADC 前级）
else
    bH1 = 1;  aH1 = 1;   bH2 = 1;  aH2 = 1;                   % 关 → 全通
end
if enableLpf
    % SOS(二阶节)形式实现——16 阶直接型 b/a 在 fc≈0.49×fs2/2 处数值不稳定,
    % 必须逐节滤波(每节 2 个状态, 跨 chirp 连续), 响应与理论巴特沃斯一致
    [zL, pL, kL] = butter(lpfOrder, lpfFc/(fs2/2), 'low');
    sosLP = zp2sos(zL, pL, kL);                               % 抗混叠低通 [nSec×6]
else
    sosLP = [1 0 0 1 0 0];                                    % 关 → 全通单节
end

% ---- 相位噪声（共锁相环：8 个 TX 与本振为同一份相位噪声实现） ----
% 单边带相位噪声 L(f) 由 pnProfile 各点做 log-log 分段线性拟合：
%   f < 第一点：按 -20 dB/dec 外推（近载频 1/f² 区，经残余滤波后几乎无影响）；
%   各点之间：分段线性；f > 最后一点：按最后一段斜率继续外推。
% 时域合成（随机相位法，覆盖整帧时间轴含 chirp 间隔）：
%   phi = Re{IFFT{A(f)·e^{jθ}}}，A = sqrt(fs·M·L(f))，共轭对称保证实序列，
%   统计方差 = 2∫L(f)df（与 L 的标准定义一致）。
% dechirp 时本振与回波携带同一 φ(t)，共模对消，残余 = φ(t) − φ(t−τ)，
%   其 PSD 乘 |1−e^{−j2πfτ}|²：近端 PN 被强烈对消，残余随 τ（即距离）增大
%   → R-FFT 中目标谱两侧出现"肩部"，肩部电平随距离抬高。
pnOn  = double(enablePhaseNoise);             % 开关（false=0 → 理想无相位噪声）
M_pn  = round((nChirps-1)*Tc*fs2) + nS2 + 64;  % 整帧 PN 序列长度（过采样率）
fAx   = (0:floor(M_pn/2)).' / M_pn * fs2;      % 单边频率轴 (Hz, 至 fs2/2)
fAxC  = max(fAx, 1);                          % 防 log10(0)（DC 频偏极低可忽略）
Lfit  = interp1(log10(pnProfile(:,1)), pnProfile(:,2), log10(fAxC), 'linear');
loEnd = fAxC < pnProfile(1,1);                % 低于第一点的近端区
Lfit(loEnd) = pnProfile(1,2) - 20*(log10(pnProfile(1,1)) - log10(fAxC(loEnd)));
hiEnd = fAxC > pnProfile(end,1);              % 高于最后一点（过采样轴可达 fs2/2）
slpEnd = (pnProfile(end,2)-pnProfile(end-1,2)) / ...
         (log10(pnProfile(end,1))-log10(pnProfile(end-1,1)));   % 末段斜率 dB/dec
Lfit(hiEnd) = pnProfile(end,2) + slpEnd*(log10(fAxC(hiEnd)) - log10(pnProfile(end,1)));
Lphi  = 10.^(Lfit/10);                        % L(f) 线性值 (1/Hz)；积分 rms² = 2∫L df
% PLL 环路滤波器（低通）整形：输出 PN 的 PSD 在环路带宽处被陡峭截断。
% 环路内（< loopBW）按拟合点；环路外按 |H_loop|² 陡降 —— R-FFT 中肩部在
% 偏移 = loopBW 处瞬间截断，且所有目标的截断偏移相同（环路特性与目标无关）。
% enableLoopFilter=false 时不乘 |H_loop|²：输出 PN 按拟合点直接外推、不截断
%（环路整形只作用于 PN 合成，故仅在 enablePhaseNoise=true 时有实际效果）。
loopOn = double(enableLoopFilter);                              % 开关（false=0 → 无环路整形）
Lphi = Lphi ./ (1 + loopOn*(fAxC/loopBW_Hz).^(2*loopOrder));    % 乘 |H_loop(f)|²
Lphi(1) = 0;                                  % DC 置零（绝对相位无意义）
A_pn  = sqrt(fs2*M_pn*Lphi);                  % 频域幅度（随机相位法方差 = 2∫L df）
Xf_pn = A_pn .* exp(1j*2*pi*rand(numel(A_pn),1));
Xf_pn(1) = 0; Xf_pn(end) = 0;                 % DC 与 Nyquist 置零
Xfull = zeros(M_pn, 1);
Xfull(1:floor(M_pn/2)+1) = Xf_pn;
Xfull(M_pn:-1:M_pn-floor(M_pn/2)+1) = conj(Xf_pn(2:end));     % 共轭对称 → 实序列
pnFull = real(ifft(Xfull)) * pnOn;            % 全帧相位噪声实现 (rad)

% 本振(=发射)与各目标回波在帧网格上的 PN 取值表
% （各天线对的 τ 差为 ns 级 ≪ PN 相干时间，按目标平均双程时延取值即可）
tauPn  = (mean(Obj_t_R_all,1) + mean(Obj_r_R_all,1)).' / c;   % 各目标平均双程时延
pnLO   = zeros(nChirps, nS2);                 % 本振/发射时刻 PN [chirp, 过采样快时间]
pnEcho = zeros(nChirps, nS2, Obj_amount);     % 回波发射时刻 PN [chirp, 快时间, 目标]
gx     = (1:M_pn).';
for n = 1:nChirps
    base = round((n-1)*Tc*fs2);               % 本 chirp 帧内起始样点（过采样）
    pnLO(n,:) = pnFull(base + (1:nS2)).';
    for k = 1:Obj_amount
        idxE = max(1, base + (1:nS2) - tauPn(k)*fs2);
        pnEcho(n,:,k) = interp1(gx, pnFull, idxE, 'linear');
    end
end
pnRms = sqrt(mean(pnFull.^2));
if printCfg, fprintf('相位噪声实现   : 积分 rms = %.2f rad (%.1f°)\n', pnRms, rad2deg(pnRms)); end

% 泄漏载波的残余 PN：泄漏路径经等效去相关延迟 leakPnTau 后与本振不完全共模
% → 泄漏线带 PN 裙边，PSD = L(f)·|1−e^{−j2πf·leakPnTau}|²（L 已含环路整形）。
% 泄漏线比目标高 ~70 dB，裙边浮于热噪之上，并在环路带宽偏移处随环路滤波器
% 陡峭截断 —— 与各目标肩部的截断偏移一致。
pnLeak = zeros(nChirps, nS2);               % 泄漏时刻的 PN 取值表（过采样）
for n = 1:nChirps
    base = round((n-1)*Tc*fs2);
    idxL = max(1, base + (1:nS2) - leakPnTau*fs2);
    pnLeak(n,:) = interp1(gx, pnFull, idxL, 'linear');
end

% ---- TX-RX 隔离泄漏（RX 直接收到 TX 能量，未经目标反射） ----
% 泄漏载波 = 幅度 sqrt(Pt·10^(iso/10)·G_rx) · DDMA 码 · 每路径固定随机相位
% leakPh（防止测角模块把泄漏误判为点目标）· 各路径延迟 leakTauP（多径表 leakPaths）。
% 关键物理：泄漏路径延迟非零 → 差拍 = slope·τ 落在各自距离 bin
%（视在距离 = 单程路径的一半），【不是直流】：
%   • 恰好直流（τ=0，仅片内直通）的成分被两级高通压死 → bin0 无峰（实测如此）
%   • 近距 bin（1~3）的码梳经高通部分压制后残留 → 可见的 8 码梳（实测如此）
%   • 非零差拍 → 该距离行慢时间序列为复数 → 无实采样共轭镜像（实测如此）
%   • 残余 PN 裙边（pnLO−pnLeak，同一份 PN 乘所有码 → 每根码线裙边一致）
% 泄漏强于目标 ~70 dB：enableHpf=false 时该码梳以全幅度出现。
%
% 闪烁噪声：泄漏能量抬升 RX 1/f^α 基带噪声，PSD = P_leak·K/f^α（正比该
% RX 总泄漏功率）。α=0.5 匹配实测形状：RX4 近距底噪 +10 dB、随距离缓降、
% 接近最大测距才与热噪合并；其余 RX 闪烁≪白噪声（隔离度再差 20 dB 的
% 通道底噪也看不出差别，与实测一致）。
lkOn     = double(enableLeakage);                 % 总开关（false=0 → 无泄漏、无闪烁）
lnOn     = lkOn * double(enableLeakNear);         % 近距多径族开关（总开关与统一开关同时开才有效）
% 泄漏多径族预计算：各路径等效延迟与电平（多径表 leakPaths，fb 单位 Hz —— 照 96 修法）
nLeakP   = size(leakPaths, 1);                    % 路径条数
leakTauP = leakPaths(:,1) / FreqSlope;            % 各路径延迟 (s) = fb/斜率（物理固定, 与采样配置无关）
A_leak   = sqrt(Pt * G_rx * 10.^(isoMat/10));     % 8×8 泄漏基准幅度（ADC 域, 不含电平修正）
leakAmpP = 10.^((leakRes_dB + leakPaths(:,2).')./20);  % 各路径电平系数（全局 leakRes + 逐路径修正）
leakPh   = 2*pi*rand(nTX, nRX, nLeakP);           % 每 TX-RX-路径 固定随机相位 (rad)
lnDC     = lkOn * double(enableLeakDC);           % 恰直流路径开关（独立于近距族）
ampDC    = 10^((leakRes_dB + leakDC_dB)/20);      % 直流路径电平系数
leakPhDC = 2*pi*rand(nTX, nRX);                   % 直流路径固定随机相位 (rad)
P_leakRX = Pt * sum(10.^(isoMat/10), 1).';        % 每 RX 总泄漏功率 (W)
% 每 RX 独立闪烁噪声（ADC 域，随机相位法合成，覆盖整帧时间轴）
fFlk  = fAx;  fFlk(1) = 1;                        % 防除零
flickADC = zeros(nRX, M_pn);                      % [RX, 帧样点]
for ii = 1:nRX
    % 低于第一个距离 bin（Δf=fs/N）的深低频分量压平为 S(Δf)：这部分物理上
    % 落在直流附近、由链路高通（enableHpf）处理；不压平且高通关闭时
    % 会全部涌入 bin0-2。
    fFlkC = max(fFlk, adcSampleRate/nSamples);     % 压平后的频率轴
    Sfl = G_rx * P_leakRX(ii) * flickerK ./ fFlkC.^flickerAlpha;  % PSD (V²/Hz)
    Sfl(1) = 0;
    Ak  = sqrt(fs2*M_pn*Sfl);
    Xk  = Ak .* exp(1j*2*pi*rand(numel(Ak),1));  Xk(1) = 0;  Xk(end) = 0;
    Xg  = zeros(M_pn, 1);
    Xg(1:floor(M_pn/2)+1) = Xk;
    Xg(M_pn:-1:M_pn-floor(M_pn/2)+1) = conj(Xk(2:end));
    flickADC(ii,:) = real(ifft(Xg)).' * lkOn;
end

% 数据缓存
s_r_normal  = zeros(nRX, nSamples);           % 原始回波（复数，未 dechirp）
zF1 = zeros(numel(aH1)-1, nRX);               % 三节 IIR 滤波器状态（列×RX, 跨 chirp 连续）
zF2 = zeros(numel(aH2)-1, nRX);
zLP = zeros(size(sosLP, 1), 2, nRX);   % 低通 SOS 各节滤波器状态 [节×2×RX], 跨 chirp 连续
sg_r_normal = zeros(nRX, chipOsr*nSamples);  % 实采样 ADC 输出（实数, ADC 域 = chipOsr×fs）
adcRawData  = zeros(chipOsr*nSamples, nChirps, nRX);  % 一帧 ADC 域数据 [快时间, 慢时间, RX]

t_fast = t(:);                                % 列向量快时间

% ---- 馈线长度色散预计算（TX/RX 16 条馈线不等长 → 通道相位随频率失配） ----
% 每条馈线相位差 φ(f) = 2π·(f − f_cal)·ΔL/c, f_cal = 全带（B + 步进频）中心
% 频率 f_center → 中心频率处 16 通道相位相等（假设已校准）, 偏离中心按线长
% 差线性计算。信号路径：发射 chirp 经 TX jj 馈线辐射（瞬时频率 f_m+slope·Δt）、
% 回波经 RX ii 馈线进入接收机（频率在传播中不变, 同一瞬时频率）→ 双路过线
% 相位合并为系数 (ΔL_tx+ΔL_rx)·2π/c；dechirp 本振不经馈线 → 不受影响。
% 注：TX-RX 泄漏路径未加馈线相位（二级小效应, 主效应为目标回路的通道失配）。
fdOn  = double(enableFeedDisp);
assert(numel(feedLenDev) == nTX + nRX, ...
    'feedLenDev 条目数(%d) 必须 = TX数+RX数(%d+%d) —— 改天线数量时同步改此旋钮表', ...
    numel(feedLenDev), nTX, nRX);
dLtx  = fdOn * feedLenScale * feedLenDev(1:nTX);           % TX 馈线长度差 (m)
dLrx  = fdOn * feedLenScale * feedLenDev(nTX+1:end);       % RX 馈线长度差 (m)
linePhCoef = 2*pi/c * (dLtx(:) + dLrx(:).');       % nTX×nRX [TX,RX] 相位系数 (rad/Hz)

% ---- 回波张量算力(GPU/CPU)与小常量预搬(张量第一遍用, 详见 chirp 循环内注释) ----
useGpu = false;
if strcmp(gpuMode, 'on') || strcmp(gpuMode, 'auto')
    try
        if gpuDeviceCount > 0
            % 新卡(算力高于自带 CUDA 库, 如 RTX50 系)必须先开前向兼容再碰设备 ——
            % 该开关是【会话级】的, 且对老卡无害; 顺序反了会触发失败的驱动调用
            parallel.gpu.enableCUDAForwardCompatibility(true);
            gpuDevice;
            useGpu = true;
        end
    catch
        useGpu = false;                                       % GPU 不可用 → 纯CPU张量
    end
end
if useGpu
    if printCfg, fprintf('回波张量算力   : GPU (%s)\n', gpuDevice(1).Name); end
    W4         = reshape(gpuArray(W), nTX, nRX, 1, []);
    linePhCoefG = gpuArray(linePhCoef);
    txPh4 = reshape(gpuArray(txRipplePhase), nTX, 1, 1, 1);
    ppTx4 = reshape(gpuArray(ppTx),          nTX, 1, 1, 1);
    rxPh4 = reshape(gpuArray(rxRipplePhase), 1, nRX, 1, 1);
    ppRx4 = reshape(gpuArray(ppRx),          1, nRX, 1, 1);
else
    if printCfg, fprintf('回波张量算力   : CPU (张量化)\n'); end
    W4          = reshape(W, nTX, nRX, 1, []);
    linePhCoefG = linePhCoef;
    txPh4 = reshape(txRipplePhase, nTX, 1, 1, 1);
    ppTx4 = reshape(ppTx,          nTX, 1, 1, 1);
    rxPh4 = reshape(rxRipplePhase, 1, nRX, 1, 1);
    ppRx4 = reshape(ppRx,          1, nRX, 1, 1);
end

for n = 1:nChirps
    t_m = (n-1)*Tc;                           % 慢时间：本 chirp 起始时刻
    f_m = startFreq + (n-1)*f_step;           % 本 chirp 起始频率（理想值）
    % 步进频跳频误差：起始频率 = 理想值 + 逐 chirp 独立正态随机偏差（不累积,
    % 以理想值为基准）。δf 同时进入发射相位与本振 → 差拍频率不变（slope·τ）,
    % 残余 = 逐 chirp 随机相位 2π·δf·τ（随时延/距离增大）; f_m 同时用于泄漏
    % 相位项（泄漏路径 τ 为 ns 级 → 影响可忽略, 与物理一致）
    f_m = f_m + double(enableHopErr) * hopErr_Hz * randn;
    % PLL 扫频跟踪误差（开关区第 8b 块）：σ_δf = pnTrackCoef×B/Ts, 逐 chirp 独立。
    % 与跳频误差同注入点 → 同机理: 残余逐 chirp 随机相位 2π·δf·τ, 随距离/带宽
    % 增长 → 目标两侧肩膀随 B 抬高（pnProfile 不含 B 依赖, 补的正是这块）。
    % 系数 0 时不抽取 randn（随机流不动 → 与旧结果逐位一致）
    if enablePllTrack && pnTrackCoef > 0
        f_m = f_m + pnTrackCoef * FreqSlope * randn;
    end
    slow_time    = t_m;
    t_fast_slow  = t2c + slow_time;           % 绝对时间（过采样栅格）

    % 发射信号相位：2π[ f_m*(t_m + t) + slope/2*t² ] + 扫频非线性 + 相位噪声
    % （相位噪声为共锁相环实现在本 chirp 各采样时刻的取值；过采样栅格）
    s_tx = exp(1j*( 2*pi*( f_m*(t_m + t2c) + 0.5*FreqSlope*t2c.^2 ...
             + nlOn*nonlinPh(t2c) ) + pnOn*pnLO(n,:).' ));

    % DDMA 相位调制：量化后的确定相位 + 移相器随机误差
    % （误差在循环内抽取 randn(nTX,1)：8 个 TX、本 chirp 各自独立的一次抽样）
    DDMA_vec = exp(1j*( DDMA_phaseQ(:, n) ...
               + double(enablePhaseErr)*deg2rad(phaseErrSigma_deg)*randn(nTX,1) ));

    % TX-TX 耦合后的实际辐射 DDMA 向量（TX n 辐射 = 自身码 + 串入的其它 TX 码）
    DDMA_eff = cTx.' * DDMA_vec;
    % PA 三阶交调：y = x·(1 + im3·(|x|²−1))——增益在标称工作点(单码 |x|=1)
    % 归一（|x|=1 时 y=x 精确成立），不整体抬高回波绝对电平；|x|² 随 chirp
    % 旋转 → 展开产生 2k_n−k_m 组合码（含空子带分量， 幅度 ≈ im3×耦合度，
    % 随 paIM3_dBc 线性可控、无饱和）
    DDMA_eff = DDMA_eff .* (1 + im3Lin*(abs(DDMA_eff).^2 - 1));

    % ---- 目标运动几何(本 chirp): 法向+切向都进几何 ----
    % 本 chirp 起点位置 p(t_m) = p0 + (vRad·û + vTan·t̂)·t_m —— 切向走动改变视角,
    % 逐 chirp 重算各天线到目标的单程距离(慢时间多普勒与阵列相位缓变都在这一项里)
    pObj3  = Obj_xy + Obj_uHat.*(Obj_v*t_m) + Obj_tHat.*(Obj_vTan*t_m);   % [3×K]
    Rt_cur = sqrt(sum((reshape(pObj3, 3, 1, []) - reshape(antenna_t, 3, [], 1)).^2, 1));  % [1×nTX×K]
    Rt_cur = reshape(Rt_cur, nTX, []);              % 发射天线→目标 单程距离 [nTX×K]
    Rr_cur = reshape(sqrt(sum((reshape(pObj3, 3, 1, []) - reshape(antenna_r, 3, [], 1)).^2, 1)), nRX, []);
    % chirp 内距离增量: dRw = |p(t)| − |p(t_m)|, t = 本 chirp 内绝对时刻(快时间栅格)。
    % 【注意】增量只算 chirp 内一段 —— 慢时间的径向累计已含在 Rt_cur/Rr_cur 里,
    % 若从帧首累计会把径向运动算两遍(多普勒翻倍)。|p(t)| 精确式含切向二阶项;
    % 纯法向时 dRw = v·(t−t_m), 与旧线性模型一致; 切向不产生多普勒(只进二阶项)
    Rm  = sqrt(Obj_R.^2 + (2*Obj_R.*Obj_v + Obj_v2.*t_m).*t_m);                    % |p(t_m)| [1×K]
    dRw = sqrt(Obj_R.^2 + (2*Obj_R.*Obj_v + Obj_v2.*t_fast_slow).*t_fast_slow) - Rm;  % [nS2×K]

    % ---- 第一遍(目标回波)：64 个 (TX,RX) 对按张量【一次算完】 ----
    % 全部元素级运算广播到 [nTX×nRX×nS2×K] —— 免去 64×512 次小阵列函数调用,
    % 且张量可整体搬 GPU(gpuMode 旋钮)。数值与原逐对循环逐步等价(仅浮点求和
    % 次序不同, 差异在 1e-15 相对量级)。物理注释见各行的原逐对版本(已吸收)。
    sgRaw = zeros(nRX, nS2);                          % [RX, 过采样快时间] 复数前端信号
    leakSigAll = zeros(nRX, nS2);                     % 泄漏单音逐 RX 缓存（旁路残余用）
    if useGpu
        Rt_cur = gpuArray(Rt_cur);  Rr_cur = gpuArray(Rr_cur);
        dRw = gpuArray(dRw);  t4base = gpuArray(t2c);
        pnE4 = gpuArray(reshape(pnEcho(n,:,:), 1, 1, nS2, []));
        DDMA4 = reshape(gpuArray(DDMA_eff), nTX, 1, 1, 1);
    else
        t4base = t2c;
        pnE4 = reshape(pnEcho(n,:,:), 1, 1, nS2, []);
        DDMA4 = reshape(DDMA_eff, nTX, 1, 1, 1);
    end
    tau4  = (reshape(Rt_cur, nTX, 1, 1, []) + reshape(Rr_cur, 1, nRX, 1, []) ...
             + 2*reshape(dRw, 1, 1, nS2, [])) / c;    % 双程时延(本chirp几何+chirp内增量)
    dt4   = reshape(t4base, 1, 1, nS2, 1) - tau4;     % 回波对应的发射时刻（过采样栅格）
    mask4 = (dt4 >= 0) & (dt4 < Ts + adcStartDelay);  % 本 chirp 扫频窗内才有差拍输出
    ph4   = 2*pi*( f_m*(t_m + reshape(t4base,1,1,nS2,1) - tau4) + 0.5*FreqSlope*dt4.^2 ...
             + nlOn*nonlinPh(dt4) ) ...
             + pnOn*pnE4;                                            % 回波时刻相位噪声
    ph4   = ph4 + reshape(linePhCoefG, nTX, nRX, 1, 1) ...
              .* (f_m - f_center + FreqSlope*dt4);                  % 馈线色散相位
    u4    = ((n-1)*f_step + FreqSlope*dt4 - min(0, f_step*(nChirps-1))) / B;  % 回波射频归一化位置; 平移到固定 span [0,1+cr](下跳频时 (n-1)f_step<0)
    g4    = (1 - rippleOn) + rippleOn * 10.^(ripplePP_dB/20 * ( ...
              rippleShape(u4, txPh4) ./ ppTx4 ...                    % TX 带内起伏(AM)
            + rippleShape(u4, rxPh4) ./ ppRx4 ));                    % RX 带内起伏(AM)
    S4    = exp(1i*ph4) .* mask4 .* g4 .* DDMA4;
    echoAll = reshape(sum(sum(S4 .* W4, 1), 4), nRX, nS2);   % 幅度加权(W)后对 TX/K 求和
    if useGpu, echoAll = gather(echoAll); end

    % ---- 第一遍(泄漏+混频)：逐 RX（量小, 顺序处理） ----
    for ii = 1:nRX
        s_rx_total = echoAll(ii, :).';                     % 本 RX 各 TX 各目标回波之和

        s_r_normal(ii, :) = s_rx_total(i0+1:osf:end).';  % 原始回波（复数, ADC 窗, 备用）

        % dechirp 混频：复数差拍信号（混频和频项被理想低通滤除，只余差频）
        sg_r = s_tx .* conj(s_rx_total);

        % TX-RX 泄漏多径族（leakPaths 表逐条累加）：携带 DDMA 与每条路径
        % 固定随机相位 + 残余 PN 裙边（pnLO−pnLeak 共用等效去相关延迟，
        % 随距离/偏移缓降，在 IF-LPF 截止处被瞬间截断）。
        % 相位按物理准确建模 = 2π·f_m·τ（逐 chirp，含步进频跳变 → 码梳随 cr
        % 产生与目标相同的距离-多普勒耦合偏移）+ 2π·slope·τ·t'（chirp 内）。
        % 主份送滤波链（受高通整形）；缓存一份给"旁路残余"（第二遍）
        leakVec = zeros(nS2, 1);
        for pp = 1:nLeakP
            beatP  = ((A_leak(:,ii)*leakAmpP(pp)) .* exp(-1j*leakPh(:,ii,pp))).' * conj(DDMA_eff);
            leakVec = leakVec + lnOn * beatP .* exp(1j*( 2*pi*f_m*leakTauP(pp) ...
                   + 2*pi*FreqSlope*leakTauP(pp)*t2c ...
                   + pnOn*(pnLO(n,:).' - pnLeak(n,:).') ));
        end
        % 恰直流泄漏路径（τ=0，独立开关 enableLeakDC）：差拍=0 → bin0；
        % τ=0 → PN 完全共模对消（φ(t)−φ(t)=0）→ 无 PN 裙边，相位只剩 DDMA 码
        % 与固定随机相位 → bin0 行慢时间序列为实数 → 实采样多普勒共轭镜像
        beatDC = ((A_leak(:,ii)*ampDC) .* exp(-1j*leakPhDC(:,ii))).' * conj(DDMA_eff);
        leakVec = leakVec + lnDC * beatDC;
        sg_r = sg_r + leakVec;
        leakSigAll(ii,:) = leakVec.';

        sgRaw(ii, :) = sg_r.';
    end

    % ---- RX-RX 收收耦合：RX n 链路输入 = Σ_m rMat(m,n)·(RX m 的前端信号)
    %（固定相位线性叠加, 位于各链噪声/滤波/采样之前; 共本振 → beat 域
    %   直接叠加即可; 串入内容包含各 RX 的目标回波与 TX-RX 泄漏载波） ----
    sgRaw = rMat.' * sgRaw;

    % ---- 第二遍：各 RX 链路（本链噪声 + IF 滤波链 + 抽取 + 实采样） ----
    for ii = 1:nRX
        sg_r = sgRaw(ii,:).';

        % 接收机白噪声与泄漏闪烁噪声：在滤波链【前】加入（物理上产生于
        % 混频器/链路前级），随后与信号一起经过完整 IF 链（两级高通+低通）：
        %   • 近距 bin 的底噪同样被高通压低几十 dB，随频率按 ~f^(2·Σn) 幂律
        %     回升（对数频率下呈抛物线），约 bin 50 处恢复 —— 与实测一致；
        %   • 热噪与闪烁受同一链路响应，两者比值不变 → RX4 抬高量保持。
        %（过采样白噪声 σ ×√osf, 经链路+抽取后带内方差与 fs 标定一致）
        sg_r = sg_r + double(enableNoise)*noiseSigmaW*sqrt(osf)*randn(nS2,1);
        idxFlk = round((n-1)*Tc*fs2) + (1:nS2);
        sg_r = sg_r + lkOn * flickADC(ii, idxFlk).';

        % IF 模拟滤波链（三节因果 IIR, 状态跨 chirp 连续 —— 真实模拟滤波器行为:
        % 因果、连续、无循环卷积环绕、无未来输入依赖, 在过采样域以真实频率作用）
        [y1, zF1(:,ii)] = filter(bH1, aH1, sg_r, zF1(:,ii));   % 高通 1
        [y2, zF2(:,ii)] = filter(bH2, aH2, y1,   zF2(:,ii));   % 高通 2
        sg_r = y2;                                          % 抗混叠低通(SOS 逐节, 状态跨 chirp)
        for ss = 1:size(sosLP, 1)
            [sg_r, zLP(ss, :, ii)] = filter(sosLP(ss,1:3), sosLP(ss,4:6), sg_r, zLP(ss, :, ii));
        end

        % 泄漏"旁路残余"：部分泄漏绕过高通耦合到后级（等效耦合 hpfFloor_dB），
        % 不经历高通整形、直接叠加在滤波输出上 → 近距码梳的来源（高度可标定）
        sg_r = sg_r + 10^(hpfFloor_dB/20) * leakSigAll(ii,:).';

        % ADC 采样（第一级抽取: 仿真域 fs2 → 芯片 ADC 域 fsAdc, ÷2）: 折叠发生在
        % 模拟滤波【之后】—— 超量程目标先吃满真实频率处的模拟衰减再混叠(与硬件一致)。
        % 后续 ÷2^nDStg 由数字半带链完成(循环外整帧做), 合计仍落在原 osf 栅格
        sg_r = sg_r(i0+1:2:end);              % ADC 采样窗：跳过 chirp 首部建立瞬态

        % ★ 实采样 ADC(取实部; 芯片噪声/量程检查/量化在 ADC 域整帧做, 见循环外) ★
        % 单路实数采样等效于取复差拍信号的实部：
        %   real{ s_tx * conj(s_rx) } = cos(φ_tx − φ_rx) —— 实混频器差频输出。
        % 实信号频谱共轭对称，每个目标会在 ±fb 出现一对峰，
        % 后续距离维只取前 nSamples/2+1 个 bin（0 ~ fs/2）即可。
        sg_r_normal(ii, :) = real(sg_r).';
    end

    adcRawData(:, n, :) = sg_r_normal.';
end

% 注意：无数字去直流。泄漏码梳的压制由 IF 链路里的两级模拟高通（enableHpf）
% 完成——恰好直流（片内直通）被压死 → bin0 无峰；各泄漏路径差拍（leakPaths
% 多径族, 近距 bin）部分压制 → 近距 bin 出码梳（复数行、无实采样镜像）；近距目标
% 幅度同被压低（防 ADC 饱和）。enableHpf=false 时码梳以全幅度出现。各路径
% 固定随机相位使泄漏跨通道不相干，避免测角模块误判为近距点目标。

%% ==================== ADC 时域量程检查（数据转 int16 counts 单位） ====================
% 硬件 ADC 输出 2 字节有符号整型：数值即 counts，范围 ±2^(adcBits-1)（int16:
% −32768 ~ +32767）。仿真按约定不做整数量化（counts 取连续值），ADC 数据统一
% 以 counts 为单位输出：counts = 电压 / (满量程/32768)。
% 时域样点必须落在 ±32768 内（超出 = 硬件削顶，此处只告警、不削波）；
% 2D-FFT 后数值超过 ±32768 属相干积累增益，正常现象、不是超限。
% 增益为固定值 G_rx_dB（无自动增益）：标定 = 无强目标时（弱目标+泄漏+底噪）
% 时域约 ±200 counts；若本场景超限，请降低 G_rx_dB 或开 enableHpf。
vPerLSB   = adcFullScale_V / 2^(adcBits-1);                % 每 LSB 对应电压 (V)
adcRawData = adcRawData / vPerLSB;                         % ★ V → int16 counts（不取整）
peakRXc  = squeeze(max(max(abs(adcRawData), [], 1), [], 2)).';  % 各 RX 全帧峰值 (counts)
adcFSmax = 2^(adcBits-1);                                  % 满量程 counts (int16: 32768)
fprintf('---- ADC 时域量程（%d bit, ±%d counts = %.2f V, 1 LSB = %.1f μV, 固定增益 %d dB）----\n', ...
        adcBits, adcFSmax, adcFullScale_V, vPerLSB*1e6, G_rx_dB);
for rx = 1:nRX
    util = peakRXc(rx)/adcFSmax*100;
    tag  = '';
    if peakRXc(rx) > adcFSmax, tag = ' ← 超限(硬件削顶)!'; end
    fprintf('  RX%d 峰值 ±%6.0f counts (%5.1f%%FS)%s\n', rx, peakRXc(rx), util, tag);
end
if max(peakRXc) > adcFSmax
    warning('ADC 超限: 峰值 %.0f counts > ±%d —— 请降低 G_rx_dB 或开启 enableHpf 压近距/泄漏', ...
            max(peakRXc), adcFSmax);
end

% ---- ADC 芯片自身噪声(数据手册 SNR 口径, 输入参考白噪) ----
% 手册 SNR = 满幅正弦 RMS / 全噪声 RMS(物理电压口径)。噪声先在电压域算好:
%   σ_V = (满量程电压/√2) / 10^(SNR/20) —— 只取决于手册 SNR 与 adcFullScale_V,
%   与 adcBits 无关(位数是数字侧网格, 不改变芯片的物理噪声)。
% 再折算到当前 counts 口径(÷vPerLSB)加到数据上。改 adcBits 时同一份物理噪声
% 自动换算, 不产生额外影响。加在量化之前、与信号无关的高斯白噪 —— 与"减小
% adcBits"有本质区别: 细网格仍在(小信号仍可分辨), 只是多了一层平坦噪声。
if enableAdcChipNoise
    if exist('adcChipSigmaV_cal', 'var') && ~isempty(adcChipSigmaV_cal)
        adcChipSigmaV = adcChipSigmaV_cal;               % 实测标定电压(优先, 与满量程无关)
        snrTxt = sprintf('实测标定 %.1f μV rms(等效 SNR %.1f dB)', ...
            adcChipSigmaV*1e6, 20*log10((adcFullScale_V/sqrt(2))/adcChipSigmaV));
    else
        adcChipSigmaV = (adcFullScale_V/sqrt(2)) / 10^(adcChipSnrDb/20);  % 手册 SNR 口径
        snrTxt = sprintf('手册 SNR %.1f dB', adcChipSnrDb);
    end
    adcChipSigma  = adcChipSigmaV / vPerLSB;                          % 折算当前 counts 口径
    adcRawData = adcRawData + adcChipSigma*randn(size(adcRawData));
    if printCfg, fprintf('ADC 芯片噪声   : 开 —— %s → σ = %.2f counts(当前位数口径)\n', ...
            snrTxt, adcChipSigma); end
end

% ---- 数字增益(纯数值缩放; 芯片噪声之后、溢出检查/量化之前) ----
% 与真实芯片 DFE/位宽打包的数字增益级对应: 信号与全部噪声同乘一个系数,
% SNR、热噪/芯片噪占比、量化底都完全不变, 只改 counts 标尺 —— 用于仿真
% rawdata 与实测对齐幅度(实测 counts 偏大但 SNR 一致的场合)。
if digGainDb ~= 0
    adcRawData = adcRawData * 10^(digGainDb/20);
    if printCfg, fprintf('数字增益       : %+g dB(纯数值缩放, 信号/噪声同比例, SNR 不变)\n', digGainDb); end
end

% ---- ADC 网格量化(物理 ADC 量化, 在 ADC 域/数字抽取前; 输出容器打包在最后) ----
if enableInt16Quant
    adcRawData = round(adcRawData);            % 在 ADC 的 2^(bits-1) 层网格上取整
end

% ---- 数字抽取链(chipOsr≥2): 逐级 半带FIR → ÷2(实数域, 每列 = 1 chirp 独立滤波) ----
% FIR 线性相位 → 群时延恒定 (L−1)/2 域内样本: 只加公共延迟, 不移距离 bin、不改
% 通道间相位差(全通道同一滤波器)。每级滤波在抽取【之前】→ 该级 ÷2 的折叠区被
% 本级阻带先压(与硬件"先滤后抽"一致)。chipOsr=1 时 nDStg=0, 本段整体跳过(=旧行为)。
for s = 1:nDStg
    for ii = 1:nRX
        adcRawData(:, :, ii) = filter(hbH{s}, 1, adcRawData(:, :, ii));
    end
    adcRawData = adcRawData(1:2:end, :, :);    % ÷2
end
% ---- 带边噪声均衡(输出域, 抽取后; 每列 = 1 chirp 独立滤波, 与半带级同口径) ----
% 把末级过渡带的噪底下垂拉平(设计见上方均衡块)。窗口沿 FIR 瞬态与半带级同性格
% (近距 bin 本就被高通压低 + 泄漏标定吸收), 对 350m 目标残余 droop ≈ −0.2dB。
if enableHbComp && nDStg > 0
    for ii = 1:nRX
        adcRawData(:, :, ii) = filter(bComp, 1, adcRawData(:, :, ii));
    end
end
if nDStg > 0
    fprintf('数字抽取链     : %d 级半带 FIR 完成, ADC 域 %d 样点/chirp → 输出 %d 样点/chirp', ...
            nDStg, chipOsr*nSamples, size(adcRawData, 1));
    if enableHbComp, fprintf(' + 带边噪声均衡 %d 抽头(噪底拉平)', numel(bComp)); end
    fprintf('\n');
end

% ---- 16bit 溢出异常模拟（补码环绕, 可开关; 位于量程检查之后、2D-FFT 之前） ----
% 上方量程检查按【未环绕】的真实峰值告警；本开关改变数据内容：超出 ±32768
% 的样点按 16 位补码环绕（+32768→−32768, +40000→−25536）——模拟寄存器按
% 补码解析、无饱和处理的真实数字链行为（不是模拟削顶）。环绕后时域波形
% "翻折到对面", 频域产生谐波+宽带杂散 → RD 底噪抬升（暗室角反 0° 现象）。
if enableInt16Wrap
    nWrap = nnz(abs(adcRawData) > adcFSmax);
    fprintf('16bit 溢出异常 : 开 —— %d 个样点超出 ±%d, 按补码环绕（+32768→−32768, 非削顶）\n', ...
            nWrap, adcFSmax);
    adcRawData = mod(adcRawData + adcFSmax, 2*adcFSmax) - adcFSmax;   % int16 补码环绕
end

% ---- ADC 整数量化（真实 2 字节数据口径, 可开关） ----
% 真实雷达 ADC 输出 2 字节整数(int16)。行业惯例: ADC 位数可以低于 16(如 12/14bit),
% 数据在 ADC 输出侧【左移 (16−bits) 位】装入 int16 —— 占满容器、下游满量程统一
% 按 ±32768 处理。左移是纯数值重贴标签, 不改变 SNR/动态范围, 但改变"1 个存储
% count 的物理电压"(存储 LSB = FS/2^15, 是 ADC 物理 LSB 的 1/2^(16−bits))——
% 标定换算时必须知道移位约定。
% 处理顺序(与硬件一致): ADC 域 2^(bits−1) 层网格取整(数字抽取前, 见上方) →
% 数字半带链 ÷2^nDStg → 16bit 补码环绕(输出寄存器) → 输出寄存器取整+左移装 int16。
% 量化附加噪声 = Δ/√12(Δ=ADC 物理 LSB), 比热噪低 ~35 dB(16bit)。
if enableInt16Quant
    adcShift = 16 - adcBits;                       % 左移位数(16bit→0)
    adcRawData = int16(round(adcRawData) * 2^adcShift);
    if printCfg, fprintf('ADC 整数量化   : 开 —— %dbit ADC 网格取整 + 左移 %d 位装 int16(存储满量程 ±32768)\n', ...
            adcBits, adcShift); end
end

% ==================== 时域绘图 ====================
% ADC 时域图不在本模块绘制(本模块只出数据): 主流程对校准遍/目标遍两次调用
% 各自用 adcTimeDraw(adcRawData, prm.fs, 图号名) 画 —— 校准遍读缓存跳过
% 射频时, 主流程用缓存数据同样能重画。见 adcTimeDraw.m。

% ---- 逐目标 IF 链路增益表(@差拍频点, LPF×两级HPF; 供 RCS 标定等下游修正已知滤形) ----
% filtGainAtR(R): 按距离查 IF 总增益(dB) —— fb = 斜率·2R/c(忽略多普勒频移, 远小于 fb),
% 供 RCS 估计对任意点云距离做同样的已知滤形修正(近距点落高通阻带)
if enableLpf, gA = @(fb) 1./sqrt(1 + (fb/lpfFc).^(2*lpfOrder));           % 模拟 LPF 增益
else,         gA = @(fb) ones(size(fb)); end
if enableHpf, gH = @(fb) (1./sqrt(1 + (hpf1Fc./max(fb,1)).^(2*hpf1Order))) .* ...
                         (1./sqrt(1 + (hpf2Fc./max(fb,1)).^(2*hpf2Order)));% 两级 HPF 级联增益
else,         gH = @(fb) ones(size(fb)); end
% 数字半带链增益(带内口径: f < fs/2 时各级直接相乘, 无折叠路径)
if nDStg > 0
    gD = @(fb) arrayfun(@(fq) prod(arrayfun(@(sI) hbMagOf(sI, fq), 1:nDStg)), fb);
else
    gD = @(fb) ones(size(fb));
end
filtGainAtR = @(Rq) 20*log10( gA(2*FreqSlope*Rq/c) .* gD(2*FreqSlope*Rq/c) ...
                              .* gH(2*FreqSlope*Rq/c) .* gC(2*FreqSlope*Rq/c) );

% 滤波链设计参数随 rfInfo 输出(图由主流程 filtChainDraw 在正确时序画, 模块不画图);
% 链路预算字段(字段名与 96/97 工程一致)供主函数 snrAnaCalc 审计: Pr 用 Pt/Gt/Gr/λ,
% PSD 用 Pn0_W_Hz = kT·NF(W/Hz, 不含带宽), 滤形用上面的 filtGainAtR
rfInfo = struct('filtGainDb', filtGainDb, 'fbHz', filtFbHz, 'filtGainAtR', filtGainAtR, ...
    'chipOsr', chipOsr, 'fsAdc', fsAdc, 'nDStg', nDStg, ...
    'hbTapsUsed', cellfun(@numel, hbH(1:max(nDStg,1))), 'B_noise', B_noise, ...
    'lpfFc', lpfFc, 'lpfOrder', lpfOrder, 'hbFc', hbFc, 'hbR', hbR, 'hbH', {hbH}, ...
    'hbCompOn', enableHbComp && nDStg > 0, 'compFreezeV', compFreezeV, ...
    'bComp', {bComp}, 'gC', {gC}, ...
    'lambda', lambda, 'Pt_dBm', Pt_dBm, 'Gt_dBi', Gt_dBi, 'Gr_dBi', Gr_dBi, ...
    'NF_dB', NF_dB, 'Pn0_W_Hz', Pn_W/B_noise, 'nTX', nTX, ...
    'chLossDb', 10*log10(mean((TX_A(:)*RX_A(:).').^2, 'all')), ...   % 通道锥削均值损耗(dB): 逐天线损耗差(TX_A/RX_A)的 64 通道平均功率 ≈−0.8, 供 snrAnaCalc 显式列项
    'pllTrackOn', enablePllTrack, 'pnTrackCoef', pnTrackCoef, ...
    'hopErrOn', enableHopErr, 'hopErrHz', hopErr_Hz);

end
