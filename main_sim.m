%% ========================================================================
%  main_sim.m —— 8T8R 车载毫米波雷达 DDMA 仿真主流程(两遍测量, 与实际工程一致)
%  实际工程流程: ①校准场放角反(+检查目标) → 测各项校准数据;
%               ②撤场布目标场景 → 冻结①的校准数据跑处理链看性能。
%  两遍共用同一套射频配置(开关/旋钮全在 rf_sim_core.m 顶部唯一定义,
%  改一处两遍同时生效), 除目标外完全一致; 校准量(cCal/AntCaliData/rcsCfg)
%  全部冻结于第一遍, 第二遍只消费不重算。
%
%  链路(按执行顺序, 图号=建图顺序, PNG 存工程目录供 PPT):
%   ---- 静态设计图(每次运行相同) ----
%    天线方向图    genAntPattern.m     图1~6
%    阵列几何      plotArrayGeom.m     图7~10
%    波形图        fmcwWaveDraw.m      图11   FMCW 频率-时间
%   ---- 第一遍: 校准测量(校准场仅放角反) ----
%    射频仿真      rf_sim_core.m       图13   ADC 时域(校准遍)
%    加窗 2D-FFT   win2dfft.m          图14/14  RD 谱(校准遍)
%    噪底估计      estNoiseFloor.m     图16/16  逐距离门(校准遍)
%    0° 校准       calib0deg.m         图18   校准检查         → cCal, calibInfo
%    导向矢量校准  genAntCaliData.m    (无图, cCal × 理论流型) → AntCaliData, rowLut
%    检测链(静默)  积累/CFAR/DBF       (无图, 供角反仿真测量幅度)  → tgtCal, detCal
%    RCS 标定      rcsCalibrate_pipe   图19   标定结果         → rcsCfg(K0)
%   ---- 第二遍: 目标场景测量(其余目标, 图号顺延) ----
%    射频仿真      rf_sim_core.m       图20   ADC 时域
%    加窗 2D-FFT   win2dfft.m          图21/21  RD 谱
%    噪底估计      estNoiseFloor.m     图23/23
%    RX 非相参     rxNonCoh_pipe.m     图25   8RX 功率平均(选阶前初查)
%    积累+选阶     jointAccum_pipe.m   图26   非相参/相参      → accumCfg
%    两级 CFAR     targetDetect_pipe.m 图27                     → detCfg, tgtOut
%    DBF 测角      doaEstimate_pipe.m  图28                     → tgtOut(含方位/俯仰)
%    OMP 超分辨    ompRefine_pipe.m    图29/28                  → tgtOut(+第二目标)
%    DML 超分辨    dmlRefine_pipe.m    图31/30                  → tgtOut(+第二目标)
%    真值/对比图   targetTruthDraw.m   图33  compareAziRDraw  图34
%    RCS           rcsCompute_pipe.m   图35   (用第一遍 rcsCfg) → tgtOut(+rcsDb)
%    链路预算审计   snrAnaCalc.m + cohStatCalc.m + cohGainDbf.m  无图, "解析 vs 仿真测量"打印
%  约定: 波形类参数只在 prm、维度类只在 dims 定义一次, 子函数形参接收;
%  变量名与 参考\ 流水线保持一致。
%  ========================================================================
clc; close all;
% 钩子变量(批量案例 runCases*.m / 门限标定 thrSweep_ompdml.m 注入, 平时不存在)
clearvars -except targetsOverride saveDBFStage saveFigDir skipFigs batchState;
rng(4); % 固定随机种子(照 96 修法): 种子只决定单实现统计起伏, 不改物理 —— 泄漏
         % 路径相位/闪烁噪声/PN 实现/量化误差等每次运行完全一致, 结果可复现、
         % 改动前后可逐位对比(换种子 = 换一次统计起伏实现, 结论不应翻转; 96 的
         % seed_sweep 多种子诊断法同样适用)。校准缓存存有 rng 快照, 命中时恢复
         % 后目标遍与全量跑仍完全一致。

%% ---- 检测配置(仿参考 AdcAndParam.mat, 变量名与参考一致) ----
load('SimAdcParam.mat');       % cfarCfg: 逐距离门 CFAR/峰判配置

%% ---- 波形与阵列参数(唯一定义处) ----
nSamples = 1024;               % 每 chirp 实采样点数。【硬件内存上限, 勿改】
                               %   推论(不可绕过): fb(350m)/奈奎斯特 = 350/(Ns·c/4B)
                               %   = 97.5%, 350m 目标永远贴着混叠边 —— 想让 LPF 缓降
                               %   只能降 B(牺牲 ΔR)或接受带边抑制变弱, 见射频模块
                               %   LPF 段注释的三选项
nChirps  = 512/2;                % 每帧 chirp 数(多普勒 FFT 长度)【硬件内存上限, 勿改】

fsAdc   = 160e6;                % 芯片 ADC 实采样率 (Hz) —— 真实芯片 160M 链; 厂商 2026-10 口径
                               %   的"40M"= 降采样后的输出数据率, 不是 ADC 物理采样率
chipOsr = 4;                    % [旋钮] 降采样倍数(1 | 2 | 4, 2 的幂)——主函数与射频函数都要用
                               %   (照 96 修法上提主函数, 唯一定义处在此、经 prm 注入射频核):
                               %   射频核用它设计数字半带链(先滤后抽, 每级 ÷2)并降采样,
                               %   本函数用它推算真正的数据率与雷达各项指标
fs = fsAdc/chipOsr;             % 降采样后输出数据率 (Hz) = 40M —— 下游全部联动:
                               %   Ts = nSamples/fs, 斜率 = B/Ts, 奈奎斯特 = fs/2,
                               %   Rmax = c·nSamples/(4B)(只由 nSamples 与 B 决定, 与 fs 无关)
B  = 214e6;                % 采样窗 Ts 内扫过的带宽 (Hz); ΔR = c/2B 与 Rmax 不随 fs 变
Tc = 37e-6;                    % chirp 周期 = Ts + 保护时间 (s); Vres/Vmax 随 Tc 自动变
cr = 2;                        % 步进频占比 = 总跳频带宽 / B
stepDir = -1;                   % [旋钮] 步进频方向: +1 上跳频(每 chirp 起始频率比上一个高, 现行波形)
winType  = 'Hanning';           % [旋钮] 2D-FFT 加窗: 'Hanning' | 'SVA'(自适应旁瓣抵消,
                               %   两遍共用; winStat 随窗类型分档标定, 见链路预算审计段)
                               %   / −1 下跳频(每 chirp 起始频率比上一个低)。
                               %   联动: 射频链路(f_m 序列)、步进频耦合方向(谱线 bin = gate_D −方向·gate_R)、
                               %   校准/RCS/案例生成的 dopC 公式全部经 prm 自动跟随
fStart = 76.6e9;               % 扫频起始频率 (Hz)
prm = struct('fs', fs, 'fsAdc', fsAdc, 'chipOsr', chipOsr, 'B', B, 'Ts', nSamples/fs, ...
             'Tc', Tc, 'cr', cr, 'fStart', fStart, 'stepDir', stepDir);
c0 = physconst('lightspeed');
fCenter = fStart + (1+cr)*B/2;   % 全帧中心频率(扫频+步进频全跨度中点)
lambda = c0 / fCenter;                        % 波长(算好传入, 子函数不自算)
TX_A_dB = [0 -0.6 -0.25 -1.0 -0.4 -0.15 -0.75 -0.3];   % TX 通道幅度差(dB, 最大=0)
RX_A_dB = [-0.35 -0.7 -0.1 0 -0.9 -0.45 -0.2 -0.55];   % RX 通道幅度差(dB)
TX_A = 10.^(TX_A_dB/20);       % 线性(最大=1; 通道只会差不会好)
RX_A = 10.^(RX_A_dB/20);
FFT_d0 = 1.865e-3;             % 阵元栅格间距 (m), 天线坐标 = 整数栅格 × FFT_d0
% ---- [旋钮] 样件版本: B1 / B2 布阵不同, 其余参数(TX_A/RX_A 通道差、编码、波形)共用
%      (B2 若另有实测通道差/隔离度, 再按样件分档)。切换样件会自动重跑校准(缓存签名含阵列)。
sampleVer = 'B2';
switch upper(sampleVer)
    case 'B1'                  % B1 样件布阵
        antenna_t = [[30;0;3],[27;0;3],[20;0;3],[33;0;9],[24;0;0],[14;0;3],[7;0;3],[0;0;3]] * FFT_d0;
        antenna_r = [[26;0;15],[29;0;18],[21;0;18],[33;0;18],[12;0;18],[4;0;18],[11;0;0],[0;0;18]] * FFT_d0;
    case 'B2'                  % B2 样件布阵(2026-10-29 用户提供)
        antenna_t = [[6;0;14],[9;0;14],[16;0;14],[0;0;8],[12;0;20],[22;0;14],[29;0;14],[36;0;14]] * FFT_d0;
        antenna_r = [[9;0;3],[6;0;0],[14;0;0],[2;0;0],[23;0;0],[31;0;0],[27;0;21],[35;0;0]] * FFT_d0;
    otherwise
        error('sampleVer 只支持 ''B1'' | ''B2''');
end
fprintf('样件           : %s 布阵(TX %d 元 / RX %d 元)\n', upper(sampleVer), ...
        size(antenna_t, 2), size(antenna_r, 2));
encodeSeq = [1 1 1 1 0 0 1 0 1 0 1 0 0 1 0 0]; % DDMA 编码(16 槽 8 用)
ddmShift  = find(encodeSeq == 1) - 1;         % 各 TX 码槽位(0 起)
prm.nGroups = numel(encodeSeq);               % DDMA 分组数, 射频模块经 prm 读取

% ---- 波形真实占用频带(法规核查): 采样窗外波形仍在扫频 ----
tWaitSample = 3e-6;                % [旋钮] 采样前等待时间 (s)
tHoldSample = 2e-6;                % [旋钮] 采样后保持时间 (s)
assert(tWaitSample + tHoldSample <= Tc - prm.Ts, ...
       '等待+保持 (%.1f us) 超出保护时隙 (%.1f us)', ...
       (tWaitSample + tHoldSample)*1e6, (Tc - prm.Ts)*1e6);
chirpSlope  = B/prm.Ts;                              % 扫频斜率 (Hz/s)
% 帧射频范围随步进方向: 上跳 [f0−w, f0+(1+cr)B+h]; 下跳 [f0−crB−w, f0+B+h]
fStartWave  = fStart + min(0, stepDir*cr*B) - chirpSlope*tWaitSample;   % 全帧最低(含等待段延伸)
fStopWave   = fStart + B + max(0, stepDir*cr*B) + chirpSlope*tHoldSample; % 全帧最高(含保持段延伸)
fprintf('波形真实频带   : %.6f ~ %.6f GHz (真实占用带宽 %.2f MHz = 全帧扫频(1+cr)·B + 等待%.1fus/保持%.1fus 斜率延伸)\n', ...
        fStartWave/1e9, fStopWave/1e9, (fStopWave - fStartWave)/1e6, ...
        tWaitSample*1e6, tHoldSample*1e6);
if fStartWave >= 76e9 && fStopWave <= 81e9
    fprintf('   在 76~81 GHz 车载典型许可频带内\n');
else
    fprintf('   !!! 超出 76~81 GHz 车载典型许可频带, 请核对法规 !!!\n');
end

% ---- 处理链公共维度(唯一定义处, 子函数经 dims 接收) ----
rngPts     = nSamples;                          % 距离维采样点数
nV         = rngPts/2 + 1;                      % 实采样共轭对称, 有效距离 bin = 半数 + DC
dopBinNum  = nChirps;                           % 速度维 bin 数
grpTotal   = numel(encodeSeq);                  % DDM 分组数 (=16)
subBins    = dopBinNum/grpTotal;                % 每子带 bin 数 (=32)
chanTotal  = numel(ddmShift)*size(antenna_r,2); % 通道数 = TX×RX (=64)
aziFftLen  = 128;                               % 方位维 FFT 点数
nLongRow = [];   % 最长行虚拟阵元数(方位维通道数), A_manifold 生成后统计回填
dims = struct('rngPts', rngPts, 'nV', nV, 'dopBinNum', dopBinNum, ...
              'grpTotal', grpTotal, 'subBins', subBins, 'chanTotal', chanTotal, ...
              'aziFftLen', aziFftLen, 'nLongRow', nLongRow, ...
              'rngUse', nV - 1, 'B', B, ...
              'DeltaR', c0/(2*B), 'Vres', c0/(2*fCenter*Tc*nChirps));
%% ---- 目标场景(结构体数组; targets(1) 固定为 0° 校准角反, 勿动) ----
% gate_R/gate_D: 距离/速度 bin, 与 RD 图对齐; gate_D 可超 512(谱线折叠回同 bin)。
% 速度分解: gate_D×Vres = 法向速度(产生多普勒, 可测); vTan = 切向 (m/s, 可为
% 小数, 测不到, 留给跟踪推测)。表列 dir: 0=Y 向行驶 / 1=X 向横穿, vTan 按
% (gate_D, azi, ele) 矢量分解自动算出。
% 案例对分工: 角差 ≥5° 考 OMP, 1.5~6° 考 DML, 5~6° 两模块重叠;
% 同速 bin 两做法: n=0 完全同速 / n=1 真速差 512×Vres=32.7 m/s 折叠同 bin
% (须选低速组合, 否则帧内距离走动糊 bin)。
% 表列: [gate_R, gate_D, azi, ele, rcsDb, rcsPh, dir, isCorner]
% rcsPh = RCS 初相(deg): 实际目标 RCS 是复数(σ∠φ), 两个同距同速目标即使幅度
%   相同, 初相不同 → 阵列干涉图样不同 → 影响 DBF 融合/OMP 残差/DML 分离。
%   不需要时全填 0(等效原有行为)。
% ★ 表内目标只进【第二遍目标场景】; 校准遍场景固定为"角反 + 检查目标",
%   与本表无关(见下方场景拆分)。
if exist('targetsOverride', 'var') && ~isempty(targetsOverride)   % 钩子: 外部场景注入
    targets = targetsOverride;
    fprintf('目标场景       : 采用外部注入 targetsOverride(%d 个目标 —— 门限标定/验证用)\n', numel(targets));
else
targets(1) = struct('gate_R', 5, 'gate_D', 0, 'azi', 0, 'ele', 0, 'rcsDb', -5, ...
                    'rcsPhDeg', 0, 'vTan', 0, 'isCorner', 1);
T = [ ...
% ---- 目标表(每行一个目标; rcsPh 列 = 初相 deg, 默认 0) ----
 150/dims.DeltaR   0  0   0   10   0  0  0];   % 300m 目标
for k = 1:size(T, 1)
    vRadK = T(k,2)*dims.Vres;                          % 法向速度 (m/s)
    if T(k,7) == 0                                     % dir=0: Y 向行驶
        vTanK = -vRadK*tand(T(k,3))/cosd(T(k,4));
    else                                               % dir=1: X 向横穿
        vTanK =  vRadK/(tand(T(k,3))*cosd(T(k,4)));
    end
    targets(k+1) = struct('gate_R', T(k,1), 'gate_D', T(k,2), 'azi', T(k,3), ...
                          'ele', T(k,4), 'rcsDb', T(k,5), 'rcsPhDeg', T(k,6), ...
                          'vTan', vTanK, 'isCorner', T(k,8));
end
clear T;
end                                                 % else(未注入 targetsOverride)

% ---- 场景拆分(与实际工程一致): 校准遍 / 目标遍 ----
% 校准遍 = 仅角反(targets(1)): 校准场里只有角反, 没有其它目标;
% 目标遍 = 除角反外的全部目标(下方处理链全部按此场景)。
targetsCal = targets(1);                   % 第一遍(校准)场景 = 角反单独
targets    = targets(2:end);               % 第二遍(目标)场景, 变量名沿用 targets
if isempty(targets)
    error('目标场景为空: 请在目标表 T 加目标, 或注入 targetsOverride');
end
fprintf('场景拆分       : 校准遍 1 目标(角反), 目标遍 %d 目标\n', numel(targets));


% ---- RAM 核查行(照 96 修法): 全 RX 帧数据量, 每次运行打印盯住 ----
% 口径: 降采样后 nSamples×nChirps×nRX×2B(int16/样点) —— 8 路 ADC 同时采样必须全落
% 芯片 RAM(用户口径), 即上方 nSamples/nChirps 两行【硬件内存上限, 勿改】的核定依据;
% 改点数/chirp 数前先看这行是否还在预算内。
ramKB = nSamples*nChirps*size(antenna_r,2)*2/1024;    % 全 RX 帧数据量 (KB)
fprintf('RAM 核查        : 全RX帧 %d样点×%dchirp×%dRX×2B(int16) = %.0f KB = %.2f MB —— nSamples/nChirps【硬件内存上限】核定口径\n', ...
        nSamples, nChirps, size(antenna_r,2), ramKB, ramKB/1024);

% ---- cfarCfg 逐距离门表对齐当前量程(SimAdcParam.mat 按 512 距离门配置) ----
% 表长与 rngUse 不一致时: 超长裁前 rngUse 行(量程缩小时取前段配置);
% 不足则整表重复铺满(量程放大时沿用配置图案)。不改此段则 targetDetect
% 里 floorPow(×rngUse) 与 cfar 表(×512) 元素运算会维度不兼容报错。
cfarFlds = {'cfarCoarseList', 'cfarThreList', 'cfarType', 'cfarPeakList'};
for qF = 1:numel(cfarFlds)
    tbl = cfarCfg.(cfarFlds{qF});  tbl = tbl(:);
    if numel(tbl) < dims.rngUse, tbl = repmat(tbl, ceil(dims.rngUse/numel(tbl)), 1); end
    cfarCfg.(cfarFlds{qF}) = tbl(1:dims.rngUse);
end

%% ---- 天线方向图(图1~6; 逐 TX/RX 复数方向图, 暂未接入任何模块) ----
% 方位/俯仰各 ±90°(1° 网), 峰值 15 dBi。方位不组阵(±70° 仍 0 dBi, 无副瓣);
% 俯仰组阵(±20° 到 0 dBi, 副瓣 −25 dB 等电平、每侧 5 个)。
% 16 个天线不再共用一张理想图: genAntPattern 用固定种子的随机模型叠加了
% 通道间不一致性 —— 幅度偏置 ±0.4dB、形状起伏 ±0.25/±0.15dB、相位偏置 ±10°、
% 相位斜率(指向微偏)与 ±5° 起伏, 模拟实测方向图的逐天线差异(换实测图时逐页替换)。
% 图1~4 看幅度(热图 + 每天线切面), 图5~6 看远场相位(复数方向图的另一半,
% 相位一致性影响测角, 这正是实测要采导向矢量而非用理论值的原因)。
% ★ 方向图仍在检查阶段: 确认形状合理前, 其他模块一律不要读 antPat。
antPat = genAntPattern(size(antenna_t, 2), size(antenna_r, 2));
[figPatTx, figPatRx, figPatTxCut, figPatRxCut, figPatTxPh, figPatRxPh] = antPatDraw(antPat);

%% ---- 图7/图8: 阵列几何(设计原始整数栅格, 居中前画, 便于对照设计查错) ----
[figPhys, figVirt] = plotArrayGeom(antenna_t, antenna_r, arrayManifold(antenna_t, antenna_r), FFT_d0);

% ---- 阵列整体居中(TX/RX 同一平移量): 只移相位参考, 不改任何相对位置 ----
antShift = [(min([antenna_t(1,:) antenna_r(1,:)]) + max([antenna_t(1,:) antenna_r(1,:)]) )/2; 0; ...
            (min([antenna_t(3,:) antenna_r(3,:)]) + max([antenna_t(3,:) antenna_r(3,:)]) )/2];
antenna_t = antenna_t - antShift;
antenna_r = antenna_r - antShift;

%% ---- 阵列流型(居中后坐标; 供校准/积累/DOA 用) ----
A_manifold = arrayManifold(antenna_t, antenna_r);   % [3 × 64], TX 内循环/RX 外循环:
                                                   % TX1RX1..TX8RX1, TX1RX2..
zRowsU = unique(A_manifold(3, :));
nInRowArr = zeros(numel(zRowsU), 1);
for rr = 1:numel(zRowsU), nInRowArr(rr) = sum(A_manifold(3, :) == zRowsU(rr)); end
nLongRow = max(nInRowArr);                  % 最长一行阵元数(当前布阵 = 36)
dims.nLongRow = nLongRow;
fprintf('阵列统计       : 最长一行 %d 个虚拟阵元(共 %d 行)\n', nLongRow, numel(zRowsU));
% 图9/图10: 居中后对照图, 应是图7/图8 的纯平移(可含半格偏移)
[figPhysC, figVirtC] = plotArrayGeom(antenna_t, antenna_r, A_manifold, FFT_d0);
set(figPhysC, 'Name', 'fig9_PhysArrayC');
set(figVirtC, 'Name', 'fig10_VirtArrayC');

%% ---- FMCW 波形图(图11: 频率-时间; 全帧 512 chirp + 前 10 个放大, 含等待/保持段) ----
figWave = fmcwWaveDraw(prm, dims, tWaitSample, tHoldSample);


%% ================== 第一遍测量: 校准(仅角反 → 各项校准数据) ==================
% 顺序同实际工程: 校准场只放角反跑射频 → RD 谱/噪底(校准遍自测量) → 0°校准
% (cCal) → 导向矢量校准表(AntCaliData) → 角反点云检测(RCS 标定要角反仿真测量
% 幅度, 不画积累/检测/DBF 图) → RCS 标定常数(rcsCfg)。射频开关/旋钮全在
% rf_sim_core.m 顶部, 与第二遍共用同一份 —— 改一处两遍同时生效。
fprintf('\n========== 第一遍测量: 校准(仅角反) ==========\n');

% ---- 校准缓存: 角反射频结果(adcCal)复用, 免每次重跑 ~11 s 仿真 ----
% calRerun = true  : 强制重跑角反射频并刷新缓存(改射频后想看新校准遍打印时用)
% calRerun = false : 优先读缓存; 缓存不存在或签名不匹配时自动重跑并保存
% 签名覆盖影响 adcCal 的全部输入(射频模块源码哈希 + 波形/阵列/通道/角反布设),
% 任一改动自动判失效 —— 缓存永远与当前配置一致, 不会用到过期校准数据。
calRerun     = false;                            % [旋钮] 校准缓存开关(见上)
calCacheFile = fullfile(fileparts(mfilename('fullpath')), 'calibCache.mat');
rfSrcFile    = fullfile(fileparts(mfilename('fullpath')), 'rf_sim_core.m');
rfD          = dir(rfSrcFile);
try                                          % 射频模块源码 MD5(JVM 不可用退化为 mtime+size)
    jMD = java.security.MessageDigest.getInstance('MD5');
    jMD.update(uint8(fileread(rfSrcFile)));
    rfHash = lower(reshape(dec2hex(typecast(uint8(jMD.digest), 'uint8')), 1, []));
catch
    rfHash = '';
end
calSig = struct('rfHash', rfHash, 'rfDate', rfD.datenum, 'rfBytes', rfD.bytes, ...
                'prm', prm, 'nSamples', nSamples, 'nChirps', nChirps, ...
                'TX_A', TX_A, 'RX_A', RX_A, ...
                'antenna_t', antenna_t, 'antenna_r', antenna_r, ...
                'ddmShift', ddmShift, 'encodeSeq', encodeSeq, 'targetsCal', targetsCal);
cacheOk = false;
if ~calRerun && exist(calCacheFile, 'file') == 2
    S = load(calCacheFile);
    if isfield(S, 'sig') && isequal(S.sig, calSig)
        adcCal = S.adcCal;  targetsCal = S.targetsCal;
        rng(S.rngState);    % 恢复校准射频结束时的随机流 → 目标遍与全量跑完全一致
        cacheOk = true;
        fprintf('校准缓存       : 命中(生成于 %s) —— 跳过角反射频仿真, 校准遍图与校准模块照常重跑\n', ...
                char(S.savedAt, 'yyyy-MM-dd HH:mm:ss'));
        if isfield(S, 'filtChain')
            hFilt = filtChainDraw(S.filtChain, prm.fs, 2*S.filtChain.chipOsr*prm.fs);  % fig12(显式编号)
        end
    else
        fprintf('校准缓存       : 签名不匹配(射频配置/波形/阵列/角反布设已改) —— 自动重跑并刷新\n');
    end
end
if ~cacheOk
    tRfC = tic;
    % printCfg=false: 射频配置与目标遍完全相同, 配置块不重复打印(目标表仍打印)
    [adcCal, rfInfoCal] = rf_sim_core(nSamples, nChirps, TX_A, RX_A, antenna_t, antenna_r, ...
                                      ddmShift, prm, targetsCal, [], false);
    fprintf('耗时           : 射频模块(校准遍) %.1f s\n', toc(tRfC));
    filtChain = struct('chipOsr', rfInfoCal.chipOsr, 'fsAdc', rfInfoCal.fsAdc, ...
                       'nDStg', rfInfoCal.nDStg, 'B_noise', rfInfoCal.B_noise, ...
                       'lpfFc', rfInfoCal.lpfFc, 'lpfOrder', rfInfoCal.lpfOrder, ...
                       'hbFc', rfInfoCal.hbFc, 'hbR', rfInfoCal.hbR, 'hbH', {rfInfoCal.hbH}, ...
                       'hbCompOn', rfInfoCal.hbCompOn, 'compFreezeV', rfInfoCal.compFreezeV, ...
                       'bComp', {rfInfoCal.bComp});
    S = struct('adcCal', adcCal, 'targetsCal', targetsCal, 'sig', calSig, ...
               'rngState', rng, 'savedAt', datetime('now'), 'filtChain', filtChain);   % rng 快照供缓存命中时恢复
    save(calCacheFile, '-struct', 'S');
    fprintf('校准缓存       : 已保存 %s(下次同配置直接读取)\n', calCacheFile);
end
if ~exist('hFilt', 'var')   % 缓存命中路径已在上面画过; 未命中路径在此补画(位置都在 fig13 之前)
    hFilt = filtChainDraw(filtChain, prm.fs, 2*filtChain.chipOsr*prm.fs);  % fig12(显式编号)
end
hAdcCal = adcTimeDraw(adcCal, prm.fs, 'fig13_ADC_time_cal');   % 图13(从数据画, 不依赖射频模块)

%% ---- 加窗 2D-FFT + 噪底(校准遍; 图14~16) ----
fftCal = win2dfft(adcCal, winType);
[fftCal, noiseCal] = estNoiseFloor(fftCal);

%% ---- 0° 校准(角反: 已知距离/速度 + 码位 + 耦合 → 各 TX bin; 图18) ----
R1 = targetsCal(1).gate_R * 3e8/(2*prm.B);   % 角反距离 = gate_R × ΔR
v1 = 0;                                      % 角反静止
[cCal, calibInfo] = calib0deg(fftCal, encodeSeq, nSamples, nChirps, prm, R1, v1, ...
                              antenna_t, antenna_r, lambda);

%% ---- 导向矢量校准表(理论流型 × 0°校准 cCal; 冻结于校准遍, 两遍共用) ----
% 俯仰网格由无模糊范围决定(λ/2×俯仰行间距), 向内取整 1°: 超出部分是栅瓣副本
dzRow  = spacingGCD(A_manifold(3, :));                   % 俯仰行基线间距(间距 GCD)
eleUa  = asind(min(1, lambda/(2*dzRow)));                % 俯仰无模糊范围 (deg)
eleRange = [-floor(eleUa), floor(eleUa), 1];
fprintf('DOA 俯仰网格    : ±%d° (无模糊 ±%.1f°, 俯仰行间距 %.3f mm)\n', ...
        floor(eleUa), eleUa, dzRow*1e3);
[AntCaliData, rowLut] = genAntCaliData(A_manifold, FFT_d0, lambda, cCal, [-75 75 1], eleRange);

% ---- 校准遍检测链(积累→CFAR→DBF, 为 RCS 标定提供角反仿真测量幅度; 静默不画图) ----
accumCal = jointAccum_pipe(fftCal, ddmShift, cCal, A_manifold, FFT_d0, dims, false);
[~, detCal, tgtCal] = targetDetect_pipe(noiseCal, accumCal, cfarCfg, dims, false);
[tgtCal, detCal]    = doaEstimate_pipe(detCal, tgtCal, AntCaliData, rowLut, lambda, dims, false);

%% ---- RCS 标定(角反 SNR 比值法; 图19) ----
% σ_dB = SNR + 40log10(R) + K0; SNR = |α|²/该门噪底, 与 IF 滤形无关
% (信号与噪底同 bin 同滤形, 无需知道滤波器形状, 只用噪底这一链路自测量)
rcsCfg = rcsCalibrate_pipe(tgtCal, detCal, AntCaliData, rowLut, dims, prm, targetsCal, calibInfo, noiseCal);

%% ================== 第二遍测量: 目标场景(消费第一遍的校准数据) ==================
% 不再重跑 calib0deg/RCS 标定 —— 积累用 cCal、测角用 AntCaliData、RCS 用 rcsCfg,
% 全部来自第一遍; 两遍射频配置唯一(均在 rf_sim_core.m 顶部), 仅目标不同。
fprintf('\n========== 第二遍测量: 目标场景(%d 个目标) ==========\n', numel(targets));
tRf = tic;
[adcRawData, rfInfo] = rf_sim_core(nSamples, nChirps, TX_A, RX_A, antenna_t, antenna_r, ddmShift, prm, targets);
fprintf('耗时           : 射频模块(目标遍) %.1f s\n', toc(tRf));
hAdcScn = adcTimeDraw(adcRawData, prm.fs, 'fig20_ADC_time');        % 图20(从数据画, 不依赖射频模块)

%% ---- 加窗 2D-FFT('Hanning' | 'SVA'; 图21/21) ----
fftCfg  = win2dfft(adcRawData, winType);

%% ---- 逐距离门噪底估计(图23/23) ----
[fftCfg, noiseCfg] = estNoiseFloor(fftCfg);

%% ---- RX 非相参积累(8RX 线性功率平均; 图25, DDMA 选阶之前的初查图) ----
% 实测排查惯例(实机排障流程, 仿真沿用): 先不做解模糊选阶, 直接把 8 个 RX 合成一张 RD 图看数据"活不活"
% (目标在不在/泄漏形状/异常条纹), 确认后再进选阶与相参积累。
% 累加方式与 TXRX 非相参一致: 线性域功率平均后转 dB(非 dB 域平均)。
figRxNc = rxNonCoh_pipe(fftCfg, dims);

%% ---- DDM 解模糊选阶 + 非相参/相参积累(64 元流型方位 FFT; 图26; 校准表 cCal 来自校准遍) ----
accumCfg = jointAccum_pipe(fftCfg, ddmShift, cCal, A_manifold, FFT_d0, dims);

%% ---- 两级 CFAR(粗检 DACFAR + 峰态精判; 图27) ----
[accumCfg, detCfg, tgtOut] = targetDetect_pipe(noiseCfg, accumCfg, cfarCfg, dims);
fprintf('两级CFAR检测   : 粗检 %d 点, 确认目标 %d 个\n', ...
        size(detCfg.coarseDetList, 1), numel(tgtOut));
if isempty(tgtOut)
    error(['两级CFAR确认 0 个目标 —— 目标表的目标可能已全部出带。', ...
           '当前量程 = rngUse×ΔR = %d×%.3f m = %.0f m(nSamples=%d, 斜率=%.2f MHz/μs), ', ...
           '请把目标表 gate_R 压进量程, 或增大 nSamples/fs。', ...
           '(出带目标差拍超奈奎斯特折叠成鬼影, 不构成有效检测)'], ...
           dims.rngUse, dims.DeltaR, dims.rngUse*dims.DeltaR, nSamples, prm.B/prm.Ts/1e12);
end

%% ---- DBF 测角(方位 DBF + 子阵俯仰 DBF; 图28; 导向矢量校准表来自校准遍) ----
[tgtOut, detCfg] = doaEstimate_pipe(detCfg, tgtOut, AntCaliData, rowLut, lambda, dims);
% DBF 阶段点云快照(OMP/DML 会改写方位, 图34 四联对比用原值)
stageDbf = struct('azi', [tgtOut.azimuth], 'rng', [tgtOut.range], ...
                  'z', [tgtOut.range].*sind([tgtOut.elevation]));
% 钩子: 预置 saveDBFStage(文件名)时存 DBF 阶段输入, 供 OMP/DML 离线扫门限
if exist('saveDBFStage', 'var') && ~isempty(saveDBFStage)
    save(saveDBFStage, 'tgtOut', 'detCfg', 'AntCaliData', 'rowLut', 'dims', 'prm', 'targets');
    fprintf('DBF 阶段快照   : 已保存 %s (OMP/DML 门限离线扫描用)\n', saveDBFStage);
end

%% ---- OMP 角度超分辨(大间隔强弱对; 图29 过程图 / 图30 点云) ----
% seedR/seedD 指定画过程图的点(按谱线实际位置含耦合); 换目标改这两行
thr2ndDb = 12;                                   % [旋钮] OMP 第二目标门限 (dB, 残差峰高出噪底中位)
                                                 % 标定: 8~14 dB 命中率同平台(89%), 16 掉到 67%;
                                                 % 10dB 强弱对残差余量天然 11~15dB, >12 会漏弱目标;
                                                 % 单目标误分裂率对门限不敏感, 靠输出门限压制
seedTgtO = min(3, numel(targets));               % OMP 过程图目标(目标表第 3 个; 不足取最后一个)
dopC3    = mod(targets(seedTgtO).gate_D + prm.stepDir*prm.cr*prm.B/(nChirps-1) * (2*targets(seedTgtO).gate_R*c0/(2*prm.B))/c0 * nChirps, nChirps);
seedR    = targets(seedTgtO).gate_R + 1;         % 距离 bin (1 起)
seedD    = round(dopC3) + 1;                     % 多普勒 bin (1 起, 含耦合)
[tgtOut, detCfg, ompCfg] = ompRefine_pipe(tgtOut, detCfg, AntCaliData, rowLut, ...
                                          dims, thr2ndDb, seedR, seedD);
fprintf('OMP 角度超分辨 : 判定出OMP第二目标 %d 个 (判定门限 %d dB), 另 %d 个超出输出范围(±20°/±25m)删除, 处理门限外跳过 %d 点\n', ...
        ompCfg.n2nd, thr2ndDb, ompCfg.n2ndDrop, ompCfg.nGate);
% OMP 阶段点云快照(图34 用; DML 还会改写方位)
stageOmp = struct('azi', [tgtOut.azimuth], 'rng', [tgtOut.range], ...
                  'z', [tgtOut.range].*sind([tgtOut.elevation]), ...
                  'isO2', logical([tgtOut.isOMP2nd]));

%% ---- DML 角度超分辨(近角对; 图31 过程图 / 图32 点云) ----
% OMP 第二目标不再进 DML; 不通过时保留 DBF 角(宁可不分, 不分错)
thrDmlDb = 3;                                    % [旋钮] DML 第二目标门限 (dB, 双角代价改善量)
                                                 % 标定: 1.5~5 dB 命中率平台, 抬门限只减误分裂;
                                                 % 实测近角对改善量可低至 3.3dB, 取 3 保灵敏度;
                                                 % 实测数据杂点多时可升到 3.5(标定证据支持)
tDml = tic;
seedTgtD = min(5, numel(targets));               % DML 过程图目标(目标表第 5 个; 不足取最后一个)
dopC     = mod(targets(seedTgtD).gate_D + prm.stepDir*prm.cr*prm.B/(nChirps-1) * (2*targets(seedTgtD).gate_R*c0/(2*prm.B))/c0 * nChirps, nChirps);
seedR    = targets(seedTgtD).gate_R + 1;         % 距离 bin (1 起)
seedD    = round(dopC) + 1;                      % 多普勒 bin (1 起, 含耦合)
[tgtOut, detCfg, dmlCfg] = dmlRefine_pipe(tgtOut, detCfg, AntCaliData, rowLut, dims, thrDmlDb, seedR, seedD);
fprintf('DML 角度超分辨 : 判定出第二目标 %d 个 (门限 %g dB); 门限跳过 %d 点; 实做 DML %d 点 (耗时 %.1f s)\n', ...
        dmlCfg.nDML2nd, thrDmlDb, dmlCfg.nGate, dmlCfg.nDone, toc(tDml));

%% ---- 真值目标图(图33; 与点云图同口径, 对比看漏检/虚假点) ----
figTruth = targetTruthDraw(targets, prm, dims);

%% ---- 方位角-距离四联对比(图34; 真值/DBF/OMP/DML 同轴同色标) ----
stageDml = struct('azi', [tgtOut.azimuth], 'rng', [tgtOut.range], ...
                  'z', [tgtOut.range].*sind([tgtOut.elevation]), ...
                  'isO2', logical([tgtOut.isOMP2nd]), 'isD1', logical([tgtOut.dmlPass]), ...
                  'isD2', logical([tgtOut.isDML2nd]));
figCmpAzi = compareAziRDraw(targets, dims, stageDbf, stageOmp, stageDml);

%% ---- RCS 估计(逐点算 RCS 写回点云; 图35; 标定常数 rcsCfg 来自第一遍校准) ----
[tgtOut, rcsCfg, figRcs] = rcsCompute_pipe(tgtOut, detCfg, AntCaliData, rowLut, dims, prm, targets, rcsCfg);

%% ---- 链路预算审计(解析 vs 仿真测量; snrAnaCalc 局部函数; 仿真测量口径 = 点云 nonCohSNR) ----
% 统一物理口径: SNR = Pr − 10lg(kT·NF) + 10lg(Ts) + 10lg(nChirps) + winStat − 3.0×(实采)
%               + 滤形 − PLL跟踪线损失 − 自扩散抬底(两项均随 R² 增长, 见 snrAnaCalc)。
% 仿真测量口径 = 点云 tgtOut.nonCohSNR —— 【TXRX 非相参积累】图(DDM 选阶后 8TX 对齐 + 8RX
% 功率平均, 64 通道全非相参)的峰单元 / 逐距离门中位噪底 floorPow, 与检测链和点云下游
% 看到的同一数字。功率平均不改均值 SNR(TX/RX 侧同理) → 解析式与单线同形, 仅 winStat
% 随度量口径重标定。对比用目标临时改近改大(用户口径, 防弱线波动); 换回弱目标后
% 单实现散布 ±2~3 dB 属口径内正常(种子已固定, 可复现)。
winStatTbl = struct('HANNING', -2.4, 'SVA', -1.4);  % [口径标定] 窗与统计净项, 随窗类型分档
                               % (2026-10-29 单 bin 口径/同种子/300m 弱目标 5 种子标定:
                               % Hanning −2.40±0.3, SVA −1.41±0.3 —— SVA 比 Hanning 好
                               % ~1.0dB = 收回矩形 ENBW 1.76dB 的大半; 强目标单实现散布
                               % ±1dB 另计)。构成: 双维窗峰损/噪底带宽/行中位统计 + 跨 bin
                               % 跨瓣 + PN/非线性肩部小量(随目标分数 bin 浮动)
winStat95 = winStatTbl.(upper(winType));        % 窗类型旋钮 → 对应档位
cohStat95 = cohStatCalc(dims.chanTotal, dims.aziFftLen, size(antenna_r,2));
                               % [理论统计值, 非凑差] 相参口径统计项: 角度维取峰使噪底
                               % 参考抬升的 max-of-bin 统计 —— 蒙特卡洛计算(cohStatCalc.m),
                               % 随雷达结构自动变(16T16R → nChan=256 重算), 勿写死常数。
                               % 审计残差(理论略优于仿真, ~2dB 口径小量)如实保留, 勿凑差
for k = 1:numel(targets)
    Rm  = targets(k).gate_R * dims.DeltaR;        % 目标距离 (m)
    % DBF 相参增益(实际阵列, cohGainDbf.m, 通道幅度归一化——算法方案口径, 幅度差均
    % 值已由 chTerm 承担): 俯仰角≠0 的 z 行失配自动进入, θel=0 时等于公式 10lg(64)
    cohGainDb = cohGainDbf(A_manifold(3, :).', lambda, targets(k).ele);
    an  = snrAnaCalc(rfInfo, prm, Rm, targets(k).rcsDb, nChirps, size(antenna_r,2), ...
                     isreal(adcRawData), winStat95, [], cohStat95, cohGainDb);
    % 期望单元(检测确认点应落处): 距离 bin + 解模糊后真多普勒 bin(选阶把谱线搬回真
    % bin, 含步进频耦合, 公式与 OMP/DML 种子同源)
    dCk = mod(targets(k).gate_D + prm.stepDir*prm.cr*prm.B/(nChirps-1) ...
              * (2*targets(k).gate_R*c0/(2*prm.B))/c0 * nChirps, nChirps);
    rEk = round(targets(k).gate_R) + 1;           % 期望距离 bin(1 起)
    dEk = mod(round(dCk), nChirps) + 1;           % 期望多普勒 bin(1 起)
    % 点云匹配: 确认点中与期望单元最近者(距离差 + 循环多普勒差/子带宽 加权); 弱目标
    % 亚 bin 偏移可能落在 ±1 bin, 容限 5
    if isempty(tgtOut)
        fprintf('链路预算审计   : 点云为空 —— 无确认目标, 跳过\n');
        break;
    end
    dAbs = abs([tgtOut.dopplerIdx] - dEk);
    dCir = min(dAbs, nChirps - dAbs);             % 多普勒循环距离(子带数单位见加权)
    distK = abs([tgtOut.rangIdx] - rEk) + dCir/dims.subBins;
    [mnD, tK] = min(distK);
    if mnD > 5
        fprintf('链路预算审计   : 目标%d @%.0f m —— 点云无匹配确认点(最近代价 %.1f), 跳过\n', k, Rm, mnD);
        continue;
    end
    meas = tgtOut(tK).nonCohSNR;                  % 仿真测量值 = 点云 nonCohSNR(TXRX 非相参积累口径)
    fprintf('链路预算审计·非相参: 目标%d @%.0f m: Pr %.1f − PSD %.1f + 10lgTs %.1f + 多普勒 %.1f + 通道锥削 %.1f + 窗统计 %.1f − 实采 %.1f + 滤形 %.1f − PLL线损 %.1f − 自扩散抬底 %.1f => 解析 %.1f dB | 仿真值(点云 nonCohSNR @bin R%d/D%d) %.1f dB, 差 %+.1f dB\n', ...
            k, Rm, an.Pr_dBm, an.psd_dBmHz, an.tTerm, an.dTerm, an.chTerm, an.winStat, an.sampP, an.filtG, ...
            an.pllTerm, an.selfFloor, an.SNR_dB, tgtOut(tK).rangIdx-1, tgtOut(tK).dopplerIdx-1, ...
            meas, meas - an.SNR_dB);
    % 相参口径对比: 理论 = 单线解析 + 通道相参(nTX·nRX) + 角度维取峰统计项(cohStat95 =
    % 理论统计值 −7.2dB), 仿真 = 点云 cohSNR(方位峰功率/相参噪底中位)。角度维取投影
    % (max) → 噪底参考抬高, 理论会略优于仿真(残差如实显示, 不凑常数)。
    measC = tgtOut(tK).cohSNR;                    % 仿真测量值 = 点云 cohSNR(TXRX 相参积累口径)
    fprintf('链路预算审计·相参  : 目标%d @%.0f m(ele=%g°): 单线 %.1f + 通道相参 %.1f(阵列DBF, 公式%.1f) + 角度维统计 %.1f(理论) = %.1f dB | 仿真值(点云 cohSNR) %.1f dB, 差 %+.1f dB\n', ...
            k, Rm, targets(k).ele, an.SNR_dB, an.cohTerm, 10*log10(numel(ddmShift)*size(antenna_r,2)), ...
            an.cohStat, an.SNRcoh_dB, measC, measC - an.SNRcoh_dB);
end

%% ---- 后续模块插入位置: 点云/跟踪(吃 tgtOut) ----


%% ==================== 保存图片(PPT 用, 每次运行覆盖) ====================
% 按图窗原尺寸渲染, 再重写 PNG 的 DPI: 插入 PPT 时默认宽度 ≈ outW_cm (16 cm)
saveDir = 'F:\Work3\CH_NEW_95_射频_DDMA';
% 钩子: saveFigDir 改存图目录(批量案例); skipFigs 为不保存的文件名
if exist('saveFigDir', 'var') && ~isempty(saveFigDir), saveDir = saveFigDir; end
if ~exist('skipFigs', 'var'), skipFigs = {}; end
outW_cm = 16;                                 % PPT 插入时的默认宽度 (cm)
res_dpi = 150;                                % 渲染密度(只决定像素数/清晰度)
figList = {figPatTx,            'fig1_AntPat_TX.png'; ...   % (图句柄, 文件名)
           figPatRx,            'fig2_AntPat_RX.png';
           figPatTxCut,         'fig3_AntPatCut_TX.png';
           figPatRxCut,         'fig4_AntPatCut_RX.png';
           figPatTxPh,          'fig5_AntPatPh_TX.png';
           figPatRxPh,          'fig6_AntPatPh_RX.png';
           figPhys,             'fig7_PhysArray.png';
           figVirt,             'fig8_VirtArray.png';
           figPhysC,            'fig9_PhysArrayC.png';
           figVirtC,            'fig10_VirtArrayC.png';
           figWave,             'fig11_FMCW_wave.png';
           % ---- 第一遍: 校准测量(角反场景; 滤波链图=链路身份卡, 放校准组最前) ----
           hFilt,               'fig12_FilterChain.png';
           hAdcCal,             'fig13_ADC_time_cal.png';
           fftCal.figRD,        'fig14_RD_8RX_cal.png';
           fftCal.figSurf,      'fig15_RD_surf_cal.png';
           noiseCal.figNoise,   'fig16_NoiseFloor_cal.png';
           noiseCal.figNoise2,  'fig17_NoiseFloor_med_cal.png';
           calibInfo.figCal,    'fig18_CalibCheck.png';
           rcsCfg.figRcsCal,    'fig19_RCSCalib.png';
           % ---- 第二遍: 目标场景测量 ----
           hAdcScn,             'fig20_ADC_time.png';
           fftCfg.figRD,        'fig21_RD_8RX.png';
           fftCfg.figSurf,      'fig22_RD_surf.png';
           noiseCfg.figNoise,   'fig23_NoiseFloor.png';
           noiseCfg.figNoise2,  'fig24_NoiseFloor_med.png';
           figRxNc,             'fig25_RXNonCoh.png';
           accumCfg.figAccum,   'fig26_JointAccum.png';
           detCfg.figDet,       'fig27_Detect.png';
           detCfg.figDOA,       'fig28_DBF.png';
           ompCfg.figSeed,      'fig29_OMP_process.png'; ...
           ompCfg.figCloud,     'fig30_DBF_OMP.png';
           dmlCfg.figSeed,      'fig31_DML_process.png';
           dmlCfg.figCloud,     'fig32_DBF_OMP_DML.png';
           figTruth,            'fig33_TargetTruth.png';
           figCmpAzi,           'fig34_AziR_Compare.png';
           figRcs,              'fig35_RCS.png'};
tExp = tic;
for kk = 1:size(figList, 1)
    fg   = figList{kk, 1};
    if isempty(fg), continue; end               % 空句柄(如种子点漏检时的过程图)跳过
    if any(strcmp(figList{kk, 2}, skipFigs)), continue; end
    set(fg, 'Name', figList{kk, 2}(1:end-4));   % 窗口名与 PNG 同名(模块内旧号统一改到此)
    fpng = fullfile(saveDir, figList{kk, 2});
    exportgraphics(fg, fpng, 'ContentType', 'image', 'Resolution', res_dpi);
    % 重写 DPI: 物理宽度 = 像素宽/DPI = outW_cm(PPT 按 DPI 定插入尺寸)
    % (PNG 用 XResolution/YResolution; 'Resolution' 是 JPEG/TIFF 参数, 这里会报错)
    img = imread(fpng);
    dpiSet = size(img, 2)/(outW_cm/2.54);
    imwrite(img, fpng, 'XResolution', dpiSet, 'YResolution', dpiSet);
end
fprintf('图片已保存到 %s：%s（原尺寸渲染, PPT默认宽 %g cm, 覆盖上次, 耗时 %.1f s）\n', ...
        saveDir, strjoin(figList(:, 2), ' / '), outW_cm, toc(tExp));
