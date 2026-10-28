%% ======================================================================
%  毫米波汽车雷达信号处理流水线(主入口)
%  ---------------------------------------------------------------------
%  数据流向:
%    rawData [距离采样点 × 接收通道 × chirp]  (来自 AdcAndParam.mat)
%      |-> 1. fft2D_fun_pip      距离/多普勒二维FFT + SVA副瓣抑制
%      |-> 2. estNoiseFloor_pipe 逐距离门噪底估计, 标记强干扰门
%      |-> 3. jointAccum_pipe    DDM解模糊 + 非相参/相参积累 + 方位维FFT
%      |-> 4. targetDetect_pipe  方位谱DACFAR粗检 + 峰值状态精判 -> 目标列表
%      |-> 5. doaEstimate_pipe   方位DBF + 子阵俯仰DBF测角
%      `-> 6. point_disp_pipe    点云可视化
%  AdcAndParam.mat 中还含 cfarCfg(逐距离门CFAR/峰判配置)与
%  AntCaliData(天线导向矢量标定数据), 分别供检测与测角使用
%  ======================================================================
clear all; close all; clc;
addpath("./pipeline");

load('./AdcAndParam.mat');               % 载入 rawData / cfarCfg / AntCaliData

%% 1. 二维FFT + SVA
%  距离维(快时间)与多普勒维(慢时间chirp)分别FFT, 再做SVA自变迹抑制副瓣
[fftCfg] = fft2D_fun_pip(rawData);

%% 2. 噪底估计 与 3. 积累
[fftCfg,noiseCfg] = estNoiseFloor_pipe(fftCfg);   % 逐距离门估计噪底功率
[accumCfg] = jointAccum_pipe(fftCfg);             % DDM解模糊选阶 + 非相参/相参积累 + 方位FFT
%% 4. cfar 检测
%  方位谱DACFAR粗检, 再结合相参/非相参峰值状态精判, 输出目标列表
[accumCfg, detCfg, tgtOut] = targetDetect_pipe(noiseCfg,accumCfg,cfarCfg);
%% 5. dbf 测角
%  对每个目标做方位DBF与子阵俯仰DBF, 测得方位/俯仰角并写回目标列表
[tgtOut, detCfg] = doaEstimate_pipe(detCfg, tgtOut);
%% 6. disp 点云显示
point_disp_pipe(tgtOut);
