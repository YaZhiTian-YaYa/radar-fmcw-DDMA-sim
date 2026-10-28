% quantProbe.m —— int16 量化(2 字节整数 ADC 数据)对 SNR 的影响
% 同一场景(案例1), 成对 A/B: 每对内 rng 种子相同 → 噪声实现完全一致,
% 唯一差异 = 是否整数量化。两组增益:
%   G_rx=66dB(当前标定, 底噪 rms≈57 counts) 与 G_rx=36dB(压增益, 底噪≈1.8 counts,
%   故意让底噪逼近 1~2 LSB, 看量化影响何时显现)。
% 指标: 量化误差 rms(counts)、逐距离门中位噪底(dB)、角反(近距强目标)峰 SNR。
close all;
S = load(fullfile('OMP和DML理论仿真', '案例1', 'rawdata.mat'));
wR = 0.5-0.5*cos(2*pi*(0:1023)'/(1023));   % 复刻 win2dfft 的 Hanning(只取指标, 不画图)
wC = 0.5-0.5*cos(2*pi*(0:511)/(511));

cfgs = { 66, false, 'G_rx=66dB 理想(不量化)';
         66, true,  'G_rx=66dB int16 量化';
         36, false, 'G_rx=36dB 理想(不量化)';
         36, true,  'G_rx=36dB int16 量化' };
res = zeros(size(cfgs,1), 3);   % [底噪中位dB, 角反峰SNR dB, 量化误差rms counts]
adcSave = cell(size(cfgs,1), 1);
for q = 1:size(cfgs, 1)
    sw = struct('G_rx_dB', cfgs{q,1}, 'enableInt16Quant', cfgs{q,2});
    rng(42);                                    % 同对同种子 → 噪声实现一致, 差异纯量化
    cmd = ['adc = rf_sim_core(1024, 512, S.TX_A, S.RX_A, S.antenna_t, ' ...
           'S.antenna_r, S.ddmShift, S.prm, S.targets, sw);'];
    evalc(cmd);  close all;
    adcSave{q} = double(adc);
    RD = fft(fft(adcSave{q}(:,:,1).*wR, [], 1).*wC, [], 2);   % RX1 的 RD(加窗口径同链路)
    P  = abs(RD).^2;
    flr = median(P(40:480, :), 'all');           % 底噪(远段距离门 × 全多普勒中位)
    pk  = max(P(4:8, :), [], 'all');             % 角反(近距 bin5±)峰
    res(q,:) = [10*log10(flr), 10*log10(pk/flr), 0];
end
% 量化误差(同增益成对)
res(2,3) = sqrt(mean((adcSave{2}(:)-adcSave{1}(:)).^2));
res(4,3) = sqrt(mean((adcSave{4}(:)-adcSave{3}(:)).^2));

fprintf('\n================ int16 量化影响(同噪声实现 A/B) ================\n');
fprintf('%-26s %12s %14s %16s\n', '配置', '底噪中位(dB)', '角反峰SNR(dB)', '量化误差rms(cnts)');
for q = 1:size(cfgs,1)
    fprintf('%-26s %12.2f %14.2f %16.4f\n', cfgs{q,3}, res(q,1), res(q,2), res(q,3));
end
fprintf('ΔSNR(量化 vs 理想): G_rx=66dB: %+.4f dB | G_rx=36dB: %+.4f dB\n', ...
        res(2,2)-res(1,2), res(4,2)-res(3,2));
fprintf('底噪抬升(量化 vs 理想): G_rx=66dB: %+.4f dB | G_rx=36dB: %+.4f dB\n', ...
        res(2,1)-res(1,1), res(4,1)-res(3,1));
