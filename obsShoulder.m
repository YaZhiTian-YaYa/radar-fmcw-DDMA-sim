% obsShoulder: 观察"目标行底噪 近线低→中段下陷→远处抬肩"现象(实测口径)
% 用法: MATLAB 命令行直接运行  >> obsShoulder
%   场景 = main_sim 默认(角反 + 50dBsm@300m 测试目标), 全损伤开启(hopErr=500Hz/pnTrack=0)。
%   为什么用强目标: 肩功率 ∝ 目标线功率, 50dBsm(线SNR≈45dB)肩才高出热噪几个 dB 可见;
%   产品 RCS(≤15dBsm)时肩埋在热噪下 —— 实测能看到肩的也是强目标/近距离目标。
%   看肩随环路带宽移动: 改 rf_sim_core.m 的 loopBW_Hz(200k/400k/800k/1200k)重跑本脚本。
addpath(fileparts(mfilename('fullpath')));
here = fileparts(mfilename('fullpath'));
saveFigDir = fullfile(here, 'obsOut');
if ~exist(saveFigDir, 'dir'), mkdir(saveFigDir); end
batchState.idx = 1;  batchState.okRun = true;
run(fullfile(here, 'main_sim.m'));

r0 = 429;  rE = 350;  rxI = 7;                    % 目标行(300m) / 空行 / RX7
rowT = squeeze(20*log10(abs(fftCfg.data_sva(r0, rxI, :)) + eps));
rowE = squeeze(20*log10(abs(fftCfg.data_sva(rE, rxI, :)) + eps));
[~, pkI] = maxk(rowT, 8);                          % 8 根 DDMA 码线
dMin = zeros(1, 256);
for b = 1:256, dMin(b) = min(abs(mod(b - pkI + 128, 256) - 128)); end
mAll = true(1,256); mAll(mod((pkI-4:pkI+3),256)+1) = false;   % 排线±4(主瓣裙)

figure('Color','w','Position',[60 60 1000 520]);
plot(0:255, rowT, '-', 'Color', [.65 .78 .95], 'LineWidth', 0.8); hold on;
plot(0:255, medfilt1(rowT, 9), 'b-', 'LineWidth', 2);
plot(0:255, medfilt1(rowE, 9), 'k--', 'LineWidth', 1.5);
% 谷底与肩位标注
dip = dMin>=4 & dMin<16;  shl = dMin>=16 & dMin<32;
plot([0 255], [median(rowE) median(rowE)], 'k:', 'HandleVisibility','off');
yL = [median(rowE)-5, median(rowE)+15];
ylim(yL); xlim([0 255]); grid on; box off;
xlabel('多普勒 bin (256)'); ylabel('幅度 (dB)');
lg = legend('目标行原始(码线被图窗截顶, 看包络)','目标行 9bin 中值包络','空行 = 纯热噪(9bin中值)', ...
    'Location','south','Orientation','horizontal');
lg.FontSize = 8;
title({sprintf('目标行多普勒剖面: 线旁主瓣裙 → 中段下陷(谷 %+.1fdB) → 远处抬肩(%+.1fdB)', ...
    median(rowT(dip))-median(rowE), median(rowT(shl))-median(rowE)); ...
    sprintf('空行热噪 %.1fdB | 全行中位(排线±4) %+.1fdB ≈ 肩位 = 实测"目标距离处不高"', ...
    median(rowE), median(rowT(mAll))-median(rowE))});
saveas(gcf, fullfile(saveFigDir, 'obsShoulder.png'));

fprintf('\n===== 距线分桶(中位 − 空行中位) =====\n');
for q = [2 4 8 16 32 64]
    mm = dMin>=q & dMin<q*2;
    if any(mm)
        fprintf('距线 %3d-%3d: 目标行中位 %7.2f | 差 %+.2f dB\n', q, q*2-1, median(rowT(mm)), median(rowT(mm))-median(rowE));
    end
end
fprintf('全行中位(排线±4): %+.2f dB | 肩位(16-31): %+.2f dB | 谷底(4-15): %+.2f dB\n', ...
    median(rowT(mAll))-median(rowE), median(rowT(shl))-median(rowE), median(rowT(dip))-median(rowE));
fprintf('现象判读: 谷底低于肩 = "先下降再抬高"; 全行中位≈肩位 = 实测"目标距离处不高"。\n');
fprintf('对比环路带宽: rf_sim_core.m 改 loopBW_Hz = 200e3/400e3/800e3/1200e3 重跑, 看 16-31/32-63 桶移动。\n');
fprintf('图已存 %s\n', fullfile(saveFigDir, 'obsShoulder.png'));
