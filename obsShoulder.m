% obsShoulder: 观察"目标行底噪 近线低→中段下陷→远处抬肩(折叠小副瓣)"现象(实测口径)
% 用法: MATLAB 命令行直接运行  >> obsShoulder
%   场景 = main_sim 默认(当前 100m/50dBsm 测试目标), 全损伤开启。目标行自动检测
%   (换场景距离/RCS 不用改本脚本); 建议强目标(线SNR≥40dB)——肩/副瓣功率∝线功率,
%   弱目标时埋热噪下(实测同理)。
%   看肩随环路带宽移动: 改 rf_sim_core.m 的 loopBW_Hz(200k/400k/800k/1200k)重跑。
addpath(fileparts(mfilename('fullpath')));
here = fileparts(mfilename('fullpath'));
saveFigDir = fullfile(here, 'obsOut');
if ~exist(saveFigDir, 'dir'), mkdir(saveFigDir); end
batchState.idx = 1;  batchState.okRun = true;
run(fullfile(here, 'main_sim.m'));

rxI = 7;
pow2D = mean(abs(fftCfg.data_sva(:, rxI, :)).^2, 3);                % [nR] 行能量(多普勒维平均)
pow2D = pow2D(:);
pow2D(1:10) = 0;                                                    % 排除近场泄漏行
[~, r0] = max(pow2D);                                               % 目标行(自动)
powE = pow2D; powE(max(1,r0-10):min(end,r0+10)) = max(powE);        % 目标±10 保护
[~, rE] = min(powE(20:end-4));  rE = rE + 19;                       % 远离目标的安静空行
rowT = squeeze(20*log10(abs(fftCfg.data_sva(r0, rxI, :)) + eps));
rowE = squeeze(20*log10(abs(fftCfg.data_sva(rE, rxI, :)) + eps));
[~, pkI] = maxk(rowT, 8);                                          % 8 根 DDMA 码线
dMin = zeros(1, 256);
for b = 1:256, dMin(b) = min(abs(mod(b - pkI + 128, 256) - 128)); end
mAll = true(1,256); mAll(mod((pkI-4:pkI+3),256)+1) = false;
dip = dMin>=4 & dMin<16;  shl = dMin>=16 & dMin<32;

% 理论折叠位: 残余相噪 4sin² 峰在 f_hump=c/(4R), 折叠回多普勒带的位置
Rm = (r0-1)*dims.DeltaR;  fHump = 3e8/(4*Rm);
dopSpan = 1/(prm.Tc*prm.nChirps);                                   % 多普勒折叠周期(Hz)
fAl = mod(fHump + dopSpan/2, dopSpan) - dopSpan/2;                  % 折叠回 ±带内
binAl = round(fAl / dopSpan * 256);
fprintf('目标行 %d (%.1f m) | 空行 %d | 4sin²峰 %.0f kHz → 折叠位 %+.0f Hz ≈ 线旁 %d bin\n', ...
    r0, Rm, rE, fHump/1e3, fAl, binAl);

figure('Color','w','Position',[60 60 1000 520]);
plot(0:255, rowT, '-', 'Color', [.65 .78 .95], 'LineWidth', 0.8); hold on;
plot(0:255, medfilt1(rowT, 9), 'b-', 'LineWidth', 2);
plot(0:255, medfilt1(rowE, 9), 'k--', 'LineWidth', 1.5);
plot([0 255], [median(rowE) median(rowE)], 'k:', 'HandleVisibility','off');
ylim([median(rowE)-5, median(rowE)+15]); xlim([0 255]); grid on; box off;
xlabel('多普勒 bin (256)'); ylabel('幅度 (dB)');
lg = legend('目标行原始(码线截顶, 看包络)','目标行 9bin 中值包络','空行 = 纯热噪(9bin中值)', ...
    'Location','south','Orientation','horizontal');
lg.FontSize = 8;
title({sprintf('目标行多普勒剖面(%.0fm): 线旁主瓣裙 → 中段下陷(谷 %+.1fdB) → 远处抬肩(%+.1fdB)', ...
    Rm, median(rowT(dip))-median(rowE), median(rowT(shl))-median(rowE)); ...
    sprintf('空行热噪 %.1fdB | 全行中位 %+.1fdB ≈ 肩位 = 实测"目标距离处不高" | 4sin²峰折叠位: 线旁 %d bin', ...
    median(rowE), median(rowT(mAll))-median(rowE), binAl)});
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
fprintf('"小副瓣对"= 4sin² 峰(c/4R=%.0fkHz)折叠像, 位置随距离/loopBW 移动; 高度 ∝ 线功率×L_vco(峰区)。\n', fHump/1e3);
fprintf('VCO 支路标定: 实测副瓣高度(已知线强度)反推 pnProfileVco 30k/100k 锚点(见报告 2026-11-08 追记)。\n');
fprintf('图已存 %s\n', fullfile(saveFigDir, 'obsShoulder.png'));
