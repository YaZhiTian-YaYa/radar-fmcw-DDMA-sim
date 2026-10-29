function figC = compareAziRDraw(targets, dims, stageDbf, stageOmp, stageDml)
%COMPAREAZIRDRAW 方位角-距离四联对比图(图24: 真值/DBF/OMP/DML 同轴同色标)
%   figC = compareAziRDraw(targets, dims, stageDbf, stageOmp, stageDml)
%   2×2 四个子图各取"方位角-距离"切面(即真值图/DBF/OMP/DML 四张点云图的第 2 子图),
%   坐标范围与颜色(高度 z)统一 —— 从左上到右下逐级看: 真值有几个目标 → DBF 融合成了
%   几个点(同单元双目标只见一个) → OMP 分出大间隔弱目标(黑圈) → DML 分出近角对
%   (红圈=DML第一目标(角度已替换), 青圈=DML第二目标)。
    % 轴范围随数据(四图统一取各阶段最大值 + 15% 余量; 全 0 保底防 xlim 非法)
    allRng = [targets.gate_R]*dims.DeltaR;
    if isfield(stageDbf, 'rng'), allRng = [allRng, stageDbf.rng(:)']; end
    if isfield(stageOmp, 'rng'), allRng = [allRng, stageOmp.rng(:)']; end
    if isfield(stageDml, 'rng'), allRng = [allRng, stageDml.rng(:)']; end
    rMax10 = max(ceil(max(allRng) * 1.15 / 10) * 10, 10);          % 距离轴
    allAzi = [targets.azi, stageDbf.azi, stageOmp.azi, stageDml.azi];
    xLimAbs = max(ceil(max(abs(allAzi)) * 1.15 / 10) * 10, 10);    % 方位轴
    allZ = [stageDbf.z, stageOmp.z, stageDml.z];
    allZ = allZ(isfinite(allZ));
    if isempty(allZ), allZ = 0; end
    zAbs = max(abs(allZ));
    cLo = -zAbs;  cHi = zAbs;                       % 色标对称(覆盖各阶段 z 范围)
    zT = [targets.gate_R] * dims.DeltaR .* sind([targets.ele]);

    screenSize = get(0, 'ScreenSize');
    figW = screenSize(3) * 0.98;  figH = screenSize(4) * 0.80;
    figPos = [(screenSize(3)-figW)/2, (screenSize(4)-figH)/2, figW, figH];
    figC = figure('Position', figPos, 'Name', 'fig34_AziR_Compare');
    tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    % ---- (1,1) 真值 ----
    nexttile;
    scatter([targets.azi], [targets.gate_R]*dims.DeltaR, 36, val2rgb(zT), 'filled', 'MarkerEdgeColor', 'k');
    grid on; xlim([-xLimAbs, xLimAbs]); ylim([0, rMax10]); clim([cLo, cHi]); colormap(parula);
    xlabel('方位角 (deg)'); ylabel('距离 (m)');
    title(sprintf('真值(%d 个目标)', numel(targets)));

    % ---- (1,2) DBF(超分辨前) ----
    nexttile;
    scatter(stageDbf.azi, stageDbf.rng, 24, val2rgb(stageDbf.z), 'filled');
    grid on; xlim([-xLimAbs, xLimAbs]); ylim([0, rMax10]); clim([cLo, cHi]);
    xlabel('方位角 (deg)'); ylabel('距离 (m)');
    title(sprintf('DBF 检测点(%d 个)', numel(stageDbf.azi)));

    % ---- (2,1) OMP 后 ----
    nexttile;
    scatter(stageOmp.azi, stageOmp.rng, 24, val2rgb(stageOmp.z), 'filled'); hold on;
    scatter(stageOmp.azi(stageOmp.isO2), stageOmp.rng(stageOmp.isO2), 44, 'ko', 'LineWidth', 1.5);
    grid on; xlim([-xLimAbs, xLimAbs]); ylim([0, rMax10]); clim([cLo, cHi]);
    xlabel('方位角 (deg)'); ylabel('距离 (m)');
    legend({'检测点', 'OMP第二目标(黑圈)'}, 'Location', 'northeast');
    title(sprintf('OMP 后: 比上一级多检出 %d 个点(黑圈)', nnz(stageOmp.isO2)));

    % ---- (2,2) DML 后 ----
    nexttile;
    scatter(stageDml.azi, stageDml.rng, 24, val2rgb(stageDml.z), 'filled'); hold on;
    scatter(stageDml.azi(stageDml.isO2), stageDml.rng(stageDml.isO2), 44, 'ko', 'LineWidth', 1.5);
    scatter(stageDml.azi(stageDml.isD1), stageDml.rng(stageDml.isD1), 44, 'ro', 'LineWidth', 1.5);
    scatter(stageDml.azi(stageDml.isD2), stageDml.rng(stageDml.isD2), 44, 'co', 'LineWidth', 1.5);
    grid on; xlim([-xLimAbs, xLimAbs]); ylim([0, rMax10]); clim([cLo, cHi]); colorbar;
    xlabel('方位角 (deg)'); ylabel('距离 (m)');
    legend({'检测点', 'OMP第二目标(黑圈)', 'DML第一目标(红圈)', 'DML第二目标(青圈)'}, ...
           'Location', 'northeast');
    title(sprintf('DML 后: 比上一级多检出 %d 个点(青圈), 另替换 %d 个母点角度(红圈)', ...
                  nnz(stageDml.isD2), nnz(stageDml.isD1)));
    sgtitle('方位角-距离逐级对比: 真值 → DBF → OMP(大间隔) → DML(近角) —— 同轴同色标');
end
