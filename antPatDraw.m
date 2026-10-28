function [figTx, figRx, figTxCut, figRxCut, figTxPh, figRxPh] = antPatDraw(antPat)
%ANTPATDRAW 画逐天线方向图(图1~6, 检查用; 方向图未接入任何模块)
%   [figTx, figRx, figTxCut, figRxCut, figTxPh, figRxPh] = antPatDraw(antPat)
%   图1/2: 8 TX / 8 RX 增益热图(2×4, dBi) —— 通道间不一致性直接可见
%   图3/4: 8 TX / 8 RX 幅度切面(2×4, 每天线一格): 实线=方位切面(el=0°),
%          虚线=俯仰切面(az=0°); 锚点线 ±70°/±20°/0dBi/副瓣电平
%   图5/6: 8 TX / 8 RX 远场相位切面(2×4): 方向图是复数, 相位含逐天线
%          偏置/指向微偏(线性斜率)/起伏 —— 将来实测导向矢量相位一致性就长这样
%   当前用随机模型模拟通道差异(rng 固定可复现), 换实测图后同版式直接对比。
    screenSize = get(0, 'ScreenSize');
    figPos = [screenSize(3)*0.02, screenSize(4)*0.06, screenSize(3)*0.96, screenSize(4)*0.84];
    az = antPat.azDeg;  el = antPat.elDeg;
    [~, i0] = min(abs(el));                  % el=0° 行
    [~, j0] = min(abs(az));                  % az=0° 列

    % ---- 图1/图2: 增益热图 ----
    figTx = drawMap(antPat.txGainC, az, el, antPat, 'TX', 1, figPos);
    figRx = drawMap(antPat.rxGainC, az, el, antPat, 'RX', 2, figPos);
    % ---- 图3~图6: 切面 ----
    figTxCut = drawCuts(antPat.txGainC, az, el, i0, j0, antPat, 'TX', 3, figPos);
    figRxCut = drawCuts(antPat.rxGainC, az, el, i0, j0, antPat, 'RX', 4, figPos);
    figTxPh = drawPhase(antPat.txGainC, az, el, i0, j0, 'TX', 5, figPos);
    figRxPh = drawPhase(antPat.rxGainC, az, el, i0, j0, 'RX', 6, figPos);
end

function fig = drawMap(gainC, az, el, antPat, tag, figNo, figPos)
%DRAWMAP 图1/2: 每天线增益热图(2×4)
    fig = figure('Position', figPos, 'Name', sprintf('fig%d_AntPat_%s', figNo, tag));
    tiledlayout(2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
    for i = 1:size(gainC, 3)
        nexttile;
        imagesc(az, el, 20*log10(abs(gainC(:, :, i)) + eps));
        axis xy;  clim([-40, antPat.gainPeakDb + 1]);  colormap(parula);
        title(sprintf('%s%d', tag, i)); xlabel('方位 (deg)'); ylabel('俯仰 (deg)');
    end
    colorbar;
    sgtitle(sprintf('%s 方向图幅度热图 (dBi): 主形状一致 + 通道间小幅差异(检查用, 未接入模块)', tag));
end

function fig = drawCuts(gainC, az, el, i0, j0, antPat, tag, figNo, figPos)
%DRAWCUTS 图3/4: 每天线幅度切面(2×4, 方位实线/俯仰虚线)
    fig = figure('Position', figPos, 'Name', sprintf('fig%d_AntPatCut_%s', figNo, tag));
    tiledlayout(2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
    for i = 1:size(gainC, 3)
        nexttile;
        gDb = 20*log10(abs(gainC(:, :, i)) + eps);
        plot(az, gDb(i0, :), '-', 'LineWidth', 1.1); hold on;
        plot(el, gDb(:, j0), '--', 'LineWidth', 1.1);
        if i == 1                                     % 锚点线只画第一格, 避免每格重复
            xline( antPat.azEdgeDeg, 'r:');  xline(-antPat.azEdgeDeg, 'r:');
            xline( antPat.elEdgeDeg, 'm:');  xline(-antPat.elEdgeDeg, 'm:');
            yline(0, 'k:');  yline(antPat.gainPeakDb - antPat.slLevelDb, 'g:');
            legend({'方位切面(el=0°)', '俯仰切面(az=0°)'}, 'Location', 'south', 'FontSize', 7);
        end
        grid on; xlim([-90 90]); ylim([-40, antPat.gainPeakDb + 2]);
        title(sprintf('%s%d (峰 %.1f dBi)', tag, i, max(gDb(:))));
        xlabel('角度 (deg)'); ylabel('增益 (dBi)');
    end
    sgtitle(sprintf('%s 幅度切面(2×4): 实线=方位 / 虚线=俯仰; 红点=±%g°, 品红点=±%g°, 绿点=副瓣 −%gdB', ...
        tag, antPat.azEdgeDeg, antPat.elEdgeDeg, antPat.slLevelDb));
end

function fig = drawPhase(gainC, az, el, i0, j0, tag, figNo, figPos)
%DRAWPHASE 图5/6: 每天线远场相位切面(2×4, 方位实线/俯仰虚线, 展开相位)
    fig = figure('Position', figPos, 'Name', sprintf('fig%d_AntPatPh_%s', figNo, tag));
    tiledlayout(2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
    for i = 1:size(gainC, 3)
        nexttile;
        gC = gainC(:, :, i);                           % 该天线切片(勿漏, 3D 索引会带全页)
        phAz = unwrap(angle(gC(i0, :))) * 180/pi;      % 方位切面相位 (deg)
        phEl = (unwrap(angle(gC(:, j0)))') * 180/pi;   % 俯仰切面相位 (deg)
        plot(az, phAz, '-', 'LineWidth', 1.1); hold on;
        plot(el, phEl, '--', 'LineWidth', 1.1);
        yline(0, 'k:');
        grid on; xlim([-90 90]);
        title(sprintf('%s%d', tag, i));
        xlabel('角度 (deg)'); ylabel('远场相位 (deg)');
    end
    legend({'方位切面(el=0°)', '俯仰切面(az=0°)'}, 'Location', 'southoutside', ...
           'NumColumns', 2, 'FontSize', 7);
    sgtitle(sprintf('%s 远场相位切面(2×4): 逐天线偏置/指向微偏(斜率)/起伏 —— 通道间相位不一致性(检查用)', tag));
end
