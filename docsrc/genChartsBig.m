% 大批量图表: 读 rerunD_state.mat(66 场景) → results.csv + 三批统计图
S = load('C:\Users\NJ\AppData\Local\Temp\probe95\rerunD_state.mat');
R = S.RES;  nmAll = {S.scenes.name};
figDir = 'C:\Users\NJ\AppData\Local\Temp\probe95\docfigs';
if ~exist(figDir, 'dir'), mkdir(figDir); end
tagN = {'独立点', 'OMP第一目标', 'OMP第二目标', 'DML第一目标', 'DML第二目标'};
samePt = zeros(size(R, 1), 1);
nmRow = nmAll(R(:, 1)).';                          % 逐行所属场景名(194×1)
for q = 1:size(R, 1)
    if R(q, 14) == 1 && ~isnan(R(q, 15))
        selQ = strcmp(nmRow, nmRow{q});
        sib = find(selQ & R(:, 2) ~= R(q, 2) & R(:, 14) == 1 & R(:, 15) == R(q, 15));
        if ~isempty(sib), samePt(q) = 1; end
    end
end
csvF = fopen([figDir '\results.csv'], 'w');
fprintf(csvF, 'case,tgt,truth,fit,dev,splitTag,conf,dom,R,azi,ele,seed,isPair,samePt\n');
for q = 1:size(R, 1)
    if isnan(R(q, 1)), continue; end
    ip = ~isnan(R(q, 14)) && R(q, 14) == 1;
    fprintf(csvF, '%s,%d,%.2f,%.2f,%.2f,%d,%d,%d,%.1f,%.1f,%.1f,%d,%d,%d\n', ...
        nmRow{q}, R(q, 2), R(q, 3), R(q, 4), R(q, 5), ...
        R(q, 6), R(q, 7), R(q, 8), R(q, 9), R(q, 10), R(q, 11), R(q, 12), ip, samePt(q));
end
fclose(csvF);
disp(table(readtable([figDir '\results.csv'])));
fprintf('BIG CSV DONE\n');

col = [0.20 0.40 0.75; 0.85 0.15 0.75; 0.10 0.10 0.10; 0.85 0.20 0.20; 0.10 0.70 0.80];
elig = ~isnan(R(:, 4)) & (R(:, 14) ~= 1 | (R(:, 6) >= 1 & R(:, 6) <= 4 & samePt == 0));
campF = {'I', 'fig_I_stat.png'; 'O', 'fig_O_stat.png'; 'D', 'fig_D_stat.png'};
for ci = 1:3
    sel = startsWith(nmRow, campF{ci, 1});  sel = sel(:);
    RESc = R(sel & elig, :);
    f = figure('Visible', 'off', 'Position', [40 40 980 430], 'Color', 'w');
    subplot(1, 2, 1); hold on; box on;
    hh = [];  lg = {};
    for tg = 0:4
        s = RESc(:, 6) == tg;
        if any(s)
            hh(end+1) = plot(RESc(s, 3), RESc(s, 4), 'o', 'MarkerEdgeColor', col(tg+1, :), ...
                'MarkerFaceColor', 0.35*col(tg+1, :) + 0.65, 'MarkerSize', 6); %#ok<AGROW>
            lg{end+1} = tagN{tg+1}; %#ok<AGROW>
        end
    end
    allv = [RESc(:, 3); RESc(:, 4)];  allv = allv(~isnan(allv));
    lim = [floor(min(allv)) - 1, ceil(max(allv)) + 1];
    plot(lim, lim, 'k--', 'LineWidth', 1);
    axis equal;  xlim(lim);  ylim(lim);  grid on;
    xlabel('真值 \sigma (dBsm)');  ylabel('拟合 \sigma (dBsm)');
    title('拟合 vs 真值 (虚线 y = x)', 'FontSize', 10);
    if ~isempty(hh), legend(hh, lg, 'Location', 'northwest', 'FontSize', 8); end
    subplot(1, 2, 2); box on; hold on;
    dv = R(sel & elig, 5);  dv = dv(~isnan(dv));
    if ~isempty(dv)
        histogram(dv, 'BinWidth', 0.5, 'FaceColor', [0.55 0.65 0.85], 'EdgeColor', [0.3 0.4 0.6]);
        plot(mean(dv) + [0 0], ylim, 'r-', 'LineWidth', 1.5);
        grid on;  xlabel('\Delta\sigma = 拟合-真值 (dB)');  ylabel('点数');
        ttl = sprintf('n=%d   \\mu=%+.2f dB   std=%.2f dB   [%+.2f, %+.2f]', ...
            numel(dv), mean(dv), std(dv), min(dv), max(dv));
    else
        ttl = '无匹配点';
    end
    nu = sum(R(sel & (R(:, 14) == 1) & ~isnan(R(:, 4)) & ~elig));
    nmis = sum(sel & isnan(R(:, 4)));
    if nu > 0, ttl = sprintf('%s   (未分离 %d)', ttl, nu); end
    if nmis > 0, ttl = sprintf('%s   (未检出 %d)', ttl, nmis); end
    title(ttl, 'FontSize', 9.5);
    exportgraphics(f, [figDir '\' campF{ci, 2}], 'Resolution', 160);  close(f);
end
fprintf('BIG FIGS DONE\n');
