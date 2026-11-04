% 第一章示意图生成: 总流程 / 五分类路由 / 分离度路由(MATLAB 盒图; 局部函数在文件尾)
figDir = 'C:\Users\NJ\AppData\Local\Temp\probe95\docfigs';
if ~exist(figDir, 'dir'), mkdir(figDir); end
bk = [0.93 0.93 0.96];  yl = [1 0.9 0.35];  gn = [0.75 0.88 0.78];  rd = [0.95 0.8 0.8];

% ---- 图 1-1 总流程(两遍) ----
f = figure('Visible', 'off', 'Position', [40 40 1000 430], 'Color', 'w'); hold on; axis off;
text(0.25, 0.97, '校准遍(暗箱, 只放角反)', 'FontSize', 11, 'FontWeight', 'bold');
drawBox(0.015, 0.72, 0.15, 0.13, {'角反 3.5m', '-5 dBsm(已知)'}, yl, 9);
drawBox(0.225, 0.72, 0.16, 0.13, {'0°校准表 cCal', '(通道相位掰齐)'}, bk, 9);
drawBox(0.445, 0.72, 0.155, 0.13, {'噪底曲线+滤形', '重建(档位锚定)'}, bk, 9);
drawBox(0.66, 0.72, 0.115, 0.13, {'K0', '(alpha域锚定)'}, bk, 9);
drawBox(0.835, 0.72, 0.115, 0.13, {'K0Coh', '(相参域)'}, bk, 9);
arrow(0.165, 0.785, 0.225, 0.785); arrow(0.385, 0.785, 0.445, 0.785);
arrow(0.60, 0.785, 0.66, 0.785); arrow(0.775, 0.785, 0.835, 0.785);
plot([0.01 0.99], [0.66 0.66], 'k--');
text(0.25, 0.615, '目标遍(冻结标定量, 换目标场景)', 'FontSize', 11, 'FontWeight', 'bold');
drawBox(0.02, 0.44, 0.16, 0.13, {'点云点', '(角色标记)'}, yl, 9);
drawBox(0.24, 0.53, 0.20, 0.17, {'独立点支路:', '相参图读数', '峰邻域+峰值插值', '+俯仰补+PLL补'}, gn, 9);
drawBox(0.24, 0.32, 0.20, 0.13, {'配对支路: alpha域联合权重', '(OMP/DML 四类成员)'}, gn, 9);
drawBox(0.52, 0.44, 0.16, 0.13, {'sigma = SNR+40lgR', '+ K(支路各配)'}, yl, 9);
drawBox(0.75, 0.44, 0.22, 0.13, {'写回点云 rcsDb', '+ 置信标志 rcsConf'}, bk, 9);
arrow(0.18, 0.505, 0.24, 0.615); arrow(0.18, 0.505, 0.24, 0.385);
arrow(0.44, 0.615, 0.52, 0.51); arrow(0.44, 0.385, 0.52, 0.49);
arrow(0.68, 0.505, 0.75, 0.505);
exportgraphics(f, [figDir '\fig_s1_flow.png'], 'Resolution', 170); close(f);

% ---- 图 1-2 五分类路由 ----
f = figure('Visible', 'off', 'Position', [40 40 1000 400], 'Color', 'w'); hold on; axis off;
drawBox(0.02, 0.40, 0.15, 0.16, {'点云点', 'splitTag 判角色'}, yl, 10);
labs = {'独立点 (0)', 'OMP第一目标 (1)', 'OMP第二目标 (2)', 'DML第一目标 (3)', 'DML第二目标 (4)'};
for k = 1:5
    drawBox(0.26, 0.78 - (k - 1) * 0.17, 0.20, 0.11, labs{k}, bk, 9);
    arrow(0.17, 0.48, 0.26, 0.835 - (k - 1) * 0.17);
end
drawBox(0.54, 0.62, 0.24, 0.16, {'相参图读数支路', '(仅独立点)', 'K0Coh 尺'}, gn, 10);
drawBox(0.54, 0.20, 0.24, 0.24, {'alpha 域联合权重支路', '(四类配对成员)', 'K0 尺'}, gn, 10);
arrow(0.46, 0.835, 0.54, 0.68, '独立点');
arrow(0.46, 0.30, 0.54, 0.32, '配对成员');
drawBox(0.84, 0.40, 0.14, 0.16, {'sigma 写回', '点云'}, yl, 10);
arrow(0.78, 0.70, 0.84, 0.52); arrow(0.78, 0.32, 0.84, 0.46);
exportgraphics(f, [figDir '\fig_s2_route.png'], 'Resolution', 170); close(f);

% ---- 图 1-3 分离度路由 ----
f = figure('Visible', 'off', 'Position', [40 40 1000 330], 'Color', 'w'); hold on; axis off;
plot([0.05 0.95], [0.52 0.52], 'k-', 'LineWidth', 2);
text(0.05, 0.60, '双目标角间隔 ->', 'FontSize', 10);
xx = [0.16, 0.52, 0.78];
for k = 1:3
    plot(xx(k), 0.52, 'ks', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
end
drawBox(0.07, 0.68, 0.30, 0.17, {'OMP 区 (约 >=2°)', '残差过 12 dB 门限', '-> OMP点1 / OMP点2'}, gn, 9);
drawBox(0.44, 0.68, 0.22, 0.17, {'DML 区 (约 0.6°~1.2°)', '双原子联合拟合', '残差改善 > 3 dB', '(还需弱成员电平)', '-> DML点1 / DML点2'}, yl, 8);
drawBox(0.69, 0.68, 0.24, 0.17, {'不可分区 (约 <0.5°)', '联合拟合也无法分离', '-> 输出单点(融合值)'}, rd, 9);
arrow(0.16, 0.68, 0.16, 0.56); arrow(0.52, 0.68, 0.52, 0.56); arrow(0.78, 0.68, 0.78, 0.56);
text(0.16, 0.44, '第一原子不吸收第二目标', 'FontSize', 8, 'Color', [0.3 0.3 0.3], 'HorizontalAlignment', 'center');
text(0.52, 0.42, {'第一原子吸收双目标能量,', 'OMP 残差塌掉'}, 'FontSize', 8, 'Color', [0.3 0.3 0.3], 'HorizontalAlignment', 'center');
text(0.80, 0.44, '条件数极限', 'FontSize', 8, 'Color', [0.3 0.3 0.3]);
exportgraphics(f, [figDir '\fig_s3_sep.png'], 'Resolution', 170); close(f);
fprintf('SCHEMATICS DONE\n');

function drawBox(x, y, w, h, txt, fc, fs)
    rectangle('Position', [x, y, w, h], 'FaceColor', fc, 'EdgeColor', [0.3 0.3 0.3], 'LineWidth', 1);
    text(x + w/2, y + h/2, txt, 'HorizontalAlignment', 'center', 'FontSize', fs, 'Color', [0.1 0.1 0.1]);
end
function arrow(x1, y1, x2, y2, txt)
    plot([x1 x2], [y1 y2], '-', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.2);
    th = atan2(y2 - y1, x2 - x1);
    patch(x2 + [0, -0.014*cos(th) - 0.006*sin(th), -0.014*cos(th) + 0.006*sin(th)], ...
          y2 + [0, -0.014*sin(th) + 0.006*cos(th), -0.014*sin(th) - 0.006*cos(th)], ...
          [0.3 0.3 0.3], 'EdgeColor', 'none');
    if nargin > 4 && ~isempty(txt)
        text((x1 + x2)/2, (y1 + y2)/2 + 0.015, txt, 'HorizontalAlignment', 'center', 'FontSize', 8, 'Color', [0.25 0.25 0.25]);
    end
end
