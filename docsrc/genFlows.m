% 三条支路的 RCS 计算流程图(独立点/OMP双点/DML双点) —— 垂直流程 + 判决分支 + 关键步骤注"为什么"
figDir = 'C:\Users\NJ\AppData\Local\Temp\probe95\docfigs';
if ~exist(figDir, 'dir'), mkdir(figDir); end
bk = [0.93 0.93 0.96];  yl = [1 0.9 0.35];  gn = [0.75 0.88 0.78];  rd = [0.95 0.8 0.8];

% ============ 图 F-A: 独立点(相参图支路) ============
f = figure('Visible', 'off', 'Position', [40 40 900 1120], 'Color', 'w'); hold on; axis off;
text(0.5, 0.985, '独立点 RCS 计算（相参图支路, 尺=K0Coh）', 'FontSize', 13, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
xb = 0.30; wb = 0.40; hh = 0.058; y0 = 0.905; dy = 0.082;
steps = {
 '点云点（角色=独立点）+ 定位 (距离bin, 多普勒bin)', yl;
 'TXRX 相参积累图取方位峰：峰 bin ±1 邻域取最大', bk;   % ↓为什么
 '角维三点抛物线定峰（三点非凹→回退原格读数）', bk;
 '距离剖面 hannPk3 + 8 根码线多普勒恢复（精确叶损）', bk;
 '俯仰损失补回：θ̂ 取同簇中位, 夹 ±15°（产品口径）', bk;
 'PLL 线损补回 4.343σφ²（σφ: 芯片PLL参数 + 自估距离）', bk;
 'SNR = 相参功率读数 − 噪底参考', bk;
 'σ = SNR + 40lgR + K0Coh（角反相参域同款读法标定）', gn;
 '写回点云 rcsDb（+ rcsConf 置信标志）', yl};
for k = 1:size(steps,1)
    drawBox(xb, y0-(k-1)*dy, wb, hh, steps{k,1}, steps{k,2}, 9.5);
    if k < size(steps,1), arrow(xb+wb/2, y0-(k-1)*dy, xb+wb/2, y0-(k-1)*dy-dy+hh); end
end
notes = {
 0.74, 0.905-1*dy+0.012, '为什么±1: 角栅格1.4°/bin, 目标卡bin间有±0.7dB纹波';
 0.74, 0.905-3*dy+0.012, '为什么hannPk3: 普通抛物线半bin处低估1.18dB(实测−1.58)';
 0.74, 0.905-4*dy+0.012, '为什么簇中位: 单点θ̂有±5°野值(实测多补+3.5dB)';
 0.74, 0.905-5*dy+0.012, '为什么显式补: 角反3.5m几乎无损, 常数吸不掉(损失∝R)';
 0.74, 0.905-7*dy+0.012, '同口径相消: 标定(角反)与估计(目标)走同一条支路'};
for q = 1:size(notes, 1)
    text(notes{q,1}, notes{q,2}, notes{q,3}, 'FontSize', 8, 'Color', [0.35 0.35 0.35]);
end
exportgraphics(f, [figDir '\fig_f_a.png'], 'Resolution', 165); close(f);

% ============ 图 F-B: OMP 双点(α域, 分离良好) ============
f = figure('Visible', 'off', 'Position', [40 40 900 1180], 'Color', 'w'); hold on; axis off;
text(0.5, 0.985, 'OMP 双点 RCS 计算（α域联合权重, 尺=K0, 角隔约≥2°）', 'FontSize', 13, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
xb = 0.30; wb = 0.42; hh = 0.055; y0 = 0.905; dy = 0.075;
steps = {
 '检测点（簇）→ 36 通道方位行快拍 x', yl;
 '原子1 = 导向表插值 @ DBF 方位角', bk;
 'LS 扣除第一原子 → 残差匹配曲线', bk;
 '判决: 残差峰 − 中位 > 12 dB ?', yl;
 '第二角度 = 残差峰抛物线亚格定位', bk;
 '双原子联合 LS（κ≤4 pinv, 共线限增益）', bk;
 '母点取 W1, 追加点取 W2（追加进点云, tag=OMP第二）', bk;
 '逐点叶损恢复（距离 hannPk3 + 8 码线多普勒均值）', bk;
 'σ_i = (W_i − 噪底参考) + 40lgR + K0, 两成员各写回', gn};
yB = y0; yPg = {};
for k = 1:size(steps,1)
    if k == 4
        drawBox(xb, yB, wb, hh+0.012, steps{k,1}, steps{k,2}, 9.5);
    else
        drawBox(xb, yB, wb, hh, steps{k,1}, steps{k,2}, 9.5);
    end
    yPg{end+1} = yB; %#ok<SAGROW>
    yB = yB - dy;
end
for k = 1:size(steps,1)-1
    arrow(xb+wb/2, yPg{k}, xb+wb/2, yPg{k+1}+hh);
end
% "否"分支 → 独立点
arrow(xb, yPg{4}+hh/2, xb-0.14, yPg{4}+hh/2, '否');
drawBox(0.02, yPg{4}-0.02, 0.14, 0.06, {'该点按独立点', '走相参图支路'}, rd, 8.5);
text(xb+wb+0.02, yPg{4}+hh/2, '是', 'FontSize', 9, 'Color', [0.2 0.2 0.2]);
notes = {
 0.76, yPg{3}+0.01, '12dB门: 分离良好时扣掉强目标后弱目标仍显著';
 0.76, yPg{6}+0.01, '为什么联合而非单原子: 互投影对消, 两成员对称';
 0.76, yPg{8}+0.01, '旧口径(插值距离+统计常数)有±2dB幻影, 已弃用'};
for q = 1:size(notes, 1)
    text(notes{q,1}, notes{q,2}, notes{q,3}, 'FontSize', 8, 'Color', [0.35 0.35 0.35]);
end
exportgraphics(f, [figDir '\fig_f_b.png'], 'Resolution', 165); close(f);

% ============ 图 F-C: DML 双点(α域, 亚波束) ============
f = figure('Visible', 'off', 'Position', [40 40 900 1240], 'Color', 'w'); hold on; axis off;
text(0.5, 0.985, 'DML 双点 RCS 计算（α域联合权重, 尺=K0, 角隔约 0.7°~1.2°）', 'FontSize', 13, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
xb = 0.30; wb = 0.42; hh = 0.052; y0 = 0.905; dy = 0.071;
steps = {
 'OMP 未出第二目标的母点 → 36 通道快拍 x', yl;
 '单原子基线: 网格 J1 + 亚bin精化(抛物线角重投影)', bk;
 '残差峰候选(主峰 ±8°内, 多起点)', bk;
 '交替坐标连续搜索(0.05° 步进 ×3 轮) → 双角最优', bk;
 '四道门: 残差改善>3dB 且 净改善imprS>1dB', yl;
 '      且 间隔≥0.7° 且 条件数 κ_G≤16', yl;
 '母点方位←ang1, 追加点←ang2(tag=DML第二)', bk;
 '双原子联合 LS(κ≤4 pinv) → W1 / W2', bk;
 '逐点叶损恢复 → σ_i = (W_i−噪底)+40lgR+K0', gn};
yB = y0; yPg = {};
for k = 1:size(steps,1)
    if k == 5 || k == 6
        drawBox(xb, yB, wb, hh, steps{k,1}, steps{k,2}, 9.5);
        if k == 5, yGate = yB; end
    else
        drawBox(xb, yB, wb, hh, steps{k,1}, steps{k,2}, 9.5);
    end
    yPg{end+1} = yB; %#ok<SAGROW>
    yB = yB - dy;
end
% 5/6 两行画成一个大判决框边框
rectangle('Position', [xb-0.008, yPg{6}-0.006, wb+0.016, 2*hh+0.02], 'EdgeColor', [0.6 0.2 0.2], 'LineWidth', 1.4);
for k = [1 2 3 4]
    arrow(xb+wb/2, yPg{k}, xb+wb/2, yPg{k+1}+hh);
end
arrow(xb+wb/2, yPg{6}, xb+wb/2, yPg{7}+hh);
for k = 7:8
    arrow(xb+wb/2, yPg{k}, xb+wb/2, yPg{k+1}+hh);
end
arrow(xb, yPg{5}+hh, xb-0.14, yPg{5}+hh, '否→');
drawBox(0.02, yPg{5}-0.05, 0.15, 0.075, {'按设计输出融合单点', '(独立支路读数', '≈两目标相干和)'}, rd, 8.5);
text(xb+wb+0.02, yPg{6}, '是', 'FontSize', 9, 'Color', [0.2 0.2 0.2]);
notes = {
 0.76, yPg{2}+0.005, '为什么亚bin基线: J1只在1°格点取峰会把"亚bin主瓣恢复"';
 0.76, yPg{2}-0.012, '偷算成改善量 → imprS剥掉该偏置(假0.3~0.5/真等强3.0)';
 0.76, yPg{5}+0.005, '间隔/条件数门: 拦"两原子塌到强目标一侧"的虚假分离;';
 0.76, yPg{5}-0.012, '强弱悬殊(差≥8dB)差分低于噪声×自由度 → 单帧不可分(已对账)'};
for q = 1:size(notes, 1)
    text(notes{q,1}, notes{q,2}, notes{q,3}, 'FontSize', 8, 'Color', [0.35 0.35 0.35]);
end
exportgraphics(f, [figDir '\fig_f_c.png'], 'Resolution', 165); close(f);
fprintf('FLOWS DONE\n');

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
        text((x1 + x2)/2, (y1 + y2)/2 + 0.012, txt, 'HorizontalAlignment', 'center', 'FontSize', 8, 'Color', [0.25 0.25 0.25]);
    end
end
