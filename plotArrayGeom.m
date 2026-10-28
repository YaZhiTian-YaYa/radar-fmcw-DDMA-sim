function [figPhys, figVirt] = plotArrayGeom(antenna_t, antenna_r, A_manifold, FFT_d0)
%PLOTARRAYGEOM 阵列几何绘图(图7 物理阵列 / 图8 虚拟阵列)
%   [figPhys, figVirt] = plotArrayGeom(antenna_t, antenna_r, A_manifold, FFT_d0)
%   输入:
%     antenna_t/r [3×n] 实际坐标(m, 已乘 FFT_d0)——建议传【居中前的设计原始值】
%     A_manifold  [3×nTX*nRX] 虚拟阵元坐标(TX+RX 矢量和, TX 内循环/RX 外循环)
%     FFT_d0      阵列间距(m)
%   说明: 内部除回 FFT_d0 取整, 用【整数栅格】(×FFT_d0 之前)绘图, 便于对照设计;
%         网格按单位 1 画线(每 5 格标注一个数); 只用第 1 行(x)与第 3 行(z)。
    tInt = antenna_t / FFT_d0;              % 物理阵元栅格坐标(设计值为整数; 不做 round
    rInt = antenna_r / FFT_d0;              %   —— 居中平移可能是半整数, 取整会伪造 0.5 格
    vInt = A_manifold / FFT_d0;             %   位置偏移, 平移图必须保持相对位置严格不变)
    nTX = size(tInt, 2);
    nRX = size(rInt, 2);

    screenSize = get(0, 'ScreenSize');
    figPos = [screenSize(3)*0.08, screenSize(4)*0.05, screenSize(3)*0.84, screenSize(4)*0.85];

    % ---- 天线封装尺寸(单位: FFT_d0 栅格; 图7 上以天线位置为中心画矩形, ----
    % ---- 用于观察天线之间有没有位置上的干涉) ----
    antSizeX = 3;                 % [旋钮] 天线 x 方向尺寸(格)
    antSizeZ = 6.5;               % [旋钮] 天线 y(图上 z)方向尺寸(格)

    % ---- 单位 1 网格设置(网格线每 1 格一条, 数字只在 5 的倍数标注) ----
    setUnitGrid = @(xv, zv) dealTicks(xv, zv);

    % ==================== 图7: 物理阵列 ====================
    % TX 红三角、RX 蓝方块 + 各自封装矩形(观察干涉); x-z 平面(整数栅格, 单位 1 网格)
    figPhys = figure('Position', figPos, 'Name', 'fig7_PhysArray');
    for i = 1:nTX                                  % TX 封装矩形(红)
        rectangle('Position', [tInt(1,i)-antSizeX/2, tInt(3,i)-antSizeZ/2, ...
                               antSizeX, antSizeZ], 'EdgeColor', 'r', 'LineWidth', 1);
    end
    hold on;
    for i = 1:nRX                                  % RX 封装矩形(蓝)
        rectangle('Position', [rInt(1,i)-antSizeX/2, rInt(3,i)-antSizeZ/2, ...
                               antSizeX, antSizeZ], 'EdgeColor', 'b', 'LineWidth', 1);
    end
    scatter(tInt(1, :), tInt(3, :), 120, 'r^', 'filled');
    scatter(rInt(1, :), rInt(3, :), 120, 'bs', 'filled');
    text(tInt(1, :), tInt(3, :) + 1.5, compose('TX%d', 1:nTX), ...
         'Color', 'r', 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    text(rInt(1, :), rInt(3, :) + 1.5, compose('RX%d', 1:nRX), ...
         'Color', 'b', 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    axis equal; grid on;
    setUnitGrid([tInt(1,:) rInt(1,:)], [tInt(3,:) rInt(3,:)]);
    xlabel('x (×FFT_d0)');
    ylabel('z (×FFT_d0)');
    title('物理阵列（设计原始整数栅格, 单位 1 网格, 矩形=天线封装(观察干涉)）');
    legend({'TX', 'RX'}, 'Location', 'best');

    % ==================== 图8: 虚拟阵列 ====================
    % 64 个虚拟阵元(TX+RX 相位中心矢量和), 每元标注 TX几RX几;
    % 列顺序与 A_manifold 一致: TX 内循环 / RX 外循环
    idx = 0;
    lab = cell(1, nTX*nRX);
    for j = 1:nRX
        for i = 1:nTX
            idx = idx + 1;
            lab{idx} = sprintf('T%dR%d', i, j);
        end
    end
    figVirt = figure('Position', figPos, 'Name', 'fig8_VirtArray');
    scatter(vInt(1, :), vInt(3, :), 50, 'ko', 'filled');
    % 每一行(z 相同)内按 x 排序, 相邻阵元的标注一上一下交替, 避免横向字重叠;
    % 水平方向居中对齐(不做 X 偏移)
    dy = zeros(1, nTX*nRX);
    zRows = unique(vInt(3, :));
    for r = 1:numel(zRows)
        rowIdx = find(vInt(3, :) == zRows(r));
        [~, ord] = sort(vInt(1, rowIdx));
        rowIdx = rowIdx(ord);                        % 行内按 x 排序
        upDown = 1 - 2*mod(0:numel(rowIdx)-1, 2);    % +1 上 / -1 下 交替
        dy(rowIdx) = upDown;
    end
    text(vInt(1, :), vInt(3, :) + dy*1.2, lab, 'FontSize', 5.5, ...
         'HorizontalAlignment', 'center');
    axis equal; grid on;
    setUnitGrid(vInt(1, :), vInt(3, :));
    xlabel('x (×FFT_d0)');
    ylabel('z (×FFT_d0)');
    title('虚拟阵列（设计原始整数栅格, TX×RX 矢量和, 单位 1 网格）');
end

% ----------------------------------------------------------------
% 单位 1 网格: 每 1 格一条网格线; 刻度数字只在 5 的倍数处显示(防拥挤)
% ----------------------------------------------------------------
function dealTicks(xv, zv)
    gx = floor(min(xv)) : ceil(max(xv));
    gz = floor(min(zv)) : ceil(max(zv));
    lx = arrayfun(@(v) sprintf('%d', v), gx, 'UniformOutput', false);
    lz = arrayfun(@(v) sprintf('%d', v), gz, 'UniformOutput', false);
    lx(mod(gx, 5) ~= 0) = {''};
    lz(mod(gz, 5) ~= 0) = {''};
    set(gca, 'XTick', gx, 'XTickLabel', lx, 'YTick', gz, 'YTickLabel', lz, ...
             'FontSize', 8, 'GridAlpha', 0.35);
end
