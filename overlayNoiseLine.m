function overlayNoiseLine(figH, specs)
%OVERLAYNOISELINE 积累图曲面叠加"逐距离门噪底"黑线(各图用本域数据的中位, 尺度天然一致)
%   overlayNoiseLine(figH, specs)
%   specs: 结构体 cell 数组(字段集可不同), 每元素:
%     .match   标题包含此串的子图应用该曲线
%     .exclude (可选) 标题包含此串则跳过(如'非相参'与'相参'区分)
%     .nf      [nR×1] 逐距离门噪底(dB, 本域)
%     .x2      后缘横坐标(数据边缘; 前缘恒为 X 轴起点 0)
%   画法: 前缘(x=0)黑色粗实线 + 后缘(x=x2)黑色细虚线, 沿距离轴(Y)延伸;
%   每子图左下角黑字图例。噪底取"本图数据逐距离门中位"(对目标占据少数列稳健,
%   与 estNoiseFloor 的对数域中位口径一致)。
if isempty(figH) || ~ishandle(figH), return; end
if ~iscell(specs), specs = {specs}; end               % 单结构体调用归一化为 cell
for ax = reshape(findall(figH, 'Type', 'axes'), 1, [])
    ttl = ''; try, ttl = ax.Title.String; catch, end
    if isempty(ttl), continue; end
    for sC = specs
        s = sC{1};                                        % specs 为 cell: 解包出结构体
        if ~contains(ttl, s.match), continue; end
        if isfield(s, 'exclude') && contains(ttl, s.exclude), continue; end
        nR  = numel(s.nf);
        yr  = ax.YLim;
        hold(ax, 'on');
        plot3(ax, zeros(1, nR), yr(1) + [0:nR-1], s.nf, 'k-',  'LineWidth', 1.6);
        plot3(ax, s.x2*ones(1, nR), yr(1) + [0:nR-1], s.nf, 'k--', 'LineWidth', 0.8);
        text(ax, 0.02, 0.04, '黑线 = 逐距离门噪底(本域中位)', 'Units', 'normalized', ...
            'FontSize', 9, 'Color', [0 0 0], 'BackgroundColor', [1 1 1 0.6]);
        break;                                          % 每子图只应用第一条匹配
    end
end
end
