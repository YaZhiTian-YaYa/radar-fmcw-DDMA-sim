function overlayNoiseSurf(figH, floorDbCh)
%OVERLAYNOISESURF RD 3D 曲面各 RX 子图叠加"噪底 vs 距离"红线(与 95 仿真同款, 独立函数版)
%   overlayNoiseSurf(figH, floorDbCh)   floorDbCh: [nR×nRX] 逐距离门噪底(dB, estNoiseFloor* 的 floorDbCh 口径)
%   前(多普勒0)后(多普勒末)两缘各画一条: 前缘粗红实线, 后缘细红虚线; 子图按标题"RX n"匹配通道。
if isempty(figH) || ~ishandle(figH), return; end
nRx = size(floorDbCh, 2);
nR  = size(floorDbCh, 1);
x   = 0:nR-1;
for ax = reshape(findall(figH, 'Type', 'axes'), 1, [])
    ttl = ''; try, ttl = ax.Title.String; catch, end
    if isempty(regexp(ttl, 'RX\s*\d+', 'once')), continue; end
    rxI = str2double(regexp(ttl, 'RX\s*(\d+)', 'tokens', 'once'));
    if ~isscalar(rxI) || isnan(rxI) || rxI > nRx, continue; end
    yd  = ax.YLim;
    hold(ax, 'on');
    plot3(ax, x, zeros(1, nR),  floorDbCh(:, rxI), 'r-',  'LineWidth', 1.6);
    plot3(ax, x, yd(2)*ones(1, nR), floorDbCh(:, rxI), 'r--', 'LineWidth', 0.8);
end
axs = reshape(findall(figH, 'Type', 'axes'), 1, []);
if ~isempty(axs)
    text(axs(end), 0.02, 0.02, '红线 = 逐距离门噪底(中位)  距离bin→', 'Units', 'normalized', ...
        'FontSize', 9, 'Color', [0.75 0.05 0.05], 'BackgroundColor', [1 1 1 0.6]);
end
end
