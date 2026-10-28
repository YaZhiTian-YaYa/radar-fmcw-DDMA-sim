function point_disp_pipe(tgtOut)
%POINT_DISP_PIPE 点云显示: 十字坐标系 / 角度-距离 / 速度-距离三视图
%   距离与速度均以 bin 为单位(未乘物理标定系数)
%   空目标列表直接返回
    if nargin < 1 || isempty(tgtOut) || size(tgtOut, 1) < 1
        return;
    end
    % 距离 = 距离bin索引 + 亚bin偏移
    rngBin = [tgtOut.rangIdx] + [tgtOut.range_alpha];
    % 多普勒折叠到[-256, 255](以0为中心, 区分正负速度)
    dopBin = mod([tgtOut.dopplerIdx] - 1 + 256, 512) - 256 + [tgtOut.doppler_alpha];
    azi    = [tgtOut.azimuth];                 % 方位角(度)
    xPos   = rngBin .* sind(azi);              % 前向直角坐标(以bin为单位)
    yPos   = rngBin .* cosd(azi);

    % 视图1: 十字坐标系(自车位于原点)
    figure(2001); clf;
    hold on;
    scatter(xPos, yPos, 16, 'k', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.3);
    xline(0, ':'); yline(0, ':');
    axis equal tight;
    xlabel('x/bin'); ylabel('y/bin'); grid on;
    title(['十字坐标系   N = ', num2str(numel(rngBin))], 'FontWeight', 'bold');

    % 视图2: 方位角-距离
    figure(2002); clf;
    hold on;
    scatter(azi, rngBin, 16, 'k', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.3);
    xlabel('azimuth/°'); ylabel('range/bin'); grid on;
    title('角度-距离', 'FontWeight', 'bold');

    % 视图3: 速度-距离
    figure(2003); clf;
    hold on;
    scatter(dopBin, rngBin, 16, 'k', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.3);
    xlabel('doppler/bin'); ylabel('range/bin'); grid on;
    title('速度-距离', 'FontWeight', 'bold');

    drawnow;
end
