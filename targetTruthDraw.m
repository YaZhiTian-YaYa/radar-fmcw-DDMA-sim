function figTruth = targetTruthDraw(targets, prm, dims)
%TARGETTRUTHDRAW 画【真值】目标图(图23; 与 DBF/OMP/DML 三张点云图同口径, 供对比观察)
%   figTruth = targetTruthDraw(targets, prm, dims)
%   用场景放入的真实目标信息(gate_R/gate_D/azi/ele/rcsDb)画 6 子图, 版式与
%   fig28_DBF / fig30_DBF_OMP / fig32_DBF_OMP_DML 完全一致 —— 对比这三张点云图
%   即可看出: 哪些目标检测出来了、哪些漏检、点云里哪些是虚假点。
%   真值量: 距离 = gate_R×ΔR; 速度谱线 = mod(gate_D + 步进频耦合, 512) —— 与
%   检测点云同基准(耦合每距离 bin 约 +1 个多普勒 bin); 角度 = 设置值。
%   x-y 子图上标注目标序号(与《目标场景设计说明.docx》编号一致); 红星 = 校准角反。
    nTgtT = numel(targets);
    R0   = [targets.gate_R] * dims.DeltaR;                     % 真值距离 (m)
    dopC = zeros(1, nTgtT);
    for k = 1:nTgtT                                            % 谱线位置(含步进频耦合)
        dopC(k) = mod(targets(k).gate_D + prm.cr*prm.B/(dims.dopBinNum-1) * ...
                      (2*targets(k).gate_R*3e8/(2*prm.B))/3e8 * dims.dopBinNum, ...
                      dims.dopBinNum);
    end
    aziT = [targets.azi];  eleT = [targets.ele];
    rcsT = [targets.rcsDb];  isC = [targets.isCorner] == 1;
    zT = R0 .* sind(eleT);
    xT = R0 .* cosd(eleT) .* sind(aziT);
    yT = R0 .* cosd(eleT) .* cosd(aziT);
    rngT  = [targets.gate_R] * dims.DeltaR;                % 真值距离(用于轴范围)
    rMax10 = max(ceil(max(rngT) * 1.15 / 10) * 10, 10);    % 随数据, 15% 余量, 10 倍数
    sym10  = @(v) max(ceil(max(abs(v)) / 10) * 10, 10);   % 全 0 时保底 ±10, 防 xlim([0,0]) 非法

    screenSize = get(0, 'ScreenSize');
    figW = screenSize(3) * 0.98;  figH = screenSize(4) * 0.80;
    figPos = [(screenSize(3)-figW)/2, (screenSize(4)-figH)/2, figW, figH];
    figTruth = figure('Position', figPos, 'Name', 'fig33_TargetTruth');
    tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile;
    scatter(xT, yT, 24, zT, 'filled'); hold on;
    scatter(xT(isC), yT(isC), 80, 'r*', 'LineWidth', 1);       % 校准角反标红星
    text(xT, yT, arrayfun(@num2str, 1:nTgtT, 'UniformOutput', false), ...
         'FontSize', 7);                                       % 目标序号(对照场景文档)
    grid on; colorbar;
    xlim([-sym10(xT), sym10(xT)]); ylim([0, rMax10]);
    xlabel('x (m)'); ylabel('y (m)'); title('x-y 位置 (颜色=高度z, 标注=目标序号)');
    nexttile;
    scatter(aziT, R0, 24, zT, 'filled'); hold on;
    scatter(aziT(isC), R0(isC), 80, 'r*', 'LineWidth', 1);
    grid on; colorbar;
    xlim([-sym10(aziT), sym10(aziT)]); ylim([0, rMax10]);
    xlabel('方位角 (deg)'); ylabel('距离 (m)'); title('方位角-距离');
    nexttile;
    scatter(dopC, [targets.gate_R], 24, zT, 'filled'); grid on; colorbar;
    xlim([0, dims.dopBinNum-1]); ylim([0, dims.rngUse]);
    xlabel('多普勒 bin'); ylabel('距离 bin'); title('多普勒-距离 (满bin, 含耦合)');
    nexttile;
    scatter(zT, yT, 24, zT, 'filled'); grid on; colorbar;
    xlim([-sym10(zT), sym10(zT)]); ylim([0, rMax10]);
    xlabel('z (m)'); ylabel('y (m)'); title('z-y 位置');
    nexttile;
    scatter(eleT, R0, 24, zT, 'filled'); grid on; colorbar;
    xlim([-sym10(eleT), sym10(eleT)]); ylim([0, rMax10]);
    xlabel('俯仰角 (deg)'); ylabel('距离 (m)'); title('俯仰角-距离');
    nexttile;
    scatter3(xT, yT, zT, 24, zT, 'filled'); grid on;
    xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)'); view(3);
    title('目标 xyz (可旋转)');
    sgtitle(sprintf('真值目标图: 共 %d 个(颜色=高度z, 红星=校准角反, 序号同《目标场景设计说明》) —— 与 fig17/19/21 点云对比看漏检与虚假点', nTgtT));
end
