function [tgtOut, detCfg] = doaEstimate_pipe(detCfg, tgtOut, AntCaliData, rowLut, lambda, dims, drawFigs)
%DOAESTIMATE_PIPE 方位DBF + 子阵俯仰DBF测角, 结果写回目标列表
%   (移植自 参考\rsp_pipeline_ReWrite\pipeline\doaEstimate_pipe.m, 变量名与参考一致;
%    两处适配: AntCaliData 由主函数 genAntCaliData 生成后传入(不再 load mat);
%    加密的 calValidChannel_pipe(子阵通道选择)用我们自己的阵列流型 rowLut 替代 ——
%    参考阵列方位最长行 30 元/俯仰 12 行, 我们的最长行 36 元/俯仰 9 行, 均读数据得出)
%   输入: detCfg.doaChanIn -- 各目标64通道×8chirp复数快照(每目标一行)
%         tgtOut          -- 检测输出的目标列表
%         AntCaliData     -- 虚拟导向矢量表(genAntCaliData 输出, 结构体名同参考)
%         rowLut          -- 通道俯仰行查找表 [行数 × 64]
%         lambda          -- 工作波长(主函数算好传入)
%         dims            -- 公共维度
%   drawFigs: false = 只测角不画(校准遍只取角反方位, 主流程两遍测量用)
%   输出: tgtOut 中补充 azimuth/elevation/doaMag/azimuthList/sideLobe
    if nargin < 7 || isempty(drawFigs), drawFigs = true; end   % 默认画图
    aziStartIdx = 1;                              % 方位导向矢量起始索引(对应网格起点 -75°)
    subArrTotal = size(rowLut, 1);                % 子阵总数(读数据: 我们 = 9 行, 参考 = 12 行)
    chanTotal = dims.chanTotal;
    aziRange = AntCaliData.aziSteerVecRange;
    eleRange = AntCaliData.eleSteerVecRange;
    aziGrid = aziRange(1) : aziRange(3) : aziRange(2);
    eleGrid = eleRange(1) : eleRange(3) : eleRange(2);   % 方位/俯仰搜索网格(1°步进)
    subArrUsed = subArrTotal;                               % 俯仰测角实际使用的子阵数
    % 方位DBF所用子阵: 取【元素最多的一行】(= 方位最长行, 我们 = 36 元, 读数据选)
    [~, aziRowBest] = max(sum(rowLut, 2));
    % ---------- 预计算各子阵的俯仰组导向矢量 ----------
    % 取该子阵所有通道的俯仰导向矢量求和后归一化, 作为子阵级俯仰导向矢量
    eleGrpSteerVec = zeros(subArrTotal, length(eleGrid));      % 各子阵的俯仰组导向矢量
    for subIdx = 1 : subArrTotal
        chanSel = find(rowLut(subIdx, :));                     % (替代 calValidChannel_pipe)
        grpSum = sum(AntCaliData.eleSteerVec(chanSel,:), 1);
        eleGrpSteerVec(subIdx, :) = grpSum ./ abs(grpSum);
    end
    %% ---------- 方位DBF: 对全部目标一次性做波束形成 ----------
    % 取目标快照中方位最长行子阵的通道, 与方位导向矩阵相乘,
    % 得到每个目标在方位网格上的DBF谱(dB)
    chanSel = find(rowLut(aziRowBest, :));                     % (替代 calValidChannel_pipe(1,64))
    aziDbfOut = detCfg.doaChanIn(:, chanSel) * 1 / length(chanSel) * ...
                    AntCaliData.aziSteerVec(chanSel, aziStartIdx : aziStartIdx + length(aziGrid) - 1);
    doaRes.aziDbfSpec = db(abs(aziDbfOut));       % [目标数 × 方位格] (内部调试/分析量, 不作为函数输出)
    doaRes.eleDbfSpec = zeros(size(tgtOut, 1), length(eleGrid));
    doaRes.item = struct('aziPos', [], 'elePos', [], 'aziSubOff', [], 'eleSubOff', [], 'aziPow', [], 'sideLobeLv', []);
    doaRes.item = repmat(doaRes.item, size(tgtOut, 1), 1);
    % ---------- 逐目标在方位谱上找峰(前4个), 并记录亚格偏移与副瓣电平 ----------
    for tgtIdx = 1 : size(tgtOut, 1)
        [pkPow, pkPos, pkSubOff] = searchPeaksPwr_pipe(doaRes.aziDbfSpec(tgtIdx, :).', 4);
        doaRes.item(tgtIdx).aziPos = single(pkPos);         % 前4个方位峰的格索引
        doaRes.item(tgtIdx).aziSubOff = single(pkSubOff);   % 各峰的亚格偏移
        doaRes.item(tgtIdx).aziPow = single(pkPow);         % 各峰功率
        doaRes.item(tgtIdx).sideLobeLv = single(max(pkPow) - pkPow);   % 相对主峰的副瓣电平
    end
    % ---------- 逐目标俯仰测角: 先子阵合成, 再子阵级俯仰DBF ----------
    for tgtIdx = 1 : size(tgtOut, 1)

        subArrOut = zeros(1, subArrUsed);
        % 每个子阵: 用目标主方位峰处的导向矢量对本子阵通道做匹配滤波(子阵合成),
        % 把相位参考统一到目标方位上
        for subIdx = 1 : subArrUsed
            steerPos = doaRes.item(tgtIdx).aziPos(1) + aziStartIdx - 1;    % 主峰对应的导向矢量索引
            chanSel = find(rowLut(subIdx, :));                             % (替代 calValidChannel_pipe)
            steerVec = AntCaliData.aziSteerVec(chanSel, steerPos).';
            subArrOut(subIdx) = steerVec * detCfg.doaChanIn(tgtIdx, chanSel).';
        end

        % 子阵输出归一(单位模, 只留相位; 2026-10-30): 俯仰谱 = 等幅孔径口径。
        % 数据侧与导向侧按行归一数学等价(同一标量乘同一项, 实测差 <1e-14), 按约定放数据侧。
        % 不归一时元素最多的行(36 元主行)的自然幅度会把俯仰副瓣压低约 6 dB
        % (实测 B1 阵 @(0°,0°): -2.8 dB vs 等幅 -8.9 dB), 方向图/副瓣类评估口径偏乐观。
        subArrOut = subArrOut ./ max(abs(subArrOut), eps);

        % 子阵输出 × 子阵俯仰导向矢量 -> 俯仰谱, 取第1个峰(幅度谱findpeaks)
        eleDbfBuf = subArrOut * eleGrpSteerVec(1:subArrUsed,:);
        pkPow = zeros(4, 1);
        pkPos = zeros(4, 1);
        pkSubOff = zeros(size(eleDbfBuf, 1), 1);

        elePow = db(abs(eleDbfBuf));
        [pkPow, pkPos, pkSubOff] = searchPeaksPwr_pipe(elePow.', 1);

        doaRes.item(tgtIdx).elePos = single(pkPos);         % 俯仰峰格索引
        doaRes.item(tgtIdx).eleSubOff = single(pkSubOff);   % 俯仰亚格偏移
    end
    % 方位谱逐目标归一化(减去各自最大值), 便于对比
    doaRes.aziDbfSpec = single(doaRes.aziDbfSpec - repmat(max(doaRes.aziDbfSpec, [], 2), 1, size(doaRes.aziDbfSpec, 2)));

    % ---------- 结果写回目标列表 ----------
    for tgtIdx = 1 : size(tgtOut, 1)
        % 方位/俯仰角 = 网格角度 + 亚格偏移(一阶精化)
        tgtOut(tgtIdx).azimuth = aziGrid(doaRes.item(tgtIdx).aziPos(1)) + doaRes.item(tgtIdx).aziSubOff(1);
        tgtOut(tgtIdx).elevation = eleGrid(doaRes.item(tgtIdx).elePos(1)) + doaRes.item(tgtIdx).eleSubOff(1);

        tgtOut(tgtIdx).doaMag = doaRes.item(tgtIdx).aziPow(1);              % 主方位峰功率
        tgtOut(tgtIdx).azimuthList = aziGrid(doaRes.item(tgtIdx).aziPos).'; % 前4个候选方位列表
        tgtOut(tgtIdx).sideLobe = doaRes.item(tgtIdx).sideLobeLv.';         % 对应副瓣电平
    end

    %% ---------- 【新增】目标点迹五子图(颜色=高度z) ----------
    DeltaR = 3e8 / (2 * dims.B);                         % 距离分辨率(B 由主函数经 dims 传入)
    % 点迹坐标用【亚bin插值后】的距离/速度(tgtOut.range/speed), 不再用整数 bin
    rngEst = cell2mat(arrayfun(@(t) double(t.range), tgtOut, 'UniformOutput', false));              % 插值后距离 (m)
    rngBinEst = rngEst / DeltaR;                          % 插值后距离 bin(小数 = 亚bin, 0~512)
    dopEst = cell2mat(arrayfun(@(t) double(t.dopplerIdx) - 1 - double(t.doppler_alpha), tgtOut, 'UniformOutput', false)); % 插值后多普勒 bin(0基小数, 光标直读=真bin号)
    aziEst = cell2mat(arrayfun(@(t) double(t.azimuth),  tgtOut, 'UniformOutput', false));
    eleEst = cell2mat(arrayfun(@(t) double(t.elevation), tgtOut, 'UniformOutput', false));
    zEst   = rngEst .* sind(eleEst);                     % 高度 z = 距离 × sin(俯仰)
    xPos   = rngEst .* cosd(eleEst) .* sind(aziEst);
    yPos   = rngEst .* cosd(eleEst) .* cosd(aziEst);
    detCfg.tgtXYZ = [xPos; yPos; zEst];                  % 供后续使用

    % 坐标轴设置: 横轴左右对称(角度) / 纵轴满量程, 均取 10 的倍数;
    % 距离/速度子图按【量程】画: 速度 bin 0~511, 距离 bin 0~512
    rngMax    = dims.rngUse * 3e8 / (2 * dims.B);        % 最大量程 (m, = 512×ΔR ≈ 299)
    rMax10    = max(ceil(max(rngEst) * 1.15 / 10) * 10, 10);  % 纵轴随检测数据, 15% 余量
    sym10     = @(v) max(ceil(max(abs(v)) / 10) * 10, 10);  % 横轴对称限(10 的倍数, 全 0 保底 ±10 防 xlim 非法)

    detCfg.figDOA = [];                                  % 不画图时句柄为空(figList 自动跳过)
    if drawFigs
    screenSize = get(0, 'ScreenSize');
    figW = screenSize(3) * 0.98;  figH = screenSize(4) * 0.80;          % 与图5(RD 3D)同尺寸
    figPos = [(screenSize(3)-figW)/2, (screenSize(4)-figH)/2, figW, figH];
    detCfg.figDOA = figure('Position', figPos, 'Name', 'fig28_DBF');
    tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile;   % 子图1: x-y 位置
    scatter(xPos, yPos, 24, zEst, 'filled'); grid on; colorbar;
    xlim([-sym10(xPos), sym10(xPos)]); ylim([0, rMax10]);
    xlabel('x (m)'); ylabel('y (m)');
    title('目标 x-y 位置 (颜色=高度z)');

    nexttile;   % 子图2: 方位角-距离
    scatter(aziEst, rngEst, 24, zEst, 'filled'); grid on; colorbar;
    xlim([-sym10(aziEst), sym10(aziEst)]); ylim([0, rMax10]);
    xlabel('方位角 (deg)'); ylabel('距离 (m)');
    title('方位角-距离 (颜色=高度z)');

    nexttile;   % 子图3: 多普勒-距离 (量程满 bin: 速度 0~511, 距离 0~512; 点=亚bin插值位置)
    scatter(dopEst, rngBinEst, 24, zEst, 'filled'); grid on; colorbar;
    xlim([0, dims.dopBinNum-1]); ylim([0, dims.rngUse]);
    xlabel('多普勒 bin'); ylabel('距离 bin');
    title('多普勒-距离 (颜色=高度z)');

    nexttile;   % 子图4: z-y 位置
    scatter(zEst, yPos, 24, zEst, 'filled'); grid on; colorbar;
    xlim([-sym10(zEst), sym10(zEst)]); ylim([0, rMax10]);
    xlabel('z (m)'); ylabel('y (m)');
    title('目标 z-y 位置 (颜色=高度z)');

    nexttile;   % 子图5: 俯仰角-距离
    scatter(eleEst, rngEst, 24, zEst, 'filled'); grid on; colorbar;
    xlim([-sym10(eleEst), sym10(eleEst)]); ylim([0, rMax10]);
    xlabel('俯仰角 (deg)'); ylabel('距离 (m)');
    title('俯仰角-距离 (颜色=高度z)');

    nexttile;   % 子图6: 点云 xyz 三维视图(鼠标可旋转)
    scatter3(xPos, yPos, zEst, 24, zEst, 'filled'); grid on; colorbar;
    xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
    title('目标点云 xyz (颜色=高度z, 鼠标可旋转)');
    view(3); axis equal; axis vis3d;
    rotate3d on;                                         % 激活鼠标旋转(桌面交互)
    sgtitle('DOA 测角结果点迹 (颜色=高度z)');
    end                                             % if drawFigs
end
