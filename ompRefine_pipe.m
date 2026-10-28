function [tgtOut, detCfg, ompCfg] = ompRefine_pipe(tgtOut, detCfg, AntCaliData, rowLut, dims, thr2ndDb, seedR, seedD, drawFigs)
%OMPREFINE_PIPE 对全部过两级 CFAR 的点云逐点做 OMP 角度超分辨(最多两个目标)
%   [tgtOut, detCfg, ompCfg] = ompRefine_pipe(tgtOut, detCfg, AntCaliData, rowLut, ...
%                                             dims, thr2ndDb, seedR, seedD)
%   输入:
%     tgtOut      两级 CFAR 确认的点云列表(每点含 azimuth = DOA 插值方位角)
%     detCfg      doaChanIn: 各点 64 通道快照(每行一点, 与 tgtOut 行序一致)
%                 (过程图用的方位 DBF 谱在本模块内部按 DOA 模块同式自算 —— 纯调试/分析量不跨模块传递)
%     AntCaliData 虚拟导向矢量表(1° 间隔, 工程采集格式)
%     rowLut      通道俯仰行查找表(取最长行通道做方位维)
%     dims        公共维度
%     thr2ndDb    【第二目标判定门限】残差匹配曲线峰高出噪底中位的 dB 数
%     seedR/seedD 指定做 OMP 过程图的目标的 距离bin/多普勒bin(1起, 与 tgtOut 同基准;
%                 容差 ±1 bin 匹配)。【可选】: 不传或传 [] → 不画过程图
%   处理(逐点依次, 最多分辨两个目标, 不判第三 —— 车载已足够):
%     0) 两级门限: 处理级(宽, ±aziProcDeg/xProcM, 性能较好范围)内才做 OMP ——
%        关心该范围内 DBF 融合目标有没有漏检; 第二目标通过判定后再过输出级
%        (窄, ±aziLimDeg/xLimM)复核, 范围外删除(不关注, 减少测错杂点);
%     1) 第一原子: 该点 DOA 插值方位角 → steerFitAt 表格插值导向矢量;
%     2) LS 扣除第一目标 → 残差在 1° 表上匹配 → 抛物线亚格 → 第二角度候选;
%     3) 判定: 残差匹配曲线峰高出噪底中位 > thr2ndDb → 有第二目标;
%     4) 有则拟合第二角度导向矢量, 两原子联合最小二乘 → 最终权重;
%        第二目标作为新点追加进点云列表(俯仰等字段继承第一目标)。
%   输出:
%     tgtOut  更新后的点云列表(全部点新增 isOMP2nd/ompGate/ompHas2nd/ompAng2/
%             ompW1Db/ompW2Db 字段; 通过判定的第二目标已追加, azimuth = 第二角度)
%     detCfg  原样返回
%     ompCfg  .figSeed   指定点的 OMP 过程图句柄(未指定/未找到时为 [])
%             .figCloud  更新点云的 6 子图 DOA 图句柄
%             .n2nd      通过判定且过输出门限、已追加点云的第二目标数
%             .n2ndDrop  通过判定但超出输出门限被删的第二目标数
%             .nGate     因处理门限(宽)跳过 OMP 的点数
%             .seedFound 指定点是否找到

    % ---------- 方位字典(最长行 36 通道, 共轭域; 与 DOA 方位维一致) ----------
    aziRange = AntCaliData.aziSteerVecRange;
    aziGrid  = aziRange(1) : aziRange(3) : aziRange(2);
    [~, aziRowBest] = max(sum(rowLut, 2));
    chanAzi  = find(rowLut(aziRowBest, :));
    A_omp    = conj(AntCaliData.aziSteerVec(chanAzi, :));    % 信号模型域字典
    nTgt  = numel(tgtOut);
    n2nd  = 0;
    n2ndDrop = 0;                             % 通过判定但超出输出门限被删的第二目标数
    nGate = 0;
    seedFound = false;
    seed = [];
    if nargin < 8, seedD = []; end
    if nargin < 7, seedR = []; end
    if nargin < 9 || isempty(drawFigs), drawFigs = true; end   % false = 批量扫描不画图
    haveSeed = ~isempty(seedR) && ~isempty(seedD);

    % 统一新增标记字段(先建字段, 再用 deal 对【全部元素】填默认值 —— 只写 tgtOut(1)
    % 时其余元素该字段为空 [], 空值在 [tgtOut.字段] 拼接时丢元素, 后级 DML 掩码会错位)
    tgtOut(1).isOMP2nd  = 0;
    tgtOut(1).ompGate   = 1;
    tgtOut(1).ompHas2nd = 0;
    tgtOut(1).omp2ndKeep = 0;
    tgtOut(1).ompAng2   = NaN;
    tgtOut(1).ompW1Db   = NaN;
    tgtOut(1).ompW2Db   = NaN;
    [tgtOut.isOMP2nd]  = deal(0);
    [tgtOut.ompGate]    = deal(1);
    [tgtOut.ompHas2nd]  = deal(0);
    [tgtOut.omp2ndKeep] = deal(0);
    [tgtOut.ompAng2]    = deal(NaN);
    [tgtOut.ompW1Db]    = deal(NaN);
    [tgtOut.ompW2Db]    = deal(NaN);

    % ---------- OMP 两级门限 ----------
    % 处理级(宽): 性能较好范围 ±50°、|x|≤30m 内的点【都做】OMP —— 关心的是这个
    %   范围内 DBF 融合的目标有没有漏检(处理范围宁大勿小)
    aziProcDeg = 50;                    % [旋钮] OMP 处理方位范围: |DBF拟合方位角| ≤ 此值才做 OMP
    xProcM     = 100;                    % [旋钮] OMP 处理 x 范围 (m)
    % 输出级(窄): OMP 第二目标通过判定后, 还须落在 ±20°、|x|≤25m 内才保留进点云;
    %   范围外真有目标也不关注, 删掉 —— 减少测错(旁瓣/噪声残余误判)产生的杂点
    aziLimDeg = 20;                     % [旋钮] OMP 第二目标输出方位门限 (deg)
    xLimM     = 25;                     % [旋钮] OMP 第二目标输出 x 门限 (m)
    DeltaR    = 3e8/(2*dims.B);         % 距离分辨率 (m)

    for t = 1:nTgt
        x    = double(detCfg.doaChanIn(t, :)).';
        x    = x(chanAzi);                                   % 方位维 36 通道
        azi1 = double(tgtOut(t).azimuth);                    % DBF 拟合后的第一目标方位角

        % ---- OMP 处理门限(宽): 性能较好范围外不做 OMP, 不产生第二目标 ----
        rngM  = double(tgtOut(t).rangIdx) * DeltaR;
        eleM  = double(tgtOut(t).elevation);
        xPosT = rngM * cosd(eleM) * sind(azi1);
        if abs(azi1) > aziProcDeg || abs(xPosT) > xProcM
            tgtOut(t).ompGate = 0;
            nGate = nGate + 1;
            continue;
        end

        % ---- 第一原子: DOA 插值方位角的导向矢量拟合 ----
        atom1 = steerFitAt(A_omp, aziGrid, azi1);

        % ---- LS 扣除第一目标 ----
        a1  = (atom1' * x) / (atom1' * atom1);
        res = x - atom1 * a1;

        % ---- 残差在 1° 表上匹配 → 第二角度候选(抛物线亚格) ----
        c2   = A_omp' * res;
        c2Db = 20*log10(abs(c2) + eps);
        pk2  = max(c2Db);  med2 = median(c2Db);
        has2 = (pk2 - med2) > thr2ndDb;

        [~, p2] = max(abs(c2));
        if p2 > 1 && p2 < numel(c2)
            d2 = (abs(c2(p2-1)) - abs(c2(p2+1))) / ...
                 (2*(abs(c2(p2-1)) + abs(c2(p2+1))) - 4*abs(c2(p2)));
        else
            d2 = 0;
        end
        ang2 = aziGrid(p2) + d2;

        % ---- 第二角度导向矢量拟合 + 联合 LS ----
        atom2 = steerFitAt(A_omp, aziGrid, ang2);
        W2    = pinv([atom1, atom2]) * x;

        % ---- 写回逐点结果 ----
        tgtOut(t).ompHas2nd = has2;
        tgtOut(t).ompAng2   = ang2;
        tgtOut(t).ompW1Db   = 20*log10(abs(a1) + eps);
        tgtOut(t).ompW2Db   = 20*log10(abs(W2(2)) + eps);

        % ---- 判定有第二目标 → 输出门限(窄)复核 → 通过才追加点云列表 ----
        % 第二目标自身位置(角度用 OMP 第二角度, 距离/俯仰继承母点)落在
        % ±aziLimDeg、|x|≤xLimM 之外 → 删除(范围外不关注, 减少测错杂点)
        keep2 = false;
        if has2
            x2    = rngM * cosd(eleM) * sind(ang2);          % OMP 第二目标 x 位置 (m)
            keep2 = abs(ang2) <= aziLimDeg && abs(x2) <= xLimM;
            if keep2
                n2nd = n2nd + 1;
                newT          = tgtOut(t);                   % 继承同距同速/俯仰等字段
                newT.ID       = max([tgtOut.ID]) + 1;
                newT.azimuth  = ang2;                        % 第二目标方位
                newT.isOMP2nd = 1;
                tgtOut(end+1) = newT;
            else
                n2ndDrop = n2ndDrop + 1;                     % 通过判定但超出输出范围 → 删
            end
        end
        tgtOut(t).omp2ndKeep = keep2;

        % ---- 指定点(R/D 序号, 容差 ±1 bin) → 保存 OMP 过程数据 ----
        % 同一 RD 窗口常有多个检测点: 首个命中先记, 之后若遇到"判出 OMP 第二目标"
        % 的命中点则升级替换 —— 过程图优先展示成功分离的点, 避免图上写
        % "无OMP第二目标"而同 RD 的点云里却有 OMP 第二目标(看似自相矛盾)
        seedHit = haveSeed && abs(tgtOut(t).rangIdx - seedR) <= 1 && ...
                  abs(tgtOut(t).dopplerIdx - seedD) <= 1;
        if seedHit && (~seedFound || (~seed.has2 && has2))
            seedFound = true;
            seed = struct('t', t, 'azi1', azi1, 'ang2', ang2, 'has2', has2, ...
                          'keep2', keep2, 'c2Db', c2Db, 'med2', med2, ...
                          'atom1', atom1, 'atom2', atom2, 'W', W2, 'x', x);
        end
    end
    if haveSeed && ~seedFound
        fprintf('OMP 提示: RD(%g, %g) 处未找到确认点, 过程图跳过\n', seedR, seedD);
    end

    % ---------- 指定点的 OMP 过程图(图22, 3 格; 未指定/未找到则为空) ----------
    screenSize = get(0, 'ScreenSize');
    figW = screenSize(3) * 0.98;  figH = screenSize(4) * 0.80;
    figPos = [(screenSize(3)-figW)/2, (screenSize(4)-figH)/2, figW, figH];
    if seedFound && drawFigs
        if seed.has2 && seed.keep2,     vTxt = '有OMP第二目标';
        elseif seed.has2,               vTxt = '有OMP第二目标(超出输出范围,已删)';
        else,                           vTxt = '无OMP第二目标'; end
        figSeed = figure('Position', figPos, 'Name', 'fig29_OMP_process');
        tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
        nexttile;   % 常规 DBF 谱(内部自算, 与 DOA 模块同式: 快照/通道数 × 方位导向矩阵)
        dbfDb = db(abs(detCfg.doaChanIn(seed.t, chanAzi) / numel(chanAzi) * ...
                       AntCaliData.aziSteerVec(chanAzi, :)));       % [1 × 方位格]
        dbfDb = dbfDb - max(dbfDb);
        plot(aziGrid, dbfDb, 'b-'); hold on;
        plot(seed.azi1, 0, 'bv', 'MarkerFaceColor', 'b');
        grid on; ylim([-50 2]);
        xlabel('方位角 (deg)'); ylabel('归一化幅度 (dB)');
        title(sprintf('常规 DBF 谱: DBF插值角 %.2f°', seed.azi1));
        nexttile;   % 残差匹配曲线(OMP第二目标)
        plot(aziGrid, seed.c2Db - max(seed.c2Db), 'r-'); hold on;
        plot(seed.ang2, 0, 'mv', 'MarkerFaceColor', 'm');
        yline(seed.med2 - max(seed.c2Db), 'k:');
        grid on; ylim([-50 2]);
        xlabel('方位角 (deg)'); ylabel('归一化幅度 (dB)');
        title(sprintf('扣除第一目标后的残差匹配: OMP第二目标 %.2f° (判定: %s)', ...
            seed.ang2, vTxt));
        nexttile;   % 联合重建检验
        xRef = max(abs(seed.x)) + eps;
        plot(20*log10(abs(seed.x)/xRef + eps), 'b-o', 'MarkerSize', 4); hold on;
        plot(20*log10(abs([seed.atom1, seed.atom2] * seed.W)/xRef + eps), 'r-*', 'MarkerSize', 4);
        grid on; ylim([-60 5]);
        xlabel('通道 (最长行 36)'); ylabel('相对快照峰值 (dB)');
        legend({'快照(归一化)', 'OMP 双原子重建(同一参考)'}, 'Location', 'best');
        title('重建检验: 快照 vs 双原子重建 (统一归一)');
        sgtitle(sprintf('OMP 角度超分辨: RD(%g,%g): DBF插值角 %.2f°, OMP第二目标 %.2f°, 判定: %s', ...
            seedR, seedD, seed.azi1, seed.ang2, vTxt));
    else
        figSeed = [];                                     % 未指定/未找到 → 无过程图
    end

    if ~drawFigs
        figSeed = [];  figCloud = [];                  % 批量扫描: 不画图
    end
    if drawFigs
    % ---------- 更新点云的 6 子图 DBF 点迹图(参考 fig28_DBF, 黑圈=OMP第二目标) ----------
    rngEst = cell2mat(arrayfun(@(t) double(t.range), tgtOut, 'UniformOutput', false));
    dopEst = cell2mat(arrayfun(@(t) double(t.dopplerIdx) - 1 - double(t.doppler_alpha), tgtOut, 'UniformOutput', false)); % 0基小数bin
    aziEst = cell2mat(arrayfun(@(t) double(t.azimuth),  tgtOut, 'UniformOutput', false));
    eleEst = cell2mat(arrayfun(@(t) double(t.elevation), tgtOut, 'UniformOutput', false));
    zEst   = rngEst .* sind(eleEst);
    xPos   = rngEst .* cosd(eleEst) .* sind(aziEst);
    yPos   = rngEst .* cosd(eleEst) .* cosd(aziEst);
    is2    = false(1, numel(tgtOut));  is2(nTgt+1:end) = true;   % 追加段 = 第二目标
    rMax10 = ceil(dims.rngUse * 3e8/(2*dims.B) / 10) * 10;
    sym10  = @(v) max(ceil(max(abs(v)) / 10) * 10, 10);   % 全 0 时保底 ±10, 防 xlim([0,0]) 非法

    figCloud = figure('Position', figPos, 'Name', 'fig30_DBF_OMP');
    tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    nexttile;
    scatter(xPos, yPos, 24, zEst, 'filled'); hold on;
    scatter(xPos(is2), yPos(is2), 40, 'ko'); grid on; colorbar;
    xlim([-sym10(xPos), sym10(xPos)]); ylim([0, rMax10]);
    xlabel('x (m)'); ylabel('y (m)'); title('x-y 位置 (黑圈=OMP第二目标)');
    nexttile;
    scatter(aziEst, rngEst, 24, zEst, 'filled'); hold on;
    scatter(aziEst(is2), rngEst(is2), 40, 'ko'); grid on; colorbar;
    xlim([-sym10(aziEst), sym10(aziEst)]); ylim([0, rMax10]);
    xlabel('方位角 (deg)'); ylabel('距离 (m)'); title('方位角-距离');
    nexttile;
    scatter(dopEst, rngEst/dims.DeltaR, 24, zEst, 'filled'); grid on; colorbar;
    xlim([0, dims.dopBinNum-1]); ylim([0, dims.rngUse]);
    xlabel('多普勒 bin'); ylabel('距离 bin'); title('多普勒-距离 (满bin)');
    nexttile;
    scatter(zEst, yPos, 24, zEst, 'filled'); hold on;
    scatter(zEst(is2), yPos(is2), 40, 'ko'); grid on; colorbar;
    xlim([-sym10(zEst), sym10(zEst)]); ylim([0, rMax10]);
    xlabel('z (m)'); ylabel('y (m)'); title('z-y 位置');
    nexttile;
    scatter(eleEst, rngEst, 24, zEst, 'filled'); hold on;
    scatter(eleEst(is2), rngEst(is2), 40, 'ko'); grid on; colorbar;
    xlim([-sym10(eleEst), sym10(eleEst)]); ylim([0, rMax10]);
    xlabel('俯仰角 (deg)'); ylabel('距离 (m)'); title('俯仰角-距离');
    nexttile;
    scatter3(xPos, yPos, zEst, 24, zEst, 'filled'); hold on;
    scatter3(xPos(is2), yPos(is2), zEst(is2), 40, 'ko'); grid on;
    xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)'); view(3);
    title('点云 xyz (可旋转)');
    sgtitle(sprintf('OMP 后点云: 原 %d 点 + OMP第二目标 %d 点 (黑圈)', nTgt, n2nd));
    end                                             % if drawFigs(云图段)

    ompCfg = struct('figSeed', figSeed, 'figCloud', figCloud, 'n2nd', n2nd, ...
                    'n2ndDrop', n2ndDrop, 'nGate', nGate, 'seedFound', seedFound);
end
