function [tgtOut, detCfg, dmlCfg] = dmlRefine_pipe(tgtOut, detCfg, AntCaliData, rowLut, dims, thrDmlDb, seedR, seedD, drawFigs, accumCfg, mcOn, robExp)
%DMLREFINE_PIPE 对全部点云(不含 OMP 第二目标)逐点做 DML 角度超分辨(最多两目标)
%   [tgtOut, detCfg, dmlCfg] = dmlRefine_pipe(tgtOut, detCfg, AntCaliData, rowLut, dims, thrDmlDb, seedR, seedD)
%   输入:
%     tgtOut      OMP 模块输出的点云列表(OMP 第二目标也在内 —— 本模块跳过它不做 DML)
%     detCfg      doaChanIn: 各点 64 通道快照(每行一点, 与 tgtOut 原始行序一致)
%     AntCaliData 虚拟导向矢量表(1° 间隔, 工程采集格式)
%     rowLut      通道俯仰行查找表(取最长行通道做方位维)
%     dims        公共维度
%     thrDmlDb    【第二目标判定门限】残差改善超过此值(dB)才输出第二目标
%     seedR/seedD 指定做 DML 过程图的目标的 距离bin/多普勒bin(1 起, 容差 ±1 bin;
%                 注意步进频耦合会使实际 D 序号偏离放置值, 按 RD 图实际位置填)
%     accumCfg    (可选, 遗留②多 cell 用) chanSnapshot/ddmOrdMap 来源; 空则跳过多 cell
%     mcOn        (可选) true=配对点生产权重改用多 cell 合成版(单 cell 值备份到 dmlW1DbSc)
%   处理(逐点依次, 算法逻辑与已验证的 dmlAngle 完全一致):
%     0) 跳过 OMP 第二目标; 门限检查(|DBF方位|≤90° 且 |x|≤300 m, 旋钮见下)超出跳过;
%     1) K=1 粗搜: 1° 表上最佳单角度(强目标);
%     2) 残差匹配曲线取候选峰(限于 DBF 主峰 ±dmlSpanDeg 内, 破 39° 边缘假峰)
%        → 各起点做交替坐标细分(±1° 内 dmlStepDeg 步进 ×3 轮, 表格插值原子);
%     3) 取最终代价最小的起点为解;
%     4) 判定: 残差改善 imprDb = 10log10(J1/J2) > thrDmlDb 且 sepMinDml(0.7°) < 间隔 ≤ sepLimDeg(6°,
%        与 OMP 衔接) 且 Gram 条件数 κ_G ≤ kapMaxDml(16, =原子κ4 与 pinv 封顶同口径)
%        【2026-11-06 收紧: 原仅 >0.2°, 强弱对角度塌缩(两原子拟到强目标一侧)会虚假分离】
%        → 通过: 第一目标方位替换为 DML 第一角度, 第二目标追加点云(继承俯仰等);
%        不通过: 该点保留 DML 之前的 DOA 结果(方位不动, 不追加)。
%   输出:
%     tgtOut  更新后的点云列表(全部原点新增 isDML2nd/dmlGate/dmlPass/dmlAng1/
%             dmlAng2/dmlW1Db/dmlW2Db/dmlW1/dmlW2 字段; 通过判定的第二目标已追加,
%             isDML2nd=1; dmlW1/dmlW2 为线性 LS 幅度, 供 RCS 推算)
%     detCfg  原样返回
%     dmlCfg  .figSeed   指定点的 DML 过程图(4 格; 未找到指定点时为 [])
%             .figCloud  更新点云的 6 子图 DOA 图(同 fig14 口径, 含 DML 第二目标)
%             .nDML2nd   判定出 DML 第二目标的点数
%             .nGate     因门限(角度/x位置)跳过 DML 的点数
%             .nDone     实际做了 DML 细分的点数
%             .seedFound 指定点是否找到

    % ---------- 方位字典(最长行 36 通道, 共轭域; 与 DOA/OMP 方位维一致) ----------
    aziRange = AntCaliData.aziSteerVecRange;
    aziGrid  = aziRange(1) : aziRange(3) : aziRange(2);
    [~, aziRowBest] = max(sum(rowLut, 2));
    chanAzi  = find(rowLut(aziRowBest, :));
    A_dml    = conj(AntCaliData.aziSteerVec(chanAzi, :));    % 信号模型域字典
    fitFun   = @(th) steerFitAt(A_dml, aziGrid, th);         % 表格插值原子(按需单向量)
    DeltaR   = 3e8/(2*dims.B);                          % 距离分辨率 (m)

    % ---------- DML 分工限定(与 OMP 无缝衔接: DML 管 1.5~6°, OMP 管 ≥5°, 5~6° 双保险) ----------
    dmlSpanDeg = 8;                     % [旋钮] DML 搜索半跨度: 第二角度候选限于 DBF 主峰 ±此值内 (deg, 覆盖 6° 分离+裕量)
    sepLimDeg  = 6;                     % [旋钮] 两目标角度间隔上限 (deg): 超出判不通过(曾为 4°, 与 OMP 的 6° 之间留 4~6° 空档, 已改 6° 衔接)
    sepMinDml  = 0.7;                   % [旋钮] DML 通过的最小角间隔 (deg)【2026-11-06 由 0.2 收紧】:
                                        %   亚半波束对(3dB宽1.88°/角bin 1.4°)分量幅度解算病态, 且实测
                                        %   强弱对发生角度塌缩(两原子都拟到强目标一侧, sep 0.35~0.6°
                                        %   仍给出>3dB假改善, D04/D05 取证) —— 拒绝分离, 输出融合单点
    thrImprSDml = 1.0;                  % [旋钮] 亚bin基线净改善门限(dB, 2026-11-07): J1 在 1° 网格
                                        %   取峰, 双原子连续搜索把"亚bin主瓣恢复"偷算进改善量 —— 剥掉
                                        %   该偏置后(dmlImprDbS): 单目标假改善实测 0.30~0.47(残差吸收
                                        %   上限), 真等强对 3.03 —— 此门关掉近凝视强目标的幻影第二点。
                                        %   强弱对真配置不在可行解中(alt 实测 imprS 0.35~0.76), 该 regime
                                        %   维持融合输出
    kapMaxDml  = 16;                    % [旋钮] Gram 条件数上限(=原子κ=4, 与下方 pinv 封顶同口径):
                                        %   超此双原子近共线, 权重解算不可信(塌缩对实测 κ_G=34~45)
    dmlStepDeg = 0.05;                  % [旋钮] 交替坐标细分步进 (deg): ±1° 内按此步进搜极小
    nSweep     = 3;                     % [旋钮] 交替坐标细分轮数
    nStart     = 20;                    % [旋钮] 残差匹配先取的峰数(再按 ±dmlSpanDeg 筛选)

    % ---------- DML 处理门限(降低虚假点: 只对近距离/小角度点做 DML) ----------
    aziLimDeg = 90;                     % [旋钮] 方位角门限: |DBF拟合方位角| 超过则跳过【2026-11-04 20→90/25→300: 验证点不因范围被删】
    xLimM     = 300;                    % [旋钮] x 位置门限 (m): |x| 超过则跳过

    nTgtAll = numel(tgtOut);            % 输入点云总数(含 OMP 第二目标)
    nDML2nd = 0;  nGate = 0;  nDone = 0;
    seedFound = false;  seed = [];
    if nargin < 8, seedD = []; end
    if nargin < 7, seedR = []; end
    if nargin < 9 || isempty(drawFigs), drawFigs = true; end   % false = 批量扫描不画图
    if nargin < 10, accumCfg = []; end                          % 遗留②多 cell 数据源
    if nargin < 11, mcOn = false; end                           % 生产权重切换到多 cell 版
    if nargin < 12, robExp = false; end                         % 远旁瓣自投影实验(诊断, 不进生产)
    haveSeed = ~isempty(seedR) && ~isempty(seedD);
    haveMc   = ~isempty(accumCfg) && isfield(accumCfg, 'chanSnapshot');

    % 统一新增标记字段(先建字段, 再用 deal 对【全部元素】填默认值; 必须在循环前 ——
    % 循环内要读 isOMP2nd。注意: 只写 tgtOut(1).字段 = 值时, 其余元素该字段是空 [],
    % 空值会在 [tgtOut.字段] 拼接时丢元素, 云图三色圈掩码就会错位)
    tgtOut(1).isDML2nd = 0;    tgtOut(1).dmlGate = 1;   tgtOut(1).dmlPass = 0;
    tgtOut(1).dmlAng1  = NaN;  tgtOut(1).dmlAng2 = NaN;
    tgtOut(1).dmlW1Db  = NaN;  tgtOut(1).dmlW2Db = NaN;
    tgtOut(1).dmlW1    = [];   tgtOut(1).dmlW2  = [];            % 线性幅度(RCS 推算用)
    tgtOut(1).dmlKappa = NaN;  tgtOut(1).dmlW1DbU = NaN;  tgtOut(1).dmlW2DbU = NaN;
    tgtOut(1).dmlSigW1 = NaN;  tgtOut(1).dmlSigW2 = NaN;         % LS 权重噪声σ(分离质量)
    tgtOut(1).dmlImprDb = NaN;                                    % 单→双原子残差改善(dB)
    tgtOut(1).dmlImprDbS = NaN;                                   % 亚bin基线改善(剥网格量化偏置)
    tgtOut(1).dmlResPm  = NaN;                                    % 拟合后残差匹配峰/中位比(dB, 真假分离判别)
    tgtOut(1).dmlAltAng = [NaN NaN];                              % 可行分离备选解角度(sep/κ 达标中 J 最优)
    tgtOut(1).dmlAltWDb = [NaN NaN];                              % 备选解权重(dB)
    tgtOut(1).dmlAltImprS = NaN;                                  % 备选解亚bin基线改善(dB)
    tgtOut(1).dmlRobAng = [NaN NaN];                              % 远旁瓣自投影实验: robust 解角度
    tgtOut(1).dmlRobImprS = NaN;                                  % robust 解亚bin基线改善(dB)
    tgtOut(1).dmlRobAltAng = [NaN NaN];                           % robust 可行备选角度
    tgtOut(1).dmlRobAltImprS = NaN;                               % robust 可行备选净改善
    tgtOut(1).dmlW1DbMc = NaN;  tgtOut(1).dmlW2DbMc = NaN;        % 多 cell 合成权重(dB, 遗留②)
    tgtOut(1).dmlW1DbSc = NaN;  tgtOut(1).dmlW2DbSc = NaN;        % 单 cell 权重备份(mcOn 时)
    tgtOut(1).dmlKeff   = NaN;                                    % 多 cell 等效快照数 Σ|h|²/|h主|²
    [tgtOut.isDML2nd] = deal(0);
    [tgtOut.dmlGate]  = deal(1);
    [tgtOut.dmlPass]  = deal(0);
    [tgtOut.dmlAng1]  = deal(NaN);  [tgtOut.dmlAng2] = deal(NaN);
    [tgtOut.dmlW1Db]  = deal(NaN);  [tgtOut.dmlW2Db] = deal(NaN);
    [tgtOut.dmlW1]    = deal([]);   [tgtOut.dmlW2]   = deal([]);
    [tgtOut.dmlKappa] = deal(NaN);  [tgtOut.dmlW1DbU] = deal(NaN);
    [tgtOut.dmlW2DbU] = deal(NaN);  [tgtOut.dmlSigW1] = deal(NaN);
    [tgtOut.dmlSigW2] = deal(NaN);  [tgtOut.dmlImprDb] = deal(NaN);
    [tgtOut.dmlImprDbS] = deal(NaN); [tgtOut.dmlResPm] = deal(NaN);
    [tgtOut.dmlW1DbMc] = deal(NaN);
    [tgtOut.dmlAltAng] = deal([NaN NaN]);  [tgtOut.dmlAltWDb] = deal([NaN NaN]);
    [tgtOut.dmlAltImprS] = deal(NaN);
    [tgtOut.dmlRobAng] = deal([NaN NaN]);  [tgtOut.dmlRobImprS] = deal(NaN);
    [tgtOut.dmlRobAltAng] = deal([NaN NaN]);  [tgtOut.dmlRobAltImprS] = deal(NaN);
    [tgtOut.dmlW2DbMc] = deal(NaN); [tgtOut.dmlW1DbSc] = deal(NaN);
    [tgtOut.dmlW2DbSc] = deal(NaN); [tgtOut.dmlKeff]   = deal(NaN);

    % ---------- 逐点 DML(跳过 OMP 第二目标; 门限内的点才做) ----------
    for t = 1:nTgtAll
        if tgtOut(t).isOMP2nd == 1
            tgtOut(t).dmlGate = 0;                       % OMP 第二目标不做 DML
            continue;
        end
        if tgtOut(t).ompHas2nd == 1
            tgtOut(t).dmlGate = 0;                       % OMP 已分辨的母点也不做 DML(2026-11-04):
            continue;                                    % 同一对被两机制重复分辨 → 图上同一目标
        end                                              % 同时叠黑圈(OMP第二)+红/青圈(DML), 标注语义失效
        azi1 = double(tgtOut(t).azimuth);                % DBF 拟合方位角(该点第一目标)
        rngM = double(tgtOut(t).rangIdx) * DeltaR;       % 距离 (m)
        eleM = double(tgtOut(t).elevation);
        xPosT = rngM * cosd(eleM) * sind(azi1);          % x 位置 (m)

        % ---- 门限检查: |方位| ≤ 90° 且 |x| ≤ 300 m; 超出 → 保留 DOA 原结果 ----
        if abs(azi1) > aziLimDeg || abs(xPosT) > xLimM
            tgtOut(t).dmlGate = 0;
            nGate = nGate + 1;
            continue;
        end

        x  = double(detCfg.doaChanIn(t, :)).';
        x  = x(chanAzi);                                 % 方位维 36 通道
        C0 = x' * x;                                     % 快照总功率
        b  = A_dml' * x;                                 % 各角度匹配向量 [M×1]
        nDone = nDone + 1;

        % ---- K=1: 最佳单角度(残差 = C0 − |b(p)|²/36, 原子单位模) ----
        [J1s, p1] = max(abs(b).^2 / 36);
        J1 = C0 - J1s;
        thS = aziGrid(p1);                               % 强目标角度(所有起点的 θ1)
        % 亚bin单原子基线角(2026-11-07 残差解剖: J1 网格量化的改善量偏置, 见写回段 dmlImprDbS)
        if p1 > 1 && p1 < numel(b)
            bbS = abs(b(p1-1:p1+1));
            d1s = 0.5*(bbS(1) - bbS(3)) / (bbS(1) - 2*bbS(2) + bbS(3) + eps);
            d1s = max(min(d1s, 0.5), -0.5);
        else
            d1s = 0;
        end
        ang1s = aziGrid(p1) + d1s;

        % ---- 残差匹配曲线 + 候选峰(第二角度起点, 限于 DBF 主峰 ±dmlSpanDeg) ----
        a1ls = (A_dml(:, p1)' * x) / (A_dml(:, p1)' * A_dml(:, p1));
        res  = x - A_dml(:, p1) * a1ls;
        c2   = A_dml' * res;
        c2Db = 20*log10(abs(c2) + eps);
        [~, pkIall] = searchPeaksPwr_pipe(c2Db, nStart);       % 多取一些峰再按范围筛选
        pkI = pkIall(abs(aziGrid(pkIall) - thS) <= dmlSpanDeg);% 限于 DBF 主峰 ±dmlSpanDeg
        if isempty(pkI), pkI = p1; end                         % 兜底: 无合格峰则用强目标自身
        nCand = numel(pkI);                              % 实际可用峰数(可能少于 nStart)

        % ---- 多起点交替坐标细分: 逻辑抽至 dmlSearch2(2026-11-07 重构, 生产行为不变; 鲁棒实验共用) ----
        % 【2026-11-06 反塌缩约束实验结论: 已回退】强制 |th1−th2|≥0.7° 的投影搜索会让第二原子
        % 被迫正交(κ_G≈1.2), 转而吸收强目标旁瓣的导向模型误差(实测单目标假改善 3.5~4.9dB,
        % 骗过 3dB 门限, C1/C2 对照 2~3 个假通过); 无约束搜索下第二原子自然滑向塌缩(改善<3dB,
        % 由判定门 sepMinDml/κ_G 拒绝)。真强弱对的差分贡献低于模型误差地板, 搜索层无解——
        % 根治在导向模型精度(标定残差进估计)与多 cell, 见报告 A.3⑤。
        [bestJ, bTh1, bTh2, bW, bSt, altJ, altTh1, altTh2, altW] = ...
            dmlSearch2(x, fitFun, thS, aziGrid(pkI), dmlStepDeg, nSweep, C0, sepMinDml, kapMaxDml, aziGrid(1), aziGrid(end));

        % ---- 角度升序整理(权重与升序角度对齐) ----
        ang1 = min(bTh1, bTh2);  ang2 = max(bTh1, bTh2);
        ampW = bW;
        if bTh1 > bTh2, ampW = ampW([2 1]); end
        % 近共线限增益(2026-11-04, 替代能量守恒缩放——LS 重构功率恒≤实测, 缩放恒为 1 无效):
        % 膨胀在"分裂"里(|W1|²+|W2|² 靠交叉对消压回重构), 根因是 (A'A)^-1 沿共线方向的特征值
        % 发散。等效岭正则: 最终解用奇异值截断的 pinv(条件数封顶 κ=4, ≈0.7° 间隔等效),
        % 共线方向的放大受限, 幅值退向"总功率按原子交叠分配"的保守解 —— 母点 +2.3dB 膨胀的修法
        A2s  = [fitFun(ang1), fitFun(ang2)];
        ampW = pinv(A2s, max(svd(A2s)) / 4) * x;

        % ---- 分离质量度量(2026-11-06, 先于判定门): Gram 条件数 κ_G + 未封顶 LS 权重 + 权重噪声σ ----
        % κ_G=κ(原子)², 下方 pinv 的 κ=4 封顶 ⇔ κ_G=16; σn²=Jf/(N-2) → cov(W)=σn²·inv(A'A)
        Gk   = A2s' * A2s;   svk = svd(Gk);
        sn2k = max(bestJ, eps) / (numel(x) - 2);
        Cwk  = sn2k * inv2x2(Gk);
        Wf   = inv2x2(Gk) * (A2s' * x);
        kapG = svk(1) / max(svk(2), eps);

        % ---- 第二目标判定(残差改善 + 亚bin净改善 + 间隔窗 + 条件数上限) ----
        imprDb = 10 * log10(max(J1, eps) / max(bestJ, eps));
        % 亚bin基线净改善: 剥掉 J1 网格量化偏置(格点抛物线插值角 + steerFitAt 重投影)
        a1sFit = fitFun(ang1s);
        J1sRes = C0 - abs(a1sFit' * x)^2 / real(a1sFit' * a1sFit);
        imprDbS = 10 * log10(max(J1sRes, eps) / max(bestJ, eps));
        dmlPass = (imprDb > thrDmlDb) && (imprDbS > thrImprSDml) && ...
                  (abs(ang1 - ang2) > sepMinDml) && ...
                  (abs(ang1 - ang2) <= sepLimDeg) && (kapG <= kapMaxDml);

        % ---- 写回逐点结果 ----
        tgtOut(t).dmlAng1 = ang1;  tgtOut(t).dmlAng2 = ang2;
        tgtOut(t).dmlW1Db = 20*log10(abs(ampW(1)) + eps);
        tgtOut(t).dmlW2Db = 20*log10(abs(ampW(2)) + eps);
        tgtOut(t).dmlW1   = ampW(1);  tgtOut(t).dmlW2 = ampW(2);   % 线性幅度(RCS 推算用)
        tgtOut(t).dmlKappa = kapG;
        tgtOut(t).dmlImprDb = imprDb;
        tgtOut(t).dmlImprDbS = imprDbS;
        % 可行分离备选解的净改善(判别实验): 备选解角度/权重/亚bin基线改善
        tgtOut(t).dmlAltAng = [min(altTh1, altTh2), max(altTh1, altTh2)];
        A2alt = [fitFun(tgtOut(t).dmlAltAng(1)), fitFun(tgtOut(t).dmlAltAng(2))];
        Walt  = pinv(A2alt, max(svd(A2alt)) / 4) * x;
        tgtOut(t).dmlAltWDb  = [20*log10(abs(Walt(1)) + eps), 20*log10(abs(Walt(2)) + eps)];
        tgtOut(t).dmlAltImprS = 10 * log10(max(J1sRes, eps) / max(altJ, eps));
        % ---- 残差结构检验(真假分离判别, 2026-11-06): 拟合后残差匹配曲线的峰/中位比 ----
        % 真对: 两原子吸收全部目标能量, 残差≈噪声 → 峰/中位比≈噪声起伏(低);
        % 假分离(第二原子吃裙边/旁瓣模型误差): 残差仍含未解释结构 → 峰/中位比显著高。
        % C1/C2 单目标对照实测标定门限(假 2 对/场), D01 真对为低侧参照。
        rRes  = x - A2s * ampW;
        mRes  = abs(A_dml' * rRes);
        tgtOut(t).dmlResPm = 20*log10(max(mRes) / (median(mRes) + eps));
        tgtOut(t).dmlW1DbU = 20*log10(abs(Wf(1)) + eps);
        tgtOut(t).dmlW2DbU = 20*log10(abs(Wf(2)) + eps);
        tgtOut(t).dmlSigW1 = sqrt(max(real(Cwk(1, 1)), 0));
        tgtOut(t).dmlSigW2 = sqrt(max(real(Cwk(2, 2)), 0));
        tgtOut(t).dmlPass = dmlPass;
        if dmlPass && tgtOut(t).splitTag == 0, tgtOut(t).splitTag = 3; end   % DML 点1(母点, 未被 OMP 标过)
        tgtOut(t).dmlSepDeg = abs(ang1 - ang2);            % DML 对角间隔(deg, 母点存, 追加点继承)

        % ---- 多 cell 相参合成(遗留②, 2026-11-07): same-RD 对在邻 cell 复增益完全相同 ----
        % 两目标共享 cell 复增益 w_c → 邻域快照与主快照平行(秩1); 主成分方向 u 的噪声
        % 降 √Keff, κ 受限的双原子幅度分裂方差 ÷Keff。幅度按主 cell 投影归一, 与单 cell
        % 口径一致(不改变角度/判定, 只提权重精度)。邻域 = 多普勒 sb±1(须同 DDM 阶, 否则
        % 该 cell 抓的是别阶目标)+ 距离 r±1; 异常 cell(PCA 残差>50%)剔除后重拟一遍。
        if haveMc
            rIdx0 = double(tgtOut(t).rangIdx);
            sb0   = mod(double(tgtOut(t).dopplerIdx) - 1, dims.subBins) + 1;
            ord0  = double(accumCfg.ddmOrdMap(rIdx0, sb0));
            candM = [0 0; -1 0; 1 0; 0 -1; 0 1];          % 主 + 距离±1 + 多普勒±1
            Xc = [];  idxM = 0;
            for qM = 1:size(candM, 1)
                rM = rIdx0 + candM(qM, 1);  sM = sb0 + candM(qM, 2);
                if rM < 1 || rM > dims.rngUse || sM < 1 || sM > dims.subBins, continue; end
                if double(accumCfg.ddmOrdMap(rM, sM)) ~= ord0, continue; end
                Xc = [Xc, double(accumCfg.chanSnapshot(chanAzi, rM, sM))]; %#ok<AGROW>
                if candM(qM, 1) == 0 && candM(qM, 2) == 0, idxM = size(Xc, 2); end
            end
            if size(Xc, 2) >= 2 && idxM > 0
                keep = true(1, size(Xc, 2));
                for itM = 1:2                              % 两轮: 第二轮用剔除后重拟
                    [Uc, ~, ~] = svd(Xc(:, keep), 0);
                    uM = Uc(:, 1);
                    hM = uM' * Xc;                         % 全部 cell 投影(含剔除的, 供 Keff 记录)
                    for qM = find(keep)
                        if qM == idxM, continue; end
                        resM = norm(Xc(:, qM) - uM * hM(qM)) / (norm(Xc(:, qM)) + eps);
                        if resM > 0.5, keep(qM) = false; end
                    end
                end
                hMain = hM(idxM);
                Wmc   = pinv(A2s, max(svd(A2s)) / 4) * (uM * hMain);
                tgtOut(t).dmlW1DbMc = 20*log10(abs(Wmc(1)) + eps);
                tgtOut(t).dmlW2DbMc = 20*log10(abs(Wmc(2)) + eps);
                tgtOut(t).dmlKeff   = sum(abs(hM(keep)).^2) / (abs(hMain)^2 + eps);
                if mcOn && dmlPass
                    tgtOut(t).dmlW1DbSc = tgtOut(t).dmlW1Db;   % 单 cell 值备份(A/B 用)
                    tgtOut(t).dmlW2DbSc = tgtOut(t).dmlW2Db;
                    tgtOut(t).dmlW1Db = tgtOut(t).dmlW1DbMc;   % 生产切多 cell 版
                    tgtOut(t).dmlW2Db = tgtOut(t).dmlW2DbMc;
                    tgtOut(t).dmlW1   = Wmc(1);  tgtOut(t).dmlW2 = Wmc(2);
                end
            end
        end

        % ---- 远旁瓣残差自投影实验(2026-11-07, robExp=true 时计算, 不进生产链) ----
        % 生产口径自查: ĥ 由【当前点自身残差】在线估计, 只用点数据+既有导向表+设计常数,
        % 无真值/无外部残差表 —— 与 cCal/thrDmlDb 同级的生产合法构造。
        % 原理: 导向模型误差(实测 −10.7dB 角结构)集中在 |θ−θ₁|>5° 远旁瓣, 而真第二目标的差分
        % 在主瓣(≤1.2°) —— 用残差远旁瓣峰重构干扰子空间 Ĥ 并投影后重跑搜索: 第二原子吃不到
        % 模型误差(假改善消失), 真差分与 Ĥ 正交(应保留)。
        if robExp
            a1r  = fitFun(ang1s);
            al1r = (a1r' * x) / real(a1r' * a1r);
            r1   = x - a1r * al1r;
            cRf  = A_dml' * r1;
            cF   = cRf;  cF(abs(aziGrid - ang1s) <= 5) = 0;     % 只留远旁瓣(主瓣保护 5°)
            [~, fk] = maxk(abs(cF), 3);
            Hh = [];
            for q = 1:numel(fk)
                aF = fitFun(aziGrid(fk(q)));
                if ~isempty(Hh), aF = aF - Hh * (Hh' * aF); end % Gram-Schmidt 正交化
                nF = norm(aF);
                if nF > 1e-3, Hh = [Hh, aF / nF]; end %#ok<AGROW>
            end
            if ~isempty(Hh)
                Pr   = eye(numel(x)) - Hh * Hh';
                xP   = Pr * x;
                fitP = @(th) Pr * fitFun(th);                    % 原子同步投影(投影空间一致 LS)
                C0p  = real(xP' * xP);
                a1p  = fitP(ang1s);
                J1p  = C0p - abs(a1p' * xP)^2 / real(a1p' * a1p);
                [rJ, rT1, rT2, ~, ~, rAJ, rAT1, rAT2] = ...
                    dmlSearch2(xP, fitP, ang1s, aziGrid(pkI), dmlStepDeg, nSweep, C0p, sepMinDml, kapMaxDml, aziGrid(1), aziGrid(end));
                tgtOut(t).dmlRobAng      = [min(rT1, rT2), max(rT1, rT2)];
                tgtOut(t).dmlRobImprS    = 10 * log10(max(J1p, eps) / max(rJ, eps));
                tgtOut(t).dmlRobAltAng   = [min(rAT1, rAT2), max(rAT1, rAT2)];
                tgtOut(t).dmlRobAltImprS = 10 * log10(max(J1p, eps) / max(rAJ, eps));
            end
        end

        % ---- 通过: 第一目标方位替换为 DML 第一角度, 追加第二目标进点云 ----
        if dmlPass
            nDML2nd = nDML2nd + 1;
            tgtOut(t).azimuth = ang1;                        % 第一目标: DBF角 → DML角
            newT          = tgtOut(t);                       % 继承同距同速/俯仰等字段
            newT.ID       = max([tgtOut.ID]) + 1;
            newT.azimuth  = ang2;                            % 第二目标方位
            newT.isDML2nd = 1;
                newT.splitTag = 4;                            % DML 点2
                newT.dmlSepDeg = abs(ang1 - ang2);            % 对角间隔(同母点)
            newT.dmlPass  = 0;                               % 追加点不是"被替换的母点", 清标记(否则计数翻倍)
            tgtOut(end+1) = newT;                            % #ok<AGROW> 追加段在云图用 is2 标记
        end

        % ---- 指定点(R/D 序号, 容差 ±1 bin) → 保存 DML 过程数据 ----
        if haveSeed && ~seedFound && ...
                abs(tgtOut(t).rangIdx - seedR) <= 1 && abs(tgtOut(t).dopplerIdx - seedD) <= 1
            seedFound = true;
            seed = struct('t', t, 'rngIdx', tgtOut(t).rangIdx, 'dopIdx', tgtOut(t).dopplerIdx, ...
                          'azi1', azi1, 'ang1', ang1, 'ang2', ang2, 'pass', dmlPass, ...
                          'imprDb', imprDb, 'nCand', nCand, 'bSt', bSt, ...
                          'bDb', 20*log10(abs(b)/max(abs(b)) + eps), 'c2Db', c2Db, ...
                          'pkI', pkI, 'x', x, 'W', ampW);
        end
    end
    if haveSeed && ~seedFound
        fprintf('DML 提示: RD(%g, %g) 处未找到确认点, 过程图跳过\n', seedR, seedD);
    end

    % ---------- 图窗公共尺寸(与 OMP/DOA 模块一致) ----------
    screenSize = get(0, 'ScreenSize');
    figW = screenSize(3) * 0.98;  figH = screenSize(4) * 0.80;
    figPos = [(screenSize(3)-figW)/2, (screenSize(4)-figH)/2, figW, figH];

    % ---------- 指定点的 DML 过程图(图24, 4 格; 未指定/未找到则为空) ----------
    figSeed = [];
    if seedFound && drawFigs
        vTxt = '不通过';  if seed.pass, vTxt = '通过'; end
        figSeed = figure('Position', figPos, 'Name', 'fig31_DML_process');
        tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

        nexttile;   % DBF 谱(单峰难分) + DML 两角度竖线
        plot(aziGrid, seed.bDb, 'b-'); hold on;
        xline(seed.ang1, 'b--', sprintf('%.2f°', seed.ang1));
        xline(seed.ang2, 'r--', sprintf('%.2f°', seed.ang2));
        grid on; xlim([min(aziGrid), max(aziGrid)]); ylim([-40 2]);
        xlabel('方位角 (deg)'); ylabel('归一化幅度 (dB)');
        title(sprintf('常规 DBF 谱 (DML 反演: %.2f° / %.2f°)', seed.ang1, seed.ang2));

        nexttile;   % 残差匹配曲线 + 全部候选起点 + 胜出起点
        plot(aziGrid, seed.c2Db - max(seed.c2Db), 'r-'); hold on;
        plot(aziGrid(seed.pkI), seed.c2Db(seed.pkI) - max(seed.c2Db), 'co', 'MarkerSize', 6);
        plot(aziGrid(seed.pkI(seed.bSt)), seed.c2Db(seed.pkI(seed.bSt)) - max(seed.c2Db), ...
            'gv', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
        yline(median(seed.c2Db) - max(seed.c2Db), 'k:');
        grid on; ylim([-40 10]);
        xlabel('方位角 (deg)'); ylabel('归一化幅度 (dB)');
        legend({'残差匹配曲线', '候选起点', '胜出起点'}, 'Location', 'best');
        title(sprintf('残差匹配: DBF主峰 ±%g°内 %d 个候选起点, 胜出 %.2f°', ...
            dmlSpanDeg, seed.nCand, aziGrid(seed.pkI(seed.bSt))));

        nexttile;   % DML 二维代价面(胜出解附近 ±2°)
        xS  = seed.x;  C0S = real(xS' * xS);
        fA  = seed.ang1 + (-2 : 0.1 : 2);  fB = seed.ang2 + (-2 : 0.1 : 2);
        Jm  = zeros(numel(fB), numel(fA));
        for ib = 1:numel(fB)
            for ia = 1:numel(fA)
                A2  = [fitFun(fA(ia)), fitFun(fB(ib))];
                aS  = inv2x2(A2' * A2) * (A2' * xS);
                Jm(ib, ia) = C0S - real(xS' * A2 * aS);
            end
        end
        JmPos = max(real(Jm), 0);            % 负值(θ1=θ2 退化线 inv2x2 近奇异)截 0, real 防复数
        JmDb  = 10*log10(JmPos / (max([JmPos(:); 0]) + eps) + eps);  % 分母追加 0: Jm 空阵时 max([]) 为 [] 会维度崩
        JmDb(~isfinite(JmDb)) = 0;
        imagesc(fA, fB, JmDb); axis xy; colorbar;
        hold on; plot(seed.ang1, seed.ang2, 'wo', 'MarkerSize', 10);
        plot(xlim, xlim, 'w:', 'LineWidth', 1);          % θ1=θ2 退化线(两原子重合, 恒为脊)
        xlabel('θ1 (deg)'); ylabel('θ2 (deg)');
        title(sprintf('DML 二维代价面 (谷底 %.2f° / %.2f°, 白虚线=θ1=θ2退化线)', seed.ang1, seed.ang2));

        nexttile;   % 重建检验
        xRef = max(abs(xS)) + eps;
        atom1R = steerFitAt(A_dml, aziGrid, seed.ang1);
        atom2R = steerFitAt(A_dml, aziGrid, seed.ang2);
        plot(20*log10(abs(xS)/xRef + eps), 'b-o', 'MarkerSize', 4); hold on;
        plot(20*log10(abs([atom1R, atom2R] * seed.W)/xRef + eps), 'r-*', 'MarkerSize', 4);
        grid on; ylim([-60 5]);
        xlabel('通道 (最长行 36)'); ylabel('相对快照峰值 (dB)');
        legend({'快照(归一化)', 'DML 双原子重建(同一参考)'}, 'Location', 'best');
        title('重建检验: 快照 vs 双原子重建');
        sgtitle(sprintf('DML 多起点超分辨: RD(%g,%g) → %.2f° / %.2f° | 残差改善 %.1f dB | 判定: %s', ...
            seedR, seedD, seed.ang1, seed.ang2, seed.imprDb, vTxt));

        fprintf('DML 指定点      : RD(%g,%g) → %.2f° / %.2f° | 残差改善 %.1f dB | 判定: %s\n', ...
            seedR, seedD, seed.ang1, seed.ang2, seed.imprDb, vTxt);
        fprintf('DML 幅度信息   : 两目标 LS 幅度 %.3g / %.3g (积累快照量纲, 供RCS推算) | 相对差 %.1f dB\n', ...
            abs(seed.W(1)), abs(seed.W(2)), 20*log10(abs(seed.W(1))/abs(seed.W(2))));
    end

    if ~drawFigs
        figCloud = [];                                   % 批量扫描: 不画图
    end
    if drawFigs
    % ---------- 更新点云的 6 子图 DBF 点迹图(同 fig28_DBF 口径, 含三种角度圈) ----------
    rngEst = cell2mat(arrayfun(@(q) double(q.range), tgtOut, 'UniformOutput', false));
    dopEst = cell2mat(arrayfun(@(q) double(q.dopplerIdx) - 1 - double(q.doppler_alpha), tgtOut, 'UniformOutput', false)); % 0基小数bin
    aziEst = cell2mat(arrayfun(@(q) double(q.azimuth),  tgtOut, 'UniformOutput', false));
    eleEst = cell2mat(arrayfun(@(q) double(q.elevation), tgtOut, 'UniformOutput', false));
    zEst   = rngEst .* sind(eleEst);
    xPos   = rngEst .* cosd(eleEst) .* sind(aziEst);
    yPos   = rngEst .* cosd(eleEst) .* cosd(aziEst);
    % 四种圈: 品红圈=OMP目标1(母点) / 黑圈=OMP目标2 / 红圈=DML目标1(DBF角已被DML角替换) / 青圈=DML目标2
    isOMP1 = logical([tgtOut.splitTag] == 1);                    % OMP 目标1(母点, 通常较强)
    isOMP2 = logical([tgtOut.isOMP2nd]);                         % OMP 第二目标
    isD1   = false(1, numel(tgtOut));
    isD1(1:nTgtAll) = logical([tgtOut(1:nTgtAll).dmlPass]);      % DML 第一目标(原点段, 角度已替换)
    isD2   = logical([tgtOut.isDML2nd]);                         % DML 第二目标(追加段)
    rMax10 = ceil(dims.rngUse * 3e8/(2*dims.B) / 10) * 10;
    sym10  = @(v) max(ceil(max(abs(v)) / 10) * 10, 10);   % 全 0 时保底 ±10, 防 xlim([0,0]) 非法

    figCloud = figure('Position', figPos, 'Name', 'fig32_DBF_OMP_DML');
    tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    nexttile;
    hP = scatter(xPos, yPos, 24, zEst, 'filled'); hold on;
    hM = scatter(xPos(isOMP1), yPos(isOMP1), 40, 'mo', 'LineWidth', 1.5);
    hR = scatter(xPos(isD1), yPos(isD1), 40, 'ro', 'LineWidth', 1.5);
    hC = scatter(xPos(isD2), yPos(isD2), 40, 'co', 'LineWidth', 1.5);
    hK = scatter(xPos(isOMP2), yPos(isOMP2), 40, 'ko', 'LineWidth', 1.5); grid on; colorbar;
    % 句柄显式传入: 无句柄时 legend 按对象创建顺序贴标签, 底图彩点(第一个对象)会被
    % 误标为"OMP第二目标(黑圈)", 黑圈(最后一个对象)反无标签 —— 整体错一位
    legend([hP hM hK hR hC], {'DBF目标(独立,颜色=z)', 'OMP目标1(母点,品红圈)', 'OMP目标2(黑圈)', ...
        'DML目标1(母点,红圈)', 'DML目标2(青圈)'}, 'Location', 'best');
    xlim([-sym10(xPos), sym10(xPos)]); ylim([0, rMax10]);
    xlabel('y 横向 (m)'); ylabel('x 纵向 (m)'); title('横向-纵向 位置(俯视: 横向在横轴)');
    nexttile;
    scatter(aziEst, rngEst, 24, zEst, 'filled'); hold on;
    scatter(aziEst(isD1), rngEst(isD1), 40, 'ro', 'LineWidth', 1.5); hold on;
    scatter(aziEst(isD2), rngEst(isD2), 40, 'co', 'LineWidth', 1.5); hold on;
    scatter(aziEst(isOMP2), rngEst(isOMP2), 40, 'ko', 'LineWidth', 1.5); grid on; colorbar;
    xlim([-sym10(aziEst), sym10(aziEst)]); ylim([0, rMax10]);
    xlabel('方位角 (deg)'); ylabel('距离 (m)'); title('方位角-距离');
    nexttile;
    scatter(dopEst, rngEst/dims.DeltaR, 24, zEst, 'filled'); grid on; colorbar;
    xlim([0, dims.dopBinNum-1]); ylim([0, dims.rngUse]);
    xlabel('多普勒 bin'); ylabel('距离 bin'); title('多普勒-距离 (满bin)');
    nexttile;
    scatter(zEst, yPos, 24, zEst, 'filled'); hold on;
    scatter(zEst(isD1), yPos(isD1), 40, 'ro', 'LineWidth', 1.5); hold on;
    scatter(zEst(isD2), yPos(isD2), 40, 'co', 'LineWidth', 1.5); hold on;
    scatter(zEst(isOMP2), yPos(isOMP2), 40, 'ko', 'LineWidth', 1.5); grid on; colorbar;
    xlim([-sym10(zEst), sym10(zEst)]); ylim([0, rMax10]);
    xlabel('z (m)'); ylabel('x 纵向 (m)'); title('高度-纵向 位置(侧视)');
    nexttile;
    scatter(eleEst, rngEst, 24, zEst, 'filled'); hold on;
    scatter(eleEst(isD1), rngEst(isD1), 40, 'ro', 'LineWidth', 1.5); hold on;
    scatter(eleEst(isD2), rngEst(isD2), 40, 'co', 'LineWidth', 1.5); hold on;
    scatter(eleEst(isOMP2), rngEst(isOMP2), 40, 'ko', 'LineWidth', 1.5); grid on; colorbar;
    xlim([-sym10(eleEst), sym10(eleEst)]); ylim([0, rMax10]);
    xlabel('俯仰角 (deg)'); ylabel('距离 (m)'); title('俯仰角-距离');
    nexttile;
    scatter3(xPos, yPos, zEst, 24, zEst, 'filled'); hold on;
    scatter3(xPos(isD1), yPos(isD1), zEst(isD1), 40, 'ro', 'LineWidth', 1.5); hold on;
    scatter3(xPos(isD2), yPos(isD2), zEst(isD2), 40, 'co', 'LineWidth', 1.5); hold on;
    scatter3(xPos(isOMP2), yPos(isOMP2), zEst(isOMP2), 40, 'ko', 'LineWidth', 1.5); grid on;
    xlabel('y 横向 (m)'); ylabel('x 纵向 (m)'); zlabel('z (m)'); view(3);
    title('点云 xyz (可旋转)');
    tagCnt = accumarray(double([tgtOut.splitTag]).' + 1, 1, [5 1]);   % 0~4 → 五类计数
    sgtitle(sprintf(['DML 后点云: 独立点 %d | OMP第一目标(母点,品红圈) %d | OMP第二目标(黑圈) %d | ' ...
        'DML第一目标(母点,红圈,DBF角已替换) %d | DML第二目标(青圈) %d'], tagCnt));
    end                                             % if drawFigs(云图段)

    dmlCfg = struct('figSeed', figSeed, 'figCloud', figCloud, 'nDML2nd', nDML2nd, ...
                    'nGate', nGate, 'nDone', nDone, 'seedFound', seedFound);
end

function [bestJ, bTh1, bTh2, bW, bSt, altJ, altTh1, altTh2, altW] = ...
        dmlSearch2(x, fitFun, thS, candAngs, dmlStepDeg, nSweep, C0, sepMinDml, kapMaxDml, thLo, thHi)
%DMLSEARCH2 多起点交替坐标双原子搜索(2026-11-07 从主链抽取, 生产行为逐句不变; 鲁棒实验共用)
%   返回: 全局最优解(best*) + 可行分离备选解(alt*: sep>sepMinDml 且 Gramκ≤kapMaxDml 中 J 最小)
altJ = inf;  altTh1 = NaN;  altTh2 = NaN;  altW = [NaN; NaN];
bestJ = inf;  bTh1 = [];  bTh2 = [];  bW = [];  bSt = 0;
for st = 1:numel(candAngs)
    th1 = thS;  th2 = candAngs(st);                   % 起点: 强目标 + 第 st 个候选峰
    for sw = 1:nSweep
        f1 = min(max(th1 + (-1:dmlStepDeg:1), thLo), thHi);  Jf1 = zeros(size(f1));
        for j1 = 1:numel(f1)
            A2 = [fitFun(f1(j1)), fitFun(th2)];
            aS = inv2x2(A2'*A2) * (A2'*x);
            Jf1(j1) = C0 - real(x'*A2*aS);
        end
        [~, j1m] = min(Jf1);  th1 = f1(j1m);
        f2 = min(max(th2 + (-1:dmlStepDeg:1), thLo), thHi);  Jf2 = zeros(size(f2));
        for j2 = 1:numel(f2)
            A2 = [fitFun(th1), fitFun(f2(j2))];
            aS = inv2x2(A2'*A2) * (A2'*x);
            Jf2(j2) = C0 - real(x'*A2*aS);
        end
        [~, j2m] = min(Jf2);  th2 = f2(j2m);
    end
    A2f = [fitFun(th1), fitFun(th2)];
    aS  = inv2x2(A2f'*A2f) * (A2f'*x);
    Jf  = C0 - real(x'*A2f*aS);
    if Jf < bestJ
        bestJ = Jf;  bTh1 = th1;  bTh2 = th2;  bW = aS;  bSt = st;
    end
    GkS = A2f'*A2f;  svS = svd(GkS);
    if abs(th1 - th2) > sepMinDml && svS(1)/max(svS(2), eps) <= kapMaxDml && Jf < altJ
        altJ = Jf;  altTh1 = th1;  altTh2 = th2;  altW = aS;
    end
end
end
