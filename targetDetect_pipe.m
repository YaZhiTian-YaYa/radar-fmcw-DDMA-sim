function [accumCfg, detCfg, tgtOut] = targetDetect_pipe(noiseCfg, accumCfg, detCfg, dims, drawFigs)
%TARGETDETECT_PIPE 方位谱DACFAR粗检 + 相参/非相参峰值状态精判, 输出目标列表
%   (移植自 参考\rsp_pipeline_ReWrite\pipeline\targetDetect_pipe.m, 变量名与参考一致;
%    硬编码维度改由主函数 dims 传入; 末尾新增检测三子图 figDet)
%   输入: noiseCfg -- 噪底估计结果(floorPow等)
%         accumCfg -- 积累结果(nonCohMap/aziFftPowDb/ddmOrdMap/chanSnapshot等)
%         detCfg   -- 来自 SimAdcParam.mat 的 cfarCfg, 含逐距离门CFAR与峰判配置
%         dims     -- 公共维度(主函数唯一定义)
%   drawFigs: false = 只检不画(校准遍只取角反点云, 主流程两遍测量用)
%   输出: tgtOut -- 确认目标列表(距离/多普勒索引, 幅度, SNR, DDM阶, 峰态等)
if nargin == 9        % (遗留代码, 本函数仅3个输入, 恒不触发)
    dummyBuf = [];
end
if nargin < 5 || isempty(drawFigs), drawFigs = true; end   % 默认画图
    % ---------- 初始化空目标列表(逐字段为空, 检测到后逐个填入) ----------
    tgtOut = struct('ID', [], ...
                       'rangIdx', [], ...
                       'dopplerIdx', [], ...
                       'range', [], ...
                       'speed', [], ...
                       'amb_speed', [], ...
                       'azimuth', [], ...
                       'elevation', [], ...
                       'xPos', [], ...
                       'yPos', [], ...
                       'zPos', [], ...
                       'nonCohMag', [], ...
                       'nonCohSNR', [], ...
                       'cohMag', [], ...
                       'cohPeakIdx', [], ...
                       'doaMag', [], ...
                       'tx_order', [], ...
                       'azimuthList', [], ...
                       'sideLobe', [], ...
                       'range_alpha', [], ...
                       'doppler_alpha', [], ...
                       'peak_type', []);


    rngBinTotal = dims.rngUse;      % 距离bin总数(参与检测)
    dopSubNum = dims.subBins;         % 多普勒子带数(每组32个bin)
    fftLen = dims.aziFftLen;           % 方位FFT长度(此处未直接使用)
    maxTgtNum = 4000;       % 单帧最大目标数


    peakPowMap = zeros(dims.rngUse, dims.dopBinNum);          % 相参方位峰功率图(512距离×512多普勒)
    peakPosMap = zeros(dims.rngUse, dims.dopBinNum);          % 相参方位峰位置图
    detCfg.detFlagMap = zeros(dims.rngUse, dims.dopBinNum);   % CFAR粗检标记图
    ratioBuf = [];                          % 过检单元的CFAR比值收集缓存
    dacfarRngCnt = dims.rngUse;                     % DACFAR参与检测的距离门数

    %% ---------- 第一步: DACFAR粗检(在方位-多普勒二维谱上做CFAR) ----------
    if(1)   % 总开关(恒开, 便于调试时整体跳过)
        % 方位谱dB转线性幅度, 并permute成dacfar要求的 [方位 × 多普勒 × 距离];
        % 粗/细门限均取配置值×1.1加严
        [detCfg.coarseDetList, peakPosRow, peakPowRow, detCfg.detFlagMap, cfarRatioMap] = dacfar_version2_pipe(permute(10.^(accumCfg.aziFftPowDb / 20), [1,3,2]), ...
                                                                    noiseCfg.floorPow, detCfg.cfarCoarseList*1.1, detCfg.cfarThreList*1.1, detCfg.guardWin, ...
                                                                    detCfg.refWin, detCfg.cfarType, dacfarRngCnt, noiseCfg.weakChirpFlag);
        detCfg.detFlagMapCoarse = detCfg.detFlagMap;        % 【新增】粗检标记图副本(检测三子图第1张用)
        cohDetMap = peakPowRow > (noiseCfg.floorPow.*detCfg.cfarCoarseList);   % (计算后未使用)

        % DACFAR输出的是32个多普勒子带的结果,
        % 按各单元DDM阶数展开回完整512个多普勒bin: dopBin = 子带 + (阶数-1)*32
        for rBin = 1:dacfarRngCnt
            dopBin = (1:dopSubNum) + (accumCfg.ddmOrdMap(rBin, :) - 1) * dopSubNum;
            peakPowMap(rBin, dopBin) = peakPowRow(rBin, :);

            peakPosMap(rBin, dopBin) = peakPosRow(rBin, :);

            % 收集过粗检单元的CFAR信噪比(门限判据: 峰功率 > 噪底×该门粗门限)
            passColIdx = find(peakPowRow(rBin, :) > (noiseCfg.floorPow(rBin)*detCfg.cfarCoarseList(rBin)));
            ratioBuf = [ratioBuf, cfarRatioMap(rBin, passColIdx)];
        end
    end


    detCfg.ratioBuf = ratioBuf;    % CFAR比值(调试/统计用)
    detCfg.peakPowMap = peakPowMap;    % 相参峰功率图

    detCfg.peakPosMap = peakPosMap;    % 相参峰位置图


    %% ---------- 第二步: 峰值状态精判(抑制CFAR漏出的非真实峰) ----------
    % 对每个CFAR过检单元, 检查其在相参谱(peakPowMap)与非相参谱(ncPwrMap)上
    % 是否构成二维局部峰, 逐距离门按配置模式pkModeTab综合判定
    ncPwrMap = accumCfg.nonCohMap;    % 非相参积累功率图 [512距离 × 32子带]
    tgtCnt = 0;
    detCfg.outBufMap = zeros(size(detCfg.detFlagMap));
    detCfg.doaChanIn = [];            % 测角用通道快照收集(每目标一行)

    for rngIdx = 1:rngBinTotal

        if rngIdx > dacfarRngCnt
            continue;                % 超出DACFAR范围的距离门不参与(当前恒不触发)
        end

        % 逐距离门选择峰判模式表(当前dacfarRngCnt=512, else分支为遗留路径)
        if(rngIdx <= dacfarRngCnt)
            pkModeTab = detCfg.cfarPeakList;
        else
            pkModeTab = detCfg.non_cfarPeakList;
        end

        for sbIdx = 1:dopSubNum

            % 该单元在完整多普勒轴上的bin序号
            dopBin = sbIdx + (accumCfg.ddmOrdMap(rngIdx, sbIdx) - 1) * dopSubNum;
            % 在粗检列表中找到该(距离, 子带)对应的行, 回填展开后的多普勒bin
            % (粗检为空时跳过 —— 空场景/目标全部出带时 coarseDetList 为 0×0,
            %  直接索引 (:,2) 会越界报错)
            if ~isempty(detCfg.coarseDetList)
                coarseRowIdx = find(detCfg.coarseDetList(:,1) == rngIdx & detCfg.coarseDetList(:,2) == mod(dopBin-1, dopSubNum)+1);
                if size(coarseRowIdx,1) > 0
                    detCfg.coarseDetList(coarseRowIdx, 5) = dopBin;
                end
            end

            % ---- 准备邻域索引: 距离/多普勒各取前后bin, 边界截断, 多普勒循环回绕 ----
            % (注意: nxt2Rng/nxt2Dop实际取+1, 与prev2侧不对称)
            prevRng = max(1, rngIdx - 1);    nxtRng = min(rngIdx+1, dims.rngUse);   % 距离邻域上界是 rngUse(勿写 dopBinNum, 原 512==512 巧合)
            prev2Rng = max(1, prevRng - 1);    nxt2Rng = min(rngIdx+1, dims.rngUse);  % 同上

            prevDop = dopBin - 1;            if prevDop == 0, prevDop = dims.dopBinNum;end
            nxtDop = dopBin + 1;            if nxtDop > dims.dopBinNum, nxtDop = 1;end

            prev2Dop = prevDop - 1;            if prev2Dop == 0, prev2Dop = dims.dopBinNum;end
            nxt2Dop = dopBin + 1;            if nxt2Dop > dims.dopBinNum, nxt2Dop = 1;end

            if(rngIdx <= dacfarRngCnt)

                % cohPkCnt: 与4个邻域(距离±1, 多普勒±1)比较的"胜出"次数;
                % 功率占优 或 方位峰位相差>3bin(容忍谱峰漂移)均计为胜出
                cohPkCnt = (peakPowMap(rngIdx, dopBin) > peakPowMap(prevRng, dopBin) | abs(peakPosMap(rngIdx, dopBin) - peakPosMap(prevRng, dopBin)) > 3)...
                        +   (peakPowMap(rngIdx, dopBin) >= peakPowMap(nxtRng, dopBin) | abs(peakPosMap(rngIdx, dopBin) - peakPosMap(nxtRng, dopBin)) >= 3)...
                        +   (peakPowMap(rngIdx, dopBin) > peakPowMap(rngIdx, prevDop) | abs(peakPosMap(rngIdx, dopBin) - peakPosMap(rngIdx, prevDop)) > 3)...
                        +   (peakPowMap(rngIdx, dopBin) >= peakPowMap(rngIdx, nxtDop) | abs(peakPosMap(rngIdx, dopBin) - peakPosMap(rngIdx, nxtDop)) >= 3);


                % aziSameCnt: 与邻域共享同一方位峰位的次数(越大越可能是同一强目标的延伸)
                aziSameCnt = ((peakPosMap(rngIdx, dopBin) == peakPosMap(prevRng, dopBin)))...
                        +   ((peakPosMap(rngIdx, dopBin) == peakPosMap(nxtRng, dopBin)))...
                        +   ((peakPosMap(rngIdx, dopBin) == peakPosMap(rngIdx, prevDop)))...
                        +   ((peakPosMap(rngIdx, dopBin) == peakPosMap(rngIdx, nxtDop)));
                % dopPkFlg: 本距离门内多普勒维的严格局部峰(左严格小, 右不大于)
                dopPkFlg = (peakPowMap(rngIdx, dopBin) > peakPowMap(rngIdx, prevDop)) && (peakPowMap(rngIdx, dopBin) >= peakPowMap(rngIdx, nxtDop));
                % cohPkCntAll: 不含峰位容忍的严格4邻域比较(仅比功率)
                cohPkCntAll = (peakPowMap(rngIdx, dopBin) > peakPowMap(prevRng, dopBin)) ...
                        +   (peakPowMap(rngIdx, dopBin) >= peakPowMap(nxtRng, dopBin)) ...
                        +   (peakPowMap(rngIdx, dopBin) > peakPowMap(rngIdx, prevDop)) ...
                        +   (peakPowMap(rngIdx, dopBin) >= peakPowMap(rngIdx, nxtDop) );

            else
                cohPkCnt = 0;
                aziSameCnt = 0;
                cohPkCntAll = 0;
                dopPkFlg = 0;
            end
            % ncPkCnt: 非相参谱上与4邻域比较的胜出次数(多普勒子带循环回绕)
            ncPkCnt = (ncPwrMap(rngIdx, mod(dopBin-1,dopSubNum)+1) > ncPwrMap(prevRng, mod(dopBin-1,dopSubNum)+1))...
                        +   (ncPwrMap(rngIdx, mod(dopBin-1,dopSubNum)+1) >= ncPwrMap(nxtRng, mod(dopBin-1,dopSubNum)+1))...
                        +   (ncPwrMap(rngIdx, mod(dopBin-1,dopSubNum)+1) > ncPwrMap(rngIdx, mod(prevDop-1,dopSubNum)+1))...
                        +   (ncPwrMap(rngIdx, mod(dopBin-1,dopSubNum)+1) >= ncPwrMap(rngIdx, mod(nxtDop-1,dopSubNum)+1));


            detPass = 0;

            % ---- 按该距离门配置的模式(0~8)组合上述峰态标志做精判 ----
            % 模式0最严, 越往后判据组合越宽松/侧重不同
            switch(pkModeTab(rngIdx))
                case 0  %%强峰值约束: 需多普勒局部峰; 无同方位邻峰时相参≥3且非相参≥2; 有则需相参全峰且非相参≥3, 或非相参全峰
                    if dopPkFlg && ((aziSameCnt == 0 && (cohPkCnt >= 3 && ncPkCnt >= 2)) ...
                                                || (aziSameCnt > 0 && ((cohPkCnt == 4 && ncPkCnt >= 3) || ncPkCnt == 4)))
                        detPass = 1;
                    end
                case 1  % 需多普勒局部峰; 无同方位邻峰时相参≥3; 有则需非相参≥3或相参全峰
                    if dopPkFlg && ((aziSameCnt == 0 && (cohPkCnt >= 3)) ...
                                                || (aziSameCnt > 0 && ((ncPkCnt >= 3) || cohPkCnt == 4)))
                        detPass = 1;
                    end
                case 2  %% 需多普勒局部峰, 且相参或非相参≥3
                    if(dopPkFlg && (cohPkCnt >= 3 || ncPkCnt >= 3))
                        detPass = 1;
                    end
                case 3  % 不要求多普勒局部峰: 相参全峰且非相参≥3, 或非相参全峰
                    if((cohPkCnt == 4 && ncPkCnt >= 3) || ncPkCnt == 4)
                        detPass = 1;
                    end
                case 4  % 仅要求非相参全峰
                    if(ncPkCnt == 4)
                        detPass = 1;
                    end
                case 5  % 相参或非相参≥3
                    if(cohPkCnt >= 3 || ncPkCnt >= 3)
                        detPass = 1;
                    end
                case 6  % 最宽松: 相参或非相参≥2
                    if(cohPkCnt >= 2 || ncPkCnt >= 2)
                        detPass = 1;
                    end
                case 7  % 相参全峰或非相参全峰
                    if(cohPkCnt >= 4 || ncPkCnt >= 4)
                        detPass = 1;
                    end
                case 8  % 仅当存在同方位邻峰时: 非相参≥3或相参全峰(强目标拖尾场景)
                    if((aziSameCnt > 0 && ((ncPkCnt >= 3) || cohPkCnt == 4)))
                        detPass = 1;
                    end
            end
            % 峰态编码: cohPkType/ncPkType, 3=严格四邻域全峰, 4=三邻域峰, 0=其他
            if cohPkCntAll == 4
                cohPkType = 3;
            elseif cohPkCntAll == 3
                cohPkType = 4;
            else
                cohPkType = 0;
            end
            if ncPkCnt == 4
                ncPkType = 3;
            elseif ncPkCnt == 3
                ncPkType = 4;
            else
                ncPkType = 0;
            end


            % ---- 最终判决: 峰态精判通过 且 CFAR粗检通过, 才确认为目标 ----
            if detPass && detCfg.detFlagMap(rngIdx, mod(dopBin-1, dopSubNum)+1) %&& ~(isRainCurv)
                tgtCnt = tgtCnt + 1;
                % 距离维亚bin插值: 取距离方向5点功率(前2/前1/当前/后1/后2)
                rngSubOff = interpolation_fun_pipe(ncPwrMap([prev2Rng, prevRng, rngIdx, nxtRng, nxt2Rng], sbIdx));

                rngOffBin = round(rngSubOff + rngIdx);    % 亚bin校正后的距离bin(用于取多普勒插值邻域)
                if(rngOffBin > dims.dopBinNum)
                    rngOffBin = dims.dopBinNum;
                elseif(rngOffBin <= 0)
                    rngOffBin = 1;
                end
                % 多普勒维亚bin插值: 在校正后距离bin处取多普勒方向5点功率
                dopSubOff = interpolation_fun_pipe(ncPwrMap(rngOffBin, mod([prev2Dop, prevDop, sbIdx, nxtDop, nxt2Dop] - 1, dopSubNum) + 1));
                % ---- 填写目标参数 ----
                tgtOut(tgtCnt).ID = tgtCnt - 1;                          % 目标编号(从0起)
                tgtOut(tgtCnt).rangIdx = rngIdx;                         % 距离bin索引
                tgtOut(tgtCnt).dopplerIdx = dopBin;                      % 多普勒bin索引(解模糊后全bin)
                tgtOut(tgtCnt).nonCohMag = 10 * log10(accumCfg.nonCohMap(rngIdx, sbIdx));            % 非相参幅度(dB)
                tgtOut(tgtCnt).nonCohSNR = 10 * log10(accumCfg.nonCohMap(rngIdx, sbIdx) / noiseCfg.floorPow(rngIdx));  % 非相参SNR(dB)

                tgtOut(tgtCnt).cohMag = 10 * log10(abs(accumCfg.cohPeakPow(rngIdx, sbIdx)));          % 相参峰幅度(dB)
                tgtOut(tgtCnt).cohPeakIdx = peakPosMap(rngIdx, dopBin);  % 方位峰bin位置

                tgtOut(tgtCnt).tx_order = accumCfg.ddmOrdMap(rngIdx, sbIdx);   % 该目标的DDM阶数

                tgtOut(tgtCnt).range_alpha = rngSubOff;                  % 距离亚bin偏移
                tgtOut(tgtCnt).doppler_alpha = dopSubOff;                % 多普勒亚bin偏移

                % ---- 亚bin插值回填(偏移约定: alpha>0 表示峰偏向【前】侧低 bin) ----
                % 距离 = (bin1起 − 1 − range_alpha)×ΔR; 速度 = (解模糊bin1起 − 1 − doppler_alpha)×Vres
                % 折叠速度 amb_speed = 子带内位置×Vres(未解模糊)
                tgtOut(tgtCnt).range = (tgtOut(tgtCnt).rangIdx - 1 - rngSubOff) * dims.DeltaR;
                tgtOut(tgtCnt).speed = (tgtOut(tgtCnt).dopplerIdx - 1 - dopSubOff) * dims.Vres;
                tgtOut(tgtCnt).amb_speed = (sbIdx - 1 - dopSubOff) * dims.Vres;

                tgtOut(tgtCnt).peak_type = cohPkType;                    % 相参峰态
                tgtOut(tgtCnt).nonpeakType = ncPkType;                   % 非相参峰态

                tgtOut(tgtCnt).dopplerPeakStatus = dopPkFlg;             % 多普勒局部峰标志
                tgtOut(tgtCnt).aziCnt = aziSameCnt;                      % 同方位邻峰计数
                tgtOut(tgtCnt).peakCnt = cohPkCnt;                       % 相参邻域胜出计数
                tgtOut(tgtCnt).nonPeakCnt = ncPkCnt;                     % 非相参邻域胜出计数


                % 收集该目标的64通道快照(1×512), 供DOA测角
                detCfg.doaChanIn = [detCfg.doaChanIn; single(reshape(accumCfg.chanSnapshot(:, rngIdx, sbIdx), 1, []))];



    if(tgtCnt >= maxTgtNum)      % 目标数达上限, 提前结束
        break;
    end
                pkDbgMap(rngIdx, sbIdx) = (cohPkCnt == 4) + 2* (ncPkCnt == 4);   % 调试图: 1=相参全峰, 2=非相参全峰, 3=两者
            end

            % detFlagMap改写为"粗检且峰态精判均通过"的最终判决
            detCfg.detFlagMap(rngIdx, mod(dopBin-1, dopSubNum)+1) = detPass;
        end
        if(tgtCnt >= maxTgtNum)
            break;
        end
    end
    tgtOut = tgtOut.';   % 转为列struct数组
    if tgtCnt == 0
        % 0 确认时清成【真正空】结构体: 初始化 struct('ID',[],...) 是 1×1 占位
        % (numel=1), 会骗过主函数的 isempty 判据 → 下游按"1 个空目标"索引崩溃
        tgtOut = repmat(tgtOut, 0, 1);
    end

    %% ---------- 【新增】检测三子图(粗检 / 细检detPass / 确认目标tgtCnt) ----------
    detCfg.figDet = [];                        % 不画图时句柄为空(figList 自动跳过)
    if drawFigs
    screenSize = get(0, 'ScreenSize');
    figW = screenSize(3) * 0.98;  figH = screenSize(4) * 0.80;          % 与图5(RD 3D)同尺寸
    figPos = [(screenSize(3)-figW)/2, (screenSize(4)-figH)/2, figW, figH];
    detCfg.figDet = figure('Position', figPos, 'Name', 'fig27_Detect');
    % 展开坐标: 每个像素的多普勒bin = 子带 + (阶数-1)*子带bin数(全部 0 基, 光标直读=真bin号)
    [sbGrid, rngGrid] = meshgrid(0:dopSubNum-1, 0:rngBinTotal-1);
    dopGrid = sbGrid + (accumCfg.ddmOrdMap(1:rngBinTotal, :) - 1) * dopSubNum;

    subplot(1,3,1);
    cdl = detCfg.coarseDetList;
    if isempty(cdl), cdl = zeros(0,5); end        % 空粗检(全漏检/目标出带)时画空图不崩
    scatter(cdl(:,5)-1, cdl(:,1)-1, 6, 'r', 'filled'); hold on;
    scatter(dopGrid(detCfg.detFlagMapCoarse==1), rngGrid(detCfg.detFlagMapCoarse==1), 2, 'b', 'filled');
    axis xy; grid on; xlim([0 dims.dopBinNum-1]); ylim([0 rngBinTotal-1]);
    xlabel('多普勒 bin'); ylabel('距离 bin');
    title(sprintf('粗CFAR: coarseDetList %d 点 | detFlagMapCoarse %d 格', ...
        size(detCfg.coarseDetList,1), nnz(detCfg.detFlagMapCoarse)));
    legend({'coarseDetList', 'detFlagMapCoarse'}, 'Location', 'best');

    subplot(1,3,2);
    scatter(dopGrid(detCfg.detFlagMap==1), rngGrid(detCfg.detFlagMap==1), 6, 'm', 'filled');
    axis xy; grid on; xlim([0 dims.dopBinNum-1]); ylim([0 rngBinTotal-1]);
    xlabel('多普勒 bin'); ylabel('距离 bin');
    title(sprintf('细CFAR通过点 (detPass, %d 个)', nnz(detCfg.detFlagMap)));

    subplot(1,3,3);
    if ~isempty(tgtOut)
        scatter([tgtOut.dopplerIdx]-1, [tgtOut.rangIdx]-1, 18, 'g', 'filled');
    end
    axis xy; grid on; xlim([0 dims.dopBinNum-1]); ylim([0 rngBinTotal-1]);
    xlabel('多普勒 bin'); ylabel('距离 bin');
    title(sprintf('确认目标 tgtCnt = %d', tgtCnt));
    sgtitle('两级CFAR检测(左:粗检 中:精判detPass 右:确认目标)');
    end                                             % if drawFigs
end
