function [coarseTargetList, maxPwrIdxMatrix, maxPwrMatrix, detMap, cfarRatio] = ...
            dacfar_version2_pipe(daFftOut, nveList, coarseThreList, fineThreList, ...
            guardWin, refWin, cfarType, rng_valid,weakStatues)
%DACFAR_VERSION2_PIPE 方位-多普勒二维谱CFAR粗检(DACFAR)
%   输入: daFftOut      -- [方位 × 多普勒 × 距离] 线性幅度谱
%         nveList       -- 逐距离门噪底功率(线性)
%         coarseThreList/fineThreList -- 逐距离门粗门限(噪底倍数)/细门限(CFAR乘性系数)
%         guardWin/refWin -- [方位, 多普勒]两维的保护单元数/参考单元数
%         cfarType      -- 逐距离门CFAR类型(1=均值型, 其他=SOGO选小型)
%         rng_valid     -- 参与检测的距离门数
%         weakStatues   -- 噪底突变门标记(当前未参与判决)
%   输出: coarseTargetList -- 粗检目标列表, 每行[距离bin, 多普勒bin, 方位峰bin]
%         maxPwrIdxMatrix/maxPwrMatrix -- 各(距离,多普勒)单元的方位峰位置/峰功率
%         detMap  -- 检测标记图 [距离 × 多普勒]
%         cfarRatio -- 各单元峰功率/局部噪底 的CFAR比值
    %% 输入daFftOut维度为[angle][chirp][range]
    coarseTargetList = [];
    aNum = size(daFftOut, 1);       % 方位bin数
    dNum = size(daFftOut, 2);       % 多普勒bin数(子带数)
    rNum = size(daFftOut, 3);       % 距离bin数

    maxPwrIdxMatrix = zeros(rNum, dNum);
    maxPwrMatrix = zeros(rNum, dNum);
    detMap = zeros(rNum, dNum);
    cfarRatio = zeros(rNum, dNum);


    % ---------- 逐距离门处理 ----------
    for rIdx = 1:rng_valid
        % 该距离门的功率图 [方位 × 多普勒]
        daMap = squeeze(abs(daFftOut(:,:, rIdx)).^2);
        % 每个多普勒bin沿方位维取最大值及其位置(DBF谱峰作为检测统计量)
        [maxPwr, maxPwrIdx] = max(daMap, [], 1);

        coarseThre = nveList(rIdx) * coarseThreList(rIdx);   % 绝对门限 = 噪底×粗门限
        fineThre = fineThreList(rIdx);                       % 乘性CFAR系数

        maxPwrIdxMatrix(rIdx, :) = maxPwrIdx;
        maxPwrMatrix(rIdx, :) = maxPwr;
        if cfarType(rIdx) == 1        % 类型1: 强制方位参考窗为6, 后续走均值型CFAR
            cfarTypeRIx = 1;
            refWin(1) = 6;
        else
            cfarTypeRIx = 0;
        end

        % ---------- 逐多普勒bin做二维CFAR ----------
        for dIdx = 1: size(daFftOut, 2)

            aIdx = maxPwrIdx(dIdx);   % 该多普勒bin的方位峰位置
if refWin(1) < 1

            % ---- 分支一: 无方位参考窗, 仅在多普勒维两侧取参考单元 ----
            UIdx = warpIdx(dIdx - guardWin(2) - refWin(2) : dIdx - guardWin(2) - 1, dNum);            % 上侧(低多普勒)参考单元
            DIdx = warpIdx(dIdx + guardWin(2) + 1 : dIdx + guardWin(2) + refWin(2), dNum);            % 下侧(高多普勒)参考单元

            if cfarTypeRIx == 0
                % SOGO(选小型): 上下两侧各自平均后取较小者, 抗邻近距离/速度目标抬升噪底
                noise00 = (sum(daMap(aIdx, UIdx)) ) / ( refWin(2));
                noise11 = (sum(daMap(aIdx, DIdx)) ) / (refWin(2));
                noiseSOGO = min(noise00, noise11);
            else
                % 均值型(CA): 上下两侧平均再取均值
                noise00 = (mean(daMap(aIdx, UIdx)) + mean(daMap(aIdx, DIdx))) / 2;

                noiseSOGO = noise00 ;
            end

else
            % ---- 分支二: 方位/多普勒两维同时取参考单元 ----
            LIdx = warpIdx(aIdx - guardWin(1) - refWin(1) : aIdx - guardWin(1) - 1, aNum);           % 方位左侧参考单元
            RIdx = warpIdx(aIdx + guardWin(1) + 1 : aIdx + guardWin(1) + refWin(1), aNum);           % 方位右侧参考单元

            UIdx = warpIdx(dIdx - guardWin(2) - refWin(2) : dIdx - guardWin(2) - 1, dNum);            % 多普勒上侧参考单元
            DIdx = warpIdx(dIdx + guardWin(2) + 1 : dIdx + guardWin(2) + refWin(2), dNum);            % 多普勒下侧参考单元


            if cfarTypeRIx == 1 && rIdx*0.85 > 80
                % 远距离门(约>94)且均值型: 多普勒侧均值 与 方位侧较大均值 再平均
                noise00 = (mean(daMap(aIdx, UIdx)) + mean(daMap(aIdx, DIdx))) / 2;
                noise11 = max(mean(daMap(LIdx, dIdx)), mean(daMap(RIdx, dIdx)));
                noiseSOGO = (noise00 + noise11) / 2;

            else
                % 常规SOGO: (方位左+多普勒上)/(参考单元总数) 与 (方位右+多普勒下)/(...) 取小
                noise00 = (sum(daMap(aIdx, UIdx)) + sum(daMap(LIdx, dIdx))) / (refWin(1) + refWin(2));
                noise11 = (sum(daMap(aIdx, DIdx)) + sum(daMap(RIdx, dIdx))) / (refWin(1) + refWin(2));
                noiseSOGO = min(noise00, noise11);
            end
end

            cfarRatio(rIdx, dIdx) = maxPwr(dIdx) / noiseSOGO;    % CFAR比值(峰功率/局部噪底)

            % 双门限判决: 高于绝对门限 且 高于CFAR相对门限 -> 记入粗检列表
            if(maxPwr(dIdx) > coarseThre) && (maxPwr(dIdx) > fineThre * noiseSOGO)
                coarseTargetList =[coarseTargetList;[rIdx, dIdx, aIdx]];
                detMap(rIdx, dIdx) = 1;
            end
        end
    end

    function outIdx = warpIdx(inputIdx, N)
        % 循环回绕: 越界索引折回到[1, N]内(方位/多普勒维均按周期处理)
        outIdx = inputIdx;
        outIdx(outIdx <= 0) = outIdx(outIdx <= 0) + N;
        outIdx(outIdx > N) = outIdx(outIdx > N) - N;
    end
end
