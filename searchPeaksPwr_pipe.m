function [peakPwr, peakIdx, alpha] = searchPeaksPwr_pipe(pwrList, peakNum)
%SEARCHPEAKSPWR_PIPE 在一维功率谱(dB)上搜索前peakNum个峰, 并做亚bin插值
%   输入: pwrList -- 一维谱(通常为dB)
%   输出: peakPwr -- 前peakNum个峰的幅值(降序)
%         peakIdx -- 对应峰的索引
%         alpha   -- 各峰的亚bin偏移(抛物线插值)

    % 找出所有严格局部峰
    [~, peakList] = findpeaks(pwrList);
    % 左端点若高于右邻也作为候选峰(端点处findpeaks不检测)
    if pwrList(1) > pwrList(2)
        peakList = [1; peakList];
    end
    if pwrList(end) >= pwrList(end - 1), peakList = [peakList; length(pwrList)]; end
    % 按峰幅值降序取前peakNum个
    [sortVal, sortIdx] = sort(pwrList(peakList), 'descend');
    peakPwr = sortVal(1 : peakNum);
    peakIdx = peakList(sortIdx(1 : peakNum));
    % 逐峰亚bin偏移: 抛物线(二阶差分)插值
    alpha = peakIdx;
    for ii = 1 : length(peakIdx)
        if peakIdx(ii) == 1 || peakIdx(ii) == length(pwrList)
            alpha(ii) = 0;      % 端点峰无法插值, 偏移置0
        else
            alpha(ii) = (pwrList(peakIdx(ii) - 1) - pwrList(peakIdx(ii) + 1)) / ...
                        (2 * (pwrList(peakIdx(ii) - 1) + pwrList(peakIdx(ii) + 1)) - 4 * pwrList(peakIdx(ii)));
        end
    end
end
