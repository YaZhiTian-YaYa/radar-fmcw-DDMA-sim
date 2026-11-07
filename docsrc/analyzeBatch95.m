% RCS 批跑结果分析: 新 28 场景(15列) vs Nov4 历史基线(13列) 逐目标对齐 + 分组统计
diary('C:\Users\NJ\AppData\Local\Temp\probe95\analyzeBatch.log'); clc; diary on;
S = load('C:\Users\NJ\AppData\Local\Temp\probe95\rerunD_state.mat');
NEW = S.RES;  nSc = numel(S.scenes);
fprintf('新批: %d 场景 ok=%d fail=%d, RES %dx%d\n', nSc, S.nOk, S.nFail, size(NEW,1), size(NEW,2));
H = load('C:\Users\NJ\AppData\Local\Temp\probe95\finalRCS_results.mat');
OLD = H.batchState.RES;
fprintf('历史(Nov4): %dx%d\n\n', size(OLD,1), size(OLD,2));
names = {S.scenes.name};
% ---- 逐目标对齐(同场景同目标号) ----
fprintf('===== 逐目标: 新dev vs 旧dev (dB) =====\n');
nCmp = 0; dAll = [];
for r = 1:size(NEW,1)
    cId = NEW(r,1); g = NEW(r,2); nm = names{cId};
    o = find(OLD(:,1)==cId & OLD(:,2)==g, 1);
    newD = NEW(r,5); oldD = NaN; if ~isempty(o), oldD = OLD(o,5); end
    newE = NEW(r,4);
    fmtE = {'%6.2f', '  miss '};  tagMiss = {'    ', 'MISS'};
    fprintf('%s g%d: rcs%+3.0f 估%s dev%+6.2f | 旧dev%+6.2f | Δ%+5.2f %s\n', nm, g, NEW(r,3), ...
        sprintf(fmtE{isnan(newE)+1}, newE), newD, oldD, newD-oldD, tagMiss{isnan(newE)+1});
    if ~isnan(newD) && ~isnan(oldD), nCmp = nCmp+1; dAll(end+1) = newD-oldD; end %#ok<AGROW>
end
fprintf('\n可比对目标 %d 个: dev位移 中位 %+.2f | 均值 %+.2f | std %.2f | min %+.2f max %+.2f\n', ...
    nCmp, median(dAll), mean(dAll), std(dAll), min(dAll), max(dAll));
% ---- 分组统计 ----
for grp = {'I','O','D'}
    idx = [];
    for k = 1:numel(names), if startsWith(names{k}, grp{1}), idx = [idx k]; end, end %#ok<AGROW>
    rows = NEW(ismember(NEW(:,1), idx), :);
    d = rows(:,5); d = d(~isnan(d));
    miss = sum(isnan(rows(:,4)));
    fprintf('%s系 %d场景 %d目标: 检出 %d miss %d | dev 中位 %+.2f std %.2f | %+6.1f~%+6.1f\n', ...
        grp{1}, numel(idx), size(rows,1), numel(d), miss, median(d), std(d), min(d), max(d));
    if ~isempty(d), fprintf('   |dev|>3dB: %d 个 (%.0f%%)\n', sum(abs(d)>3), 100*mean(abs(d)>3)); end
end
% ---- 配对角色检查(gp=1: 第一应≥第二) ----
fprintf('\n===== 配对场景角色检查(第一 vs 第二, 强者应为第一) =====\n');
gpRows = NEW(NEW(:,14)==1, :);
if ~isempty(gpRows)
    for cId = unique(gpRows(:,1))'
        pr = gpRows(gpRows(:,1)==cId, :);
        if size(pr,1)>=2
            w1 = pr(1,4); w2 = pr(2,4);
            okRole = ~isnan(w1) && ~isnan(w2) && w1 >= w2;
            tags = {'CHECK(弱为第一?)', 'OK(强为第一)'};
            fprintf('%s: 第一 %s 第二 %s %s\n', names{cId}, num2str(w1,'%.2f'), num2str(w2,'%.2f'), tags{okRole+1});
        end
    end
end
% ---- 逐场景审计差值(从日志抓) ----
fprintf('\n===== 链路预算审计差值(本轮日志) =====\n');
try
    txt = fileread('C:\Users\NJ\AppData\Local\Temp\probe95\rerunD.log');
    ca = regexp(txt, '目标\d+ @\s?(\d+) m[^\n]*?解析\s?(-?\d+\.?\d*) dB \| 仿真值\(点云 nonCohSNR[^)]*\)[^\n]*?(-?\d+\.?\d*) dB, 差 (-?[+-]?\d+\.?\d*) dB', 'tokens');
    for k = 1:numel(ca)
        fprintf('场景目标@%sm: 解析%s 仿真%s 差%s\n', ca{k}{1}, ca{k}{2}, ca{k}{3}, ca{k}{4});
    end
    caM = vertcat(ca{:});  ds = str2double(caM(:,4)).';
    fprintf('非相参差: 中位 %+.2f | 范围 %+.2f~%+.2f | n=%d\n', median(ds), min(ds), max(ds), numel(ds));
    cb = regexp(txt, '链路预算审计·相参[^\n]*?=\s?(-?\d+\.?\d*) dB \| 仿真值\(点云 cohSNR\)\s?(-?\d+\.?\d*) dB, 差 (-?[+-]?\d+\.?\d*) dB', 'tokens');
    cbM = vertcat(cb{:});  dsb = str2double(cbM(:,3)).';
    fprintf('相参差: 中位 %+.2f | 范围 %+.2f~%+.2f | n=%d\n', median(dsb), min(dsb), max(dsb), numel(dsb));
catch err
    fprintf('审计抓取异常: %s\n', err.message);
end
diary off;
