% 多距离 bin 联合多普勒诊断(步进频斜拉): 目标±8行逐行找"非自身结构"杂散
% 每行合法结构: 本行码线±2 + 空子带±2(仅目标行) + ±1行内的目标线距离裙(同多普勒位)
% 杂散 = 其余高于 本行热噪+6dB 的 bin; 输出逐行杂散数/最高电平/位置
diary('C:\Users\NJ\AppData\Local\Temp\probe95\multiRow.log'); clc; diary on;
addpath('F:\Work3\CH_NEW_95_射频_DDMA');
saveFigDir = 'C:\Users\NJ\AppData\Local\Temp\probe95\cohfigs';
if ~exist(saveFigDir,'dir'), mkdir(saveFigDir); end
batchState.idx = 1;  batchState.okRun = true;
run('F:\Work3\CH_NEW_95_射频_DDMA\main_sim.m');
rxI = 7;
prow = mean(abs(fftCfg.data_sva(:, rxI, :)).^2, 3);  prow = prow(:);
prow(1:12) = 0;  [~, r0] = max(prow);
row0 = squeeze(20*log10(abs(fftCfg.data_sva(r0, rxI, :)) + eps));
[~, pk] = maxk(row0, 8);
fprintf('MR 目标行 %d, 码线(0基) %s, cr=%g\n', r0, mat2str(sort(pk-1)), prm.cr);
% 合法掩码(目标行的码线±2); 相邻行 ±1 内同多普勒位也是目标裙(合法)
mL = false(1,256);
for k = 1:8, mL(mod((pk(k)-1-2):(pk(k)-1+2),256)+1) = true; end
for rr = r0-8 : r0+8
    if rr < 1 || rr > size(fftCfg.data_sva,1), continue; end
    rowT = squeeze(20*log10(abs(fftCfg.data_sva(rr, rxI, :)) + eps));
    rowE = squeeze(20*log10(abs(fftCfg.data_sva(350, rxI, :)) + eps));
    fMed = median(rowE);
    legit = mL;                                   % 码线±2(目标裙/斜拉产物在此附近算合法)
    if rr == r0                                   % 目标行再排除空子带±2(交调允许)
        g16 = unique(mod(pk-1,16));
        for s = (0:15)*16 + g16(1)
            if ~ismember(mod(s,256), pk-1)
                legit(mod(mod(s,256)-2:mod(s,256)+2,256)+1) = true;
            end
        end
    end
    spur = rowT > fMed+6 & ~legit;
    if any(spur)
        [sv, si] = sort(rowT(spur), 'descend');
        top = min(4, sum(spur));
        pos = '';
        for q = 1:top, pos = [pos sprintf('%d(%.0fdB) ', si(q)-1, sv(q)-fMed)]; end %#ok<AGROW>
        fprintf('MR 行%+2d(%3d): 杂散%2d个 最高%+6.1fdB | %s\n', rr-r0, rr, sum(spur), max(rowT(spur))-fMed, pos);
    else
        fprintf('MR 行%+2d(%3d): 干净\n', rr-r0, rr);
    end
end
diary off;
