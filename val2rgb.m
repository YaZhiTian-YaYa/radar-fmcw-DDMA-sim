function C = val2rgb(v)
%VAL2RGB 数值向量 → N×3 RGB 颜色(经 parula 色图映射)
%   规避 scatter/scatter3 的 N=3 陷阱: 数值颜色数据恰为 3 个元素时被 MATLAB 解释为
%   单点 RGB 三元组, 值域在 [0,1] 外(如高度/距离)即报"RGB 三元组无效"且全部点不着
%   色(2026-10-29 三目标场景踩中); 先转 N×3 RGB 矩阵则任意 N 都按逐点颜色处理。
%   配套建议: 调用侧 clim([min(v) max(v)]) 保持 colorbar 量程一致。
    v = v(:).';
    if isempty(v), C = zeros(0, 3); return; end
    vm = (v - min(v)) / (max(v) - min(v) + eps);
    cm = parula(256);
    C  = cm(max(1, min(256, round(vm*255) + 1)), :);
end
