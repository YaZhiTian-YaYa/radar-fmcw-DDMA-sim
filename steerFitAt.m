function atom = steerFitAt(A, angCoarse, angTgt)
%STEERFITAT 对指定角度做导向矢量拟合(表格插值, 按需调用)
%   atom = steerFitAt(A, angCoarse, angTgt)
%   从 1° 间隔的导向矢量表 A(通道×列) 拟合出任意角度的导向矢量(通道×1)。
%   【按需调用】: 每次只返回单个通道×1 向量, 绝不预生成更细的导向矢量表 ——
%   雷达内存没有额外空间放大表(工程约束)。
%   拟合规律(符合物理): 每阵元相位随角度线性变化(1° 内), 步进 = 表中相邻
%   两列的逐元素相位差(直接读表, 不需要阵元坐标/波长先验); 幅度取格点值。
%   输入: A        导向矢量表(通道×列; 共轭域由调用方约定, 拟合规律不变)
%         angCoarse 表的列角度 (deg, 1° 间隔, 单调)
%         angTgt    目标角度 (deg, 任意; 超出表范围按最近段外推)
%   输出: atom      拟合的导向矢量 [通道×1]
    M = numel(angCoarse);
    if isnan(angTgt)                             % NaN 角度 → 返回零向量
        atom = zeros(size(A, 1), 1);
        return;
    end
    % 最近列 + 亚列偏移(列单位, ±0.5 内) —— 算术直查(等价 interp1 nearest+限幅,
    % 但无函数调用开销: 本函数被 OMP/DML 逐角度调用千万次, interp1 是热点)
    p0 = round((angTgt - angCoarse(1)) / (angCoarse(2) - angCoarse(1))) + 1;
    p0 = min(max(p0, 1), M);                      % 限幅在表范围内(范围外就近外推)
    d  = angTgt - angCoarse(p0);                  % 亚列偏移 ∈ [-0.5, 0.5]
    if d >= 0
        pc = min(p0, M-1);                       % 段 [pc → pc+1]
        beta = angle(A(:, pc+1) .* conj(A(:, pc)));
        atom = A(:, p0) .* exp(1i * beta * d);
    else
        pc = max(p0-1, 1);                       % 段 [pc-1 → pc]
        beta = angle(A(:, p0) .* conj(A(:, pc)));
        atom = A(:, p0) .* exp(-1i * beta * (-d));
    end
end
