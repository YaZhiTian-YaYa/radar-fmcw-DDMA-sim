function Minv = inv2x2(M)
%INV2X2 2×2 复矩阵求逆 —— 伴随矩阵闭式法(嵌入式资源友好)
%   Minv = inv2x2(M)
%   M = [a b; c d] → M⁻¹ = (1/det)·[d −b; −c a]
%   资源: 4 次复乘 + 2 次复加 + 1 次复除, 无通用矩阵求逆/无 SVD/无循环
%   —— 嵌入式实现的标准做法(OMP/DML 等小规模最小二乘通用)。
%   注: 这是闭式精确解(非近似); 相比 pinv 的 SVD 分解, 计算量与内存小一个量级。
    a = M(1,1);  b = M(1,2);
    c = M(2,1);  d = M(2,2);
    detM = a*d - b*c;
    Minv = [d, -b; -c, a] / detM;
end
