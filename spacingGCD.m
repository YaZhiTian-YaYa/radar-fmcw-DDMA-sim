function dz = spacingGCD(coordsM)
%SPACINGGCD 阵列坐标的【基本间距】= 相邻阵元间距的最大公约数
%   dz = spacingGCD(coordsM)
%   例: 阵元位于 0, 3, 5 (×基线) → 间距 3, 2 → 最大公约数 1 (×基线)。
%   内部按 1 nm 量化成整数后做 gcd, 避免浮点求和尾数破坏等值判断。
    u = unique(round(coordsM(:)/1e-9));      % 1 nm 量化(整数, 去重升序)
    d = diff(u);                             % 相邻间距(1nm 整数)
    g = d(1);
    for kk = 2:numel(d)
        g = gcd(g, d(kk));
    end
    dz = g*1e-9;                             % 还原为米
end
