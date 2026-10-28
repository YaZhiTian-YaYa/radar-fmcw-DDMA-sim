function figRx = rxNonCoh_pipe(fftCfg, dims)
%RXNONCOH_PIPE RX 非相参积累图(8RX 直接合成一张 RD 图, DDMA 选阶之前的初查图)
%   figRx = rxNonCoh_pipe(fftCfg, dims)
%   背景: 实测数据排查的第一步往往不做 DDM 解模糊选阶 —— 先把 8 个 RX 直接
%   合起来看一张 RD 图: 目标在不在、泄漏/底噪形状对不对、有没有异常条纹,
%   一眼可判; 确认数据"活着"再进选阶与相参积累。本函数就是这个初查视图,
%   位于 jointAccum_pipe(TXRX 非相参 + 相参积累)之前。
%   累加方式(与 TXRX 非相参积累完全一致):
%       rxPow = mean( |data_sva|² , 第 2 维 8 个 RX )   —— 线性域功率平均,
%       显示时才 10·log10 转 dB。【不是】dB 域平均: dB 平均等价于几何平均,
%       会抬高底噪、压低峰值, 不是非相参积累的口径。范围维取前 rngUse 个
%   有效 bin(实采样共轭对称, 与噪底估计的 rxPowSum 同一口径)。
%   绘图: 2×4 曲面, 第一行 4 个视角, 第二行内容与第一行完全相同 —— 重复
%   画只为让图幅/视角与后面"TXTX 非相参+相参积累"图(2×4)一致, 便于对照;
%   视角序列也与该图完全一致([0 0; 45 45; 90 0; 0 90])。
%   输入: fftCfg.data_sva [nV×nRX×nChirps] 加窗 2D-FFT 复数结果(win2dfft 输出)
%         dims             公共维度(主函数唯一定义)
%   输出: figRx 图句柄(图序号由主流程 figList 统一管理)

    % ---- RX 非相参积累: 线性功率平均(与 TXRX 非相参同口径, 见函数头注释) ----
    rxPow = squeeze(mean(abs(fftCfg.data_sva(1:dims.rngUse, :, :)).^2, 2));  % [rngUse×dopBinNum]
    rxDb  = 10*log10(rxPow);                       % 显示时才转 dB

    % ---- 绘图(2×4, 下行重复上行; 视角与"积累图"一致) ----
    screenSize = get(0, 'ScreenSize');
    figW = screenSize(3) * 0.98;  figH = screenSize(4) * 0.80;
    figPos = [(screenSize(3)-figW)/2, (screenSize(4)-figH)/2, figW, figH];
    viewList = [0 0; 45 45; 90 0; 0 90];           % 依次 = 积累图中 视角1~4
    viewName = {'视角1', '视角2', '视角3', '视角4'};
    [dopGrid, rngGrid] = meshgrid(0:dims.dopBinNum-1, 0:dims.rngUse-1);  % 距离 0 基(与 RD 图/gate_R 一致)
    floorEst = median(rxDb(:));                      % 底噪估计(数据中位)
    zl = [max(5*floor((floorEst-10)/5), 0), 5*ceil((max(rxDb(:))+5)/5)];  % 底噪−10 ~ 峰值+5

    figRx = figure('Position', figPos, 'Name', 'fig25_RXNonCoh');
    tiledlayout(2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
    for row = 1:2                                  % 下行与上行同数据同视角(凑齐 2×4 图幅)
        for v = 1:4
            nexttile;
            surf(dopGrid, rngGrid, rxDb, 'EdgeColor', 'none');
            view(viewList(v, 1), viewList(v, 2));
            xlim([0 dims.dopBinNum-1]); zlim(zl); colormap(jet);
            title(sprintf('RX 非相参积累 %s', viewName{v}));
            xlabel('速度 bin(全带 0~511)'); ylabel('距离 bin'); zlabel('功率 (dB)');
            grid on;
        end
    end
    sgtitle('RX 非相参积累(8RX 线性功率平均→dB, 选阶前初查; 下行与上行内容相同, 图幅对齐积累图)');
end
