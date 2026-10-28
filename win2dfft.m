function fftCfg = win2dfft(adcRawData, winType)
%WIN2DFFT 加窗 2D-FFT 模块（含 FFT 域绘图：图8 RD 图 / 图9 3D 曲面）
%   fftCfg = win2dfft(adcRawData, winType)
%   输入:
%     adcRawData [nSamples × nChirps × nRX] 时域 ADC 数据(int16 counts 口径)
%     winType    'Hanning' | 'SVA' —— 二选一, 不能同时加:
%                  'Hanning': 时域加汉宁窗(距离维+多普勒维)后再做 2D-FFT
%                  'SVA'    : 矩形窗 2D-FFT 后, 频域逐 bin 三点自适应旁瓣抵消
%                             (sva_fun_pipe.m, 移植自 参考\, 两维都启用)
%   输出(结构体, 字段约定与 参考\fft2D_fun_pip 兼容):
%     fftCfg.winType   本次使用的窗类型
%     fftCfg.data_2D   矩形窗 2D-FFT 结果 [nV × nRX × nChirps] 复数(参考基准)
%     fftCfg.data_sva  加窗处理结果(复数)——【字段名沿用 data_sva, Hanning 时为
%                      加窗 FFT 结果, SVA 时为 SVA 结果; 后续噪底/CFAR/DOA 按此取数】
%     fftCfg.rdDB      20log10|data_sva| dB 谱
%     fftCfg.figRD     图4(RD 图)句柄
%     fftCfg.figSurf   图5(3D 曲面)句柄
    nS = size(adcRawData, 1);
    nV = floor(nS/2) + 1;                    % 实采样有效距离 bin 数(1024 → 513)
    nRX = size(adcRawData, 3);
    nChirps = size(adcRawData, 2);
    X  = double(permute(adcRawData, [1 3 2]));  % → [nSamples × nRX × nChirps](int16 输入时转 double, fft 不吃整型)

    % ---- 矩形窗 2D-FFT(参考基准, 始终输出) ----
    data_2D = fft(fft(X, [], 1), [], 3);
    data_2D = data_2D(1:nV, :, :);

    % ---- 按窗类型处理 ----
    switch upper(winType)
        case 'HANNING'                        % ---- 时域加窗 → 2D-FFT ----
            Xw = X .* hanning(nS);                        % 距离维(快时间)加窗
            Xw = Xw .* reshape(hanning(nChirps), 1, 1, []); % 多普勒维(慢时间)加窗
            data_sva = fft(fft(Xw, [], 1), [], 3);
            data_sva = data_sva(1:nV, :, :);
        case 'SVA'                            % ---- 2D-FFT → 频域 SVA ----
            data_sva = sva_fun_pipe(data_2D, 1, 1);        % 距离/多普勒两维 SVA
        otherwise
            error('win2dfft: winType 必须为 ''Hanning'' 或 ''SVA''');
    end

    % ---- dB 谱(全距离 bin 输出, 不做近直流丢弃 —— 实测链路不丢任何 bin) ----
    rdDB = 20*log10(abs(data_sva));           % dB 谱 [nV × nRX × nChirps]

    % ==================== FFT 域绘图(图14: 8通道RD / 图15: RD 3D曲面) ====================
    screenSize = get(0, 'ScreenSize');
    figW = screenSize(3) * 0.98;  figH = screenSize(4) * 0.80;
    figPos = [(screenSize(3)-figW)/2, (screenSize(4)-figH)/2, figW, figH];

    % ----- 图14: 距离-多普勒图（8 个 RX 通道） -----
    fig2 = figure('Position', figPos, 'Name', 'fig14_RD_8RX');
    tiledlayout(2, ceil(nRX/2), 'TileSpacing', 'compact', 'Padding', 'compact');
    floorEst = median(rdDB(isfinite(rdDB(:))));      % 底噪估计(数据中位)
    cRange = [max(5*floor((floorEst-10)/5), 0), 5*ceil((max(rdDB(:))+5)/5)];  % 底噪−10 ~ 峰值+5
    for i = 1:nRX
        nexttile;
        imagesc(0:nChirps-1, 0:nV-1, squeeze(rdDB(:, i, :)));
        axis xy;
        title(sprintf('RX %d', i));
        xlabel('多普勒 bin');
        ylabel('距离 bin');
        colorbar;
        clim(cRange);
    end
    sgtitle(sprintf('Range-Doppler Map（%s 窗, 实采样 + 模拟LPF）', upper(winType)));

    % ----- 图15: RD 幅度 3D 曲面（与图14同一数据, 可旋转, 默认视角） -----
    fig3 = figure('Position', figPos, 'Name', 'fig15_RD_surf');
    tiledlayout(2, ceil(nRX/2), 'TileSpacing', 'compact', 'Padding', 'compact');
    [dopGrid, rngGrid] = meshgrid(0:nChirps-1, 0:nV-1);
    for i = 1:nRX
        nexttile;
        surf(rngGrid, dopGrid, squeeze(rdDB(:, i, :)), 'EdgeColor', 'none');
        switch i
            case 1, view(90, 0);
            case 2, view(0, 90);
            case 3, view(45, 45);
            case 5, view(30, 60);
            case 6, view(-40, 5);
            otherwise, view(0, 0);
        end
        zlim(cRange);                     % 与图14 相同的动态范围
        title(sprintf('RX %d', i));
        xlabel('距离 bin');
        ylabel('多普勒 bin');
        zlabel('幅度 (dB)');
        grid on;
    end
    sgtitle(sprintf('RD 幅度 3D 曲面（%s 窗, 距离 bin-多普勒 bin-幅度dB）', upper(winType)));

    % ---- 输出结构体 ----
    fftCfg = struct('winType', upper(winType), 'data_2D', data_2D, ...
                    'data_sva', data_sva, 'rdDB', rdDB, ...
                    'figRD', fig2, 'figSurf', fig3);
end
