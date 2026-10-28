function figT = adcTimeDraw(adcData, fs, figName)
%ADCTIMEDRAW ADC 时域图(每 RX 画前 16 个 chirp, 一个完整 DDMA 码周期)
%   figT = adcTimeDraw(adcData, fs, figName)
%   从射频模块的绘图段抽出成独立函数: 主流程对【校准遍/目标遍】两次射频
%   各画一张, 校准遍读缓存跳过射频时也能用缓存数据重画(图只认数据不认流程)。
%   输入: adcData [nSamples×nChirps×nRX] int16 口径实采样 counts
%         fs ADC 实采样率 (Hz); figName 图窗名(含图号, 如 'fig19_ADC_time')
    plotChirpNum = 16;    % [旋钮] 绘制前 N 个 chirp(16 = 一个完整 DDMA 码周期,
                          % 可看到码周期性的包络起伏)
    plotChirpNum = min(plotChirpNum, size(adcData, 2));
    tVec = (0:size(adcData, 1)-1).' / fs;
    screenSize = get(0, 'ScreenSize');
    figW = screenSize(3) * 0.98;  figH = screenSize(4) * 0.80;
    figPos = [(screenSize(3)-figW)/2, (screenSize(4)-figH)/2, figW, figH];
    figT = figure('Position', figPos, 'Name', figName);
    tiledlayout(2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
    yl1 = 1.10 * max(abs(adcData(:, 1:plotChirpNum, :)), [], 'all');  % 各 RX 统一纵轴(留 10% 边距)
    for rx = 1:size(adcData, 3)
        nexttile;
        plot(tVec*1e6, adcData(:, 1:plotChirpNum, rx));
        ylim([-yl1, yl1]);                        % 统一显示范围(一眼看出哪个 RX 差)
        title(['RX', num2str(rx), ' 实采样时域']);
        xlabel('快时间 (\mus)');
        ylabel('幅度 (counts)');
        legend(arrayfun(@(x) sprintf('Chirp%d', x), 1:plotChirpNum, ...
                'UniformOutput', false), 'Location', 'northoutside', ...
                'NumColumns', 4, 'FontSize', 7);
        grid on;
    end
end
