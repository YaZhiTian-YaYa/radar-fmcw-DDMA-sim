function figW = fmcwWaveDraw(prm, dims, tWait, tHold)
%FMCWWAVEDRAW FMCW+步进频波形图(图11: 频率 vs 时间, 2 行 1 列)
%   figW = fmcwWaveDraw(prm, dims, tWait, tHold)
%   子图1: 全帧 512 个 chirp 的 f-t 折线(步进频整体趋势一目了然);
%   子图2: 前 10 个 chirp 放大(单 chirp 结构: 等待段[波形已扫]/ADC 采样窗/保持段/保护)。
%   波形约定(与主函数/射频模块一致):
%     每 chirp 物理发射 = [等待 tWait(扫频中, ADC 未开)] + [ADC 窗 Ts(扫 B)] + [保持 tHold(继续扫)]
%     + [保护剩余 Tc−Ts−tWait−tHold(停发)] —— 等待/保持也是发射时间, 计入占用频带
%     起始频率按步进方向步进: f_m = fStart + stepDir·(m−1)·f_step
    Ts = prm.Ts;  Tc = prm.Tc;  B = prm.B;
    slope = B/Ts;
    fstep = prm.stepDir * prm.cr*B/(dims.dopBinNum-1);
    if prm.stepDir > 0, dirTxt = '上跳频'; else, dirTxt = '下跳频'; end
    nCh   = dims.dopBinNum;
    tGuard = Tc - Ts - tWait - tHold;              % 停发保护时隙

    screenSize = get(0, 'ScreenSize');
    figPos = [screenSize(3)*0.15, screenSize(4)*0.06, screenSize(3)*0.70, screenSize(4)*0.84];
    figW = figure('Position', figPos, 'Name', 'fig11_FMCW_wave');
    tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

    % ---- 子图1: 全帧 ----
    nexttile;  hold on;
    for m = 1:nCh
        t0 = (m-1)*Tc;
        f0 = prm.fStart + (m-1)*fstep;
        % 发射段 = 等待起点 → 保持终点(连续扫频); 保护段停发不画线
        tt = t0 + [0, tWait+Ts+tHold];
        ff = f0 + slope*[0, tWait+Ts+tHold];
        plot(tt*1e6, ff/1e9, 'b-', 'LineWidth', 0.3);
    end
    grid on;
    xlabel('时间 (\mus)');  ylabel('射频频率 (GHz)');
    title(sprintf('全帧 %d 个 chirp(步进频 %s, 每 chirp 起始步进 %.3f MHz)', nCh, dirTxt, abs(fstep)/1e6));
    % ---- 子图2: 前 nShow 个 chirp 放大(画单 chirp 四段结构 + chirp 间跳变线) ----
    nexttile;  hold on;
    nShow = min(3, nCh);      % [旋钮] 放大显示的 chirp 数: 当前斜率下 10 个会把锯齿压成
                              %   近竖线(单 chirp 扫频只占周期的 6%), 3 个才能看清 FMCW 形状
    for m = 1:nShow
        t0 = (m-1)*Tc;
        f0 = prm.fStart + (m-1)*fstep;
        plot((t0+[0, tWait])*1e6, (f0+slope*[0, tWait])/1e9, 'Color', [0.85 0.5 0.1], 'LineWidth', 1.2);   % 等待段
        plot((t0+tWait+[0, Ts])*1e6, (f0+slope*(tWait+[0, Ts]))/1e9, 'b-', 'LineWidth', 1.8);              % ADC 窗(扫 B)
        plot((t0+tWait+Ts+[0, tHold])*1e6, (f0+slope*(tWait+Ts+[0, tHold]))/1e9, 'Color', [0.5 0.7 0.2], 'LineWidth', 1.2); % 保持段
        if tGuard > 0                                                        % 保护段: 停发(虚线占位, 同除 1e9 用 GHz)
            plot((t0+tWait+Ts+tHold)*1e6 + [0, tGuard]*1e6, ...
                 (f0 + slope*(tWait+Ts+tHold) + [0, 0])/1e9, 'k:', 'LineWidth', 0.8);
        end
        if m < nShow                                                          % chirp 间频率跳变(瞬时回扫, 非扫频)
            tEnd = t0 + tWait + Ts + tHold;
            fEnd = f0 + slope*(tWait+Ts+tHold);
            f1n = prm.fStart + m*fstep;
            plot([tEnd, tEnd]*1e6, [fEnd, f1n]/1e9, 'Color', [0.6 0.6 0.6], 'LineStyle', '--', 'LineWidth', 0.8);
        end
    end
    grid on;
    xlabel('时间 (\mus)');  ylabel('射频频率 (GHz)');
    legend({'等待段(ADC 未开, 已在扫频)', 'ADC 采样窗(扫 B)', '保持段(ADC 关, 继续扫)', ...
            '保护时隙(停发)', 'chirp 间跳变(瞬时)'}, 'Location', 'northwest', 'FontSize', 7);
    title(sprintf('前 %d 个 chirp 放大: 等待 %.1f / 采样 %.2f / 保持 %.1f / 保护 %.2f \\mus, 周期 Tc=%.0f \\mus, 斜率 %.2f MHz/\\mus', ...
        nShow, tWait*1e6, Ts*1e6, tHold*1e6, tGuard*1e6, Tc*1e6, slope/1e12));
    sgtitle(sprintf('FMCW+步进频波形(频率-时间): 每 chirp = 等待+采样+保持+保护, 步进 %.3f MHz/chirp(%s)', ...
        abs(fstep)/1e6, dirTxt));
end
