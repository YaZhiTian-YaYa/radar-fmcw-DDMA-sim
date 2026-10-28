function figFilt = filtChainDraw(fc, fs, fs2)
%FILTCHAINDRAW 滤波链幅频曲线图(排查用; 主流程图序 fig12, 显式编号钉死)
%   figFilt = filtChainDraw(fc, fs, fs2)
%   fc   = 滤波链设计参数(rf_sim_core 的 rfInfo 子集, 或校准缓存里的 filtChain):
%          .lpfFc .lpfOrder(模拟) .nDStg .hbH{cell} .hbR .hbFc(数字半带)
%          .chipOsr .fsAdc .B_noise
%   fs   = 输出采样率 (Hz); fs2 = 仿真过采样率 (Hz) —— 横轴上限 = fs2/2
%   图窗显式 figure(12): 编号不依赖创建顺序(主流程在校准遍时域图之前调用,
%   时序与编号天然一致; 重跑自动复用同一窗口)。
    screenSize = get(0, 'ScreenSize');                       % 与其它图统一尺寸(0.98×0.80 屏, 居中)
    figW = screenSize(3) * 0.98;  figH = screenSize(4) * 0.80;
    figPos = [(screenSize(3)-figW)/2, (screenSize(4)-figH)/2, figW, figH];
    figFilt = figure(12);                                    % ★ 显式编号(旧版语法: 编号单独传)
    set(figFilt, 'Name', 'fig12_FilterChain', 'Position', figPos);
    clf(figFilt); hold on;

    hbH = fc.hbH(1:max(fc.nDStg,1));  hbR = fc.hbR(1:max(fc.nDStg,1));
    hbFc = fc.hbFc(1:max(fc.nDStg,1));
    nnOf = @(h) -(numel(h)-1)/2 : (numel(h)-1)/2;            % FIR 对称中心轴(行向量)
    hbMagOf = @(sI, fHz) arrayfun(@(fq) abs(sum(hbH{sI} .* ...
        exp(-1j*2*pi*(fq/hbR(sI))*nnOf(hbH{sI})))), fHz);

    % 模拟 LPF(蓝, 画到仿真奈奎斯特)
    fA1 = linspace(1e5, fs2/2, 3000);
    gA = @(fb) 1./sqrt(1 + (fb/fc.lpfFc).^(2*fc.lpfOrder));
    plot(fA1/1e6, 20*log10(gA(fA1) + eps), 'b', 'LineWidth', 1.2);
    legTxt = {sprintf('模拟LPF %d 阶 @%.0fM (0.45×ADC率)', fc.lpfOrder, fc.lpfFc/1e6)};

    % 各级数字半带(红/品红, 各自域)
    hbClr = {'r', 'm'};
    for s = 1:fc.nDStg
        fAs = linspace(1e5, hbR(s)/2, 2500);
        plot(fAs/1e6, 20*log10(hbMagOf(s, fAs) + eps), hbClr{s}, 'LineWidth', 1.2);
        legTxt{end+1} = sprintf('数字半带%d: %d 抽头 @%gM 域(截止 %gM)', ...
                                s, numel(hbH{s}), hbR(s)/1e6, hbFc(s)/1e6); %#ok<AGROW>
    end

    % 级联总响应(黑粗, 带内口径到 fs) —— 数据实际看到的形状
    fAc = linspace(1e5, fs, 1500);
    casc = gA(fAc);
    for s = 1:fc.nDStg, casc = casc .* hbMagOf(s, fAc); end
    plot(fAc/1e6, 20*log10(casc + eps), 'k', 'LineWidth', 1.8);
    legTxt{end+1} = '级联总响应(数据实际看到)'; %#ok<AGROW>

    % 带边噪声均衡 + 噪底实际形状(照 96 修法; 均衡曲线为绿, 噪底形状为青虚线)
    % 噪底形状 = 白噪经链路 + 两级抽取【折叠到输出 bin 的合计幅度】× 均衡 —— 仿真口径,
    % 与芯片实测噪底-距离剖面对照, 均衡生效时应全带平(未均衡时末级过渡带下垂 ~2dB)
    if isfield(fc, 'bComp') && numel(fc.bComp) > 1
        gCe = abs(freqz(fc.bComp, 1, fAc, fs));
        plot(fAc(fAc <= fs/2)/1e6, 20*log10(gCe(fAc <= fs/2) + eps), 'g', 'LineWidth', 1.2);
        legTxt{end+1} = sprintf('带边噪声均衡 %d 抽头 @输出域(冻结 v=%.2f)', ...
                                numel(fc.bComp), fc.compFreezeV); %#ok<AGROW>
        plot(fAc/1e6, 20*log10(casc .* gCe + eps), 'k--', 'LineWidth', 1.2);
        legTxt{end+1} = '均衡后信号响应(级联×均衡)'; %#ok<AGROW>
        gCh = @(fb) gA(fb) .* arrayfun(@(fq) prod(arrayfun(@(sI) hbMagOf(sI, fq), 1:fc.nDStg)), fb);
        vG = linspace(0, 1, 513).';  fG = vG*fs/2;
        Hp = zeros(size(fG));
        for kA = 0:fc.chipOsr-1
            Hp = Hp + gCh(kA*fs + fG).^2 + gCh((kA+1)*fs - fG).^2;
        end
        noiseResp = sqrt(Hp) .* abs(freqz(fc.bComp, 1, fG, fs));
        plot(fG(2:end)/1e6, 20*log10(noiseResp(2:end) + eps), 'c-', 'LineWidth', 1.8);
        legTxt{end+1} = '噪底实际形状(折叠对合计×均衡, 应平——对齐芯片实测)'; %#ok<AGROW>
    end

    xline(fs/2/1e6, 'r--', '信号带边 fs/2', 'LabelOrientation', 'horizontal');
    xline(fs/1e6, 'k:', 'fs', 'LabelOrientation', 'horizontal');
    if fc.chipOsr > 1
        xline(fc.fsAdc/2/1e6, 'k--', 'ADC 奈奎斯特', 'LabelOrientation', 'horizontal');
    end
    grid on; xlim([0, fs2/2/1e6]); ylim([-100, 2]);
    xlabel('频率 (MHz)'); ylabel('|H| (dB)');
    legend(legTxt, 'Location', 'southeast');   % 右下角: 不挡曲线主体
    titleTxt = sprintf(['两段式滤波链 chipOsr=%d: 模拟@%.0fM + %d 级半带, ' ...
                        '输出 NEB=%.1fM, 末级悬崖 @%gM'], ...
                       fc.chipOsr, fc.lpfFc/1e6, fc.nDStg, fc.B_noise/1e6, hbFc(end)/1e6);
    if isfield(fc, 'hbCompOn') && fc.hbCompOn
        titleTxt = [titleTxt, sprintf(', 带边噪声均衡 %d 抽头', numel(fc.bComp))];
    end
    title(titleTxt);
end
