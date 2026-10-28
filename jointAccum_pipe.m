function [accumCfg] = jointAccum_pipe(fftCfg, ddmShift, cCal, A_manifold, FFT_d0, dims, drawFigs)
%JOINTACCUM_PIPE 非相参积累 + DDM解模糊选阶 + 通道重排 + 方位维FFT(64元流型版)
%   (参考 参考\rsp_pipeline_ReWrite\pipeline\jointAccum_pipe.m, 适配点:
%    ddmShift 用我们自己的 DDMA 编码(find(encodeSeq)-1);
%    chanCalib 用我们自己的 0° 校准表 cCal(64×1, TX内循环/RX外循环顺序);
%    chanPosLut 由我们自己的 64 单元阵列流型 A_manifold 的 x 坐标生成——
%    相参积累用全部 64 个虚拟阵元(方位维, 不考虑俯仰: 车载俯仰范围很小,
%    即便有俯仰误差, 64 元效果也优于只用单行 30 元, 更优于非相参积累)。)
%   accumCfg.nonCohMap / ddmOrdMap 由非相参积累给出, 其余为相参积累输出
%   输出: accumCfg.aziFftSpec -- 方位维FFT谱 [128方位 × rngUse距离 × subBins(=16)子带]
%         accumCfg.cohPeakPow / cohPeakPos -- 各单元方位谱峰值功率/位置
%         accumCfg.chanSnapshot -- 各单元64通道复数快照(供测角)
%         accumCfg.figAccum -- 积累图句柄(不画图时为空)
%   drawFigs: false = 只算不画(校准遍只取积累结果, 主流程两遍测量用)
    if nargin < 7 || isempty(drawFigs), drawFigs = true; end   % 默认画图

    %% 第一步: 非相参积累(内含DDM逐阶对齐与选阶, 得到nonCohMap/ddmOrdMap)
    accumCfg = nonCohAccumulate_fun_pipe(fftCfg, ddmShift, dims);
    chanTotal = dims.chanTotal;               % 接收通道数(TX×RX, 主函数传入)
    accumCfg.aziFftLen = dims.aziFftLen;      % 方位维FFT长度(主函数传入)

    %% 第二步: 相参积累(通道校正 -> 按虚拟阵元位置重排 -> 方位维FFT)
    % chanPosLut: 64个虚拟通道 -> 128点方位FFT中阵元位置查找表(由 A_manifold 的 x 坐标生成;
    % 同一 x 位置的多个通道用累加实现阵元合成; 平移到从 1 开始)
    xvGrid = A_manifold(1, :) / FFT_d0;                      % 虚拟阵元 x(格)
    chanPosLut = double(round(xvGrid - min(xvGrid)) + 1);    % [1×64], TX内循环/RX外循环
    chanCalib = cCal(:);                                     % 0°校准表(64×1, 顺序与快照一致)
    aziFftIn = zeros(dims.aziFftLen, dims.rngUse, dims.subBins);      % 方位FFT输入 [方位×距离×子带](第二维是距离!)
    accumCfg.chanSnapshot = zeros(dims.chanTotal, dims.rngUse, dims.subBins);   % 64通道快照(供测角, 第二维是距离)
    %% 通道重排与方位FFT: 逐(距离, 多普勒子带)单元取回通道数据并做方位FFT
    for rngIdx = 1 : dims.rngUse               % 距离维遍历(nV-1 个有效 bin; 勿写 dopBinNum —— 原 512==512 巧合)
        for dopIdx = 1 : dims.subBins
            % 该单元选中的DDM阶数
            txOrd = accumCfg.ddmOrdMap(rngIdx, dopIdx);
            % 由阶数txOrd与各TX频移ddmShift反推目标信号所在的8个多普勒bin位置:
            % bin序号 = 子带序号 + (阶数+TX频移)×子带bin数(每组dims.subBins个bin), 回绕到速度全带
            chirpPos = mod(dopIdx + (txOrd + ddmShift - 1) * dims.subBins - 1, ...
                             dims.dopBinNum) + 1;
            % 取8TX×8RX=64通道快照, 排成 TX内循环/RX外循环(与 chanCalib/A_manifold 一致)
            blk = squeeze(fftCfg.data_sva(rngIdx, :, chirpPos));   % [nRX × nTX]
            accumCfg.chanSnapshot(:, rngIdx, dopIdx) = reshape(blk.', [], 1);
            % 同一快照做通道幅相校正(0°校准表)
            aziCompBuf = reshape(blk.', [], 1) .* chanCalib;

            % 按虚拟阵元位置查找表, 把各通道数据累加进128点方位FFT输入的对应位置
            for cIdx = 1 : chanTotal
                aziFftIn(chanPosLut(cIdx), rngIdx, dopIdx) = aziFftIn(chanPosLut(cIdx), rngIdx, dopIdx) + aziCompBuf(cIdx);
            end
        end
    end
    accumCfg.aziFftBuf = aziFftIn;                            % 方位FFT输入(调试用)
    accumCfg.aziFftSpec = single(fft(aziFftIn, [], 1));        % 沿方位维128点FFT(即DBF)
    % 各单元取方位谱峰值功率与对应方位bin位置
    [accumCfg.cohPeakPow, accumCfg.cohPeakPos] = max(abs(accumCfg.aziFftSpec) .^ 2, [], 1);
    accumCfg.cohPeakPow = single(squeeze(accumCfg.cohPeakPow));   % 相参积累峰值功率 [512×32]
    accumCfg.cohPeakPos = single(squeeze(accumCfg.cohPeakPos));   % 峰值方位bin位置   [512×32]
    accumCfg.aziFftPowDb = single(db(abs(accumCfg.aziFftSpec))); % 方位谱dB形式(供CFAR检测)
    accumCfg.cohGainDb = 10*log10(accumCfg.cohPeakPow ./ accumCfg.nonCohMap);  % 相参相对非相参的积累增益(dB)

    %% 第三步: 绘图(一张大图 2×4: 上排非相参 / 下排相参, 各 4 个视角的 3D 曲面)
    % 视角依次参考 RD 3D 图(RX4/RX3/RX1/RX2 的视角), 免旋转直接观察;
    % 速度维坐标显示整个速度带(0~dopBinNum−1 bin)——数据只占 subBins(=16) 个子带位置,
    % 一眼即可看出速度维元素不满(图是按子带画的)。
    accumCfg.figAccum = [];                     % 不画图时句柄为空(figList 自动跳过)
    if drawFigs
    screenSize = get(0, 'ScreenSize');
    figW = screenSize(3) * 0.98;  figH = screenSize(4) * 0.80;
    figPos = [(screenSize(3)-figW)/2, (screenSize(4)-figH)/2, figW, figH];
    viewList = [0 0; 45 45; 90 0; 0 90];        % 依次 = RD3D 图中 RX4/RX3/RX1/RX2 的视角
    viewName = {'视角1', '视角2', '视角3', '视角4'};
    nDb = 10*log10(accumCfg.nonCohMap);          % [rngUse距离 × subBins(=16)子带] dB
    cDb = 10*log10(accumCfg.cohPeakPow);
    [dopGrid, rngGrid] = meshgrid(0:dims.subBins-1, 0:dims.rngUse-1);  % 子带位置(0基), xlim 显示全带; 距离 0 基(与 RD 图/gate_R 一致)
    floorN = median(nDb(:));
    zlN = [max(5*floor((floorN-10)/5), 0), 5*ceil((max(nDb(:))+5)/5)];  % 底噪~峰值自适应
    floorC = median(cDb(:));
    zlC = [max(5*floor((floorC-10)/5), 0), 5*ceil((max(cDb(:))+5)/5)];  % 底噪~峰值自适应

    accumCfg.figAccum = figure('Position', figPos, 'Name', 'fig26_JointAccum');
    tiledlayout(2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
    for v = 1 : 4
        nexttile;
        surf(dopGrid, rngGrid, nDb, 'EdgeColor', 'none');
        view(viewList(v, 1), viewList(v, 2));
        xlim([0 dims.dopBinNum-1]); zlim(zlN); colormap(jet);
        title(sprintf('非相参积累 %s', viewName{v}));
        xlabel('速度 bin(全带 0~511)'); ylabel('距离 bin'); zlabel('功率 (dB)');
        grid on;
    end
    for v = 1 : 4
        nexttile;
        surf(dopGrid, rngGrid, cDb, 'EdgeColor', 'none');
        view(viewList(v, 1), viewList(v, 2));
        xlim([0 dims.dopBinNum-1]); zlim(zlC); colormap(jet);
        title(sprintf('相参积累 %s', viewName{v}));
        xlabel('速度 bin(全带 0~511)'); ylabel('距离 bin'); zlabel('功率 (dB)');
        grid on;
    end
    sgtitle('积累图(上排: TXRX非相参 8TX对齐+8RX功率平均·DDM选阶 | 下排: TXRX相参 64元方位FFT·0°校准)');
    end                                             % if drawFigs
end
