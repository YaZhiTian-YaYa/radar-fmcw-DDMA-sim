function [accumCfg] = jointAccum_pipe(fftCfg)
%JOINTACCUM_PIPE 非相参积累 + DDM解模糊选阶 + 通道重排 + 方位维FFT
%   accumCfg.nonCohMap / ddmOrdMap 由非相参积累给出, 其余为相参积累输出
%   输出: accumCfg.aziFftSpec -- 方位维FFT谱 [128方位 × 512距离 × 32子带]
%         accumCfg.cohPeakPow / cohPeakPos -- 各单元方位谱峰值功率/位置
%         accumCfg.chanSnapshot -- 各单元64通道复数快照(供测角)

    %% 第一步: 非相参积累(内含DDM逐阶对齐与选阶, 得到nonCohMap/ddmOrdMap)
    accumCfg = nonCohAccumulate_fun_pipe(fftCfg);
    %% 第二步: 相参积累(通道校正 -> 按虚拟阵元位置重排 -> 方位维FFT)
    % chanPosLut: 64个物理接收通道 -> 128点方位FFT中虚拟阵元位置查找表(稀疏阵)
    chanPosLut = [59,71,56,48,42,43,39,37,54,65,51,43,37,38,34,32,52,63,49,41,35,36,32,30,47,58,44,36,30,31,27,25,44,55,41,33,27,28,24,22,34,45,31,23,17,18,14,12,31,42,28,20,14,15,11,9,24,35,21,13,6,8,4,1];
    chanTotal = 64;                           % 接收通道数
    accumCfg.aziFftLen = 128;                 % 方位维FFT长度
    aziFftIn = zeros(128,512,32);             % 方位FFT输入缓冲 [128方位 × 512距离 × 32子带]
    accumCfg.chanSnapshot = zeros(64,512,32); % 各单元64通道快照缓存(供测角)
    sigPowBuf = zeros(512,32,2);              % (预留缓冲, 当前未使用)
    noisePowBuf = zeros(512,32,2);            % (预留缓冲, 当前未使用)
    chanCalib = [-0.6563225 - 0.7544804i;0.7639087 - 0.6453244i;0.6868339 - 0.7268144i;0.9716317 + 0.2364990i;-0.2534203 - 0.9673563i;-0.4040243 + 0.9147483i;-0.6518335 + 0.7583621i;1.000000 + 0.000000i;-0.3520049 - 0.9359981i;0.9355057 - 0.3533116i;0.8959444 - 0.4441662i;0.8341747 + 0.5515003i;0.1029662 - 0.9946849i;-0.6965392 + 0.7175188i;-0.8747044 + 0.4846569i;0.9401693 + 0.3407078i;-0.6077918 + 0.7940964i;-0.7888696 - 0.6145607i;-0.8512695 - 0.5247288i;0.02896841 - 0.9995803i;-0.8981044 + 0.4397823i;0.9775999 + 0.2104717i;0.8742810 + 0.4854201i;-0.2039834 - 0.9789743i;-0.8487872 - 0.5287346i;0.5064853 - 0.8622486i;0.4426017 - 0.8967183i;0.9980158 - 0.06296471i;-0.5303631 - 0.8477706i;-0.1182307 + 0.9929861i;-0.3977887 + 0.9174771i;0.9502986 - 0.3113399i;0.1755386 - 0.9844726i;0.9855856 + 0.1691775i;0.9973145 + 0.07323819i;0.4408793 + 0.8975664i;0.5974851 - 0.8018800i;-0.9649707 + 0.2623577i;-0.9996150 - 0.02774708i;0.6375585 + 0.7704020i;-0.1952615 - 0.9807512i;0.9767846 - 0.2142239i;0.9558443 - 0.2938735i;0.7353879 + 0.6776464i;0.2601675 - 0.9655635i;-0.8036491 + 0.5951034i;-0.9407624 + 0.3390666i;0.8772300 + 0.4800704i;0.5872203 - 0.8094271i;0.8209499 + 0.5710002i;0.8637475 + 0.5039248i;-0.003665125 + 0.9999933i;0.8845810 - 0.4663866i;-0.9823200 - 0.1872096i;-0.8871720 - 0.4614388i;0.2420913 + 0.9702535i;0.6589521 + 0.7521849i;-0.7532191 + 0.6577697i;-0.6846745 + 0.7288490i;-0.9717142 - 0.2361599i;0.2599990 + 0.9656089i;0.4013085 - 0.9159430i;0.6497135 - 0.7601792i;-0.9999731 + 0.007330264i];                              % (64×1通道幅相校正系数, 乘在通道数据上)
    ddmShift = [0    14    12     9     8     5     4     1];    % 8个TX的DDM多普勒频移(单位:组, 与非相参积累中ddmList一致)
    dopBinNum = 512; grpTotal = 16;          % 多普勒总bin数 / DDM分组数
    %% 通道重排与方位FFT: 逐(距离, 多普勒子带)单元取回通道数据并做方位FFT
    for rngIdx = 1 : 512
        for dopIdx = 1 : 32
            % 该单元选中的DDM阶数
            txOrd = accumCfg.ddmOrdMap(rngIdx, dopIdx);
            % 由阶数txOrd与各TX频移ddmShift反推目标信号所在的8个chirp位置:
            % chirp序号 = 子带序号 + (阶数+TX频移)×32(每组32chirp), 循环回绕到1~512
            chirpPos = mod(dopIdx + (txOrd + ddmShift - 1) * size(accumCfg.ddmOrdMap, 2) - 1, ...
                             512) + 1;
            % 保存该单元64通道×8chirp的复数快照(展成1×512, 供后续测角)
            accumCfg.chanSnapshot(:, rngIdx, dopIdx) = reshape(squeeze(fftCfg.data_sva(rngIdx, :, chirpPos)), 1, []);
            % 同一快照做通道幅相校正
            aziCompBuf = reshape(squeeze(fftCfg.data_sva(rngIdx, :, chirpPos)), [], 1) .* chanCalib;


            % 按虚拟阵元位置查找表, 把各通道数据累加进128点方位FFT输入的对应位置
            % (同一位置可能被多个通道共享, 用累加实现阵元合成)
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

end
