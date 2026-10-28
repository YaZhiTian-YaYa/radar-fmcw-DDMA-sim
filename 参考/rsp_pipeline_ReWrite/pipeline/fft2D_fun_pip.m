function [fftCfg] = fft2D_fun_pip(rawData)
%FFT2D_FUN_PIP 距离维与多普勒维二维FFT, 并做SVA副瓣抑制
%   输入: rawData [距离采样点 × 接收通道 × chirp]
%   输出: fftCfg.data_2D -- 二维FFT结果 [513距离bin × 64通道 × 512多普勒bin]
%         fftCfg.data_sva -- SVA处理后的频谱(后续处理均基于此数据)

% 先对距离维(dim1, 快时间采样)FFT, 再对多普勒维(dim3, chirp慢时间)FFT,
% 通道维(dim2)保持不变
data_2D = fft(fft(rawData,[],1),[],3);
% 距离维为实信号FFT, 只取正频率一半(前513个bin)即可
data_2D = data_2D(1:513,:,:);
fftCfg.data_2D = data_2D;
% SVA(Spatially Variant Apodization)自变迹: 在距离维和多普勒维上
% 利用相邻bin自适应加权, 抑制副瓣且基本不展宽主瓣
fftCfg.data_sva = sva_fun_pipe(data_2D,1,1);
end
