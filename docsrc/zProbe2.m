addpath('F:\Work3\CH_NEW_95_射频_DDMA');
saveFigDir='C:\Users\NJ\AppData\Local\Temp\probe95\cohfigs'; batchState.idx = 1; batchState.okRun = true;
run('F:\Work3\CH_NEW_95_射频_DDMA\main_sim.m');
ad = double(squeeze(adcRawData(:, :, 7)));
w = hann(size(ad,1));
Zr = fft(ad .* w, [], 1);
z = Zr(429, :).';
Zf = 20*log10(abs(fft(z .* hann(numel(z)))) + eps);   % z 自身 doppler 谱(Hann)
[~, pk] = maxk(Zf, 8);
m = true(256,1); m(mod((pk-3:pk+3),256)+1) = false;
fprintf('Z2 z谱: 峰 %7.2f | 排线中位 %7.2f | 谷 %7.2f\n', max(Zf), median(Zf(m)), min(Zf));
fprintf('Z2 |z| 周期16? 前32: %s\n', mat2str(round(abs(z(1:32)).',0)));
% 相邻 bin(430)与远处(350)的 z 谱中位对照
z2 = Zr(430, :).'; Zf2 = 20*log10(abs(fft(z2.*hann(256)))+eps);
z3 = Zr(350, :).'; Zf3 = 20*log10(abs(fft(z3.*hann(256)))+eps);
[~, pk2] = maxk(Zf2, 8); m2 = true(256,1); m2(mod((pk2-3:pk2+3),256)+1) = false;
fprintf('Z2 bin430 谱: 峰 %7.2f 排线中位 %7.2f | bin350 中位 %7.2f\n', max(Zf2), median(Zf2(m2)), median(Zf3));
