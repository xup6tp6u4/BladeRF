clc; clear; close all;
tic;
%% 讀取音訊文件
file_path = 'C:\Users\F112112136\Downloads\321.wav'; 
[audio, Fs] = audioread(file_path); 
if size(audio, 2) > 1
    audio = mean(audio, 2);
end
s_t = audio / max(abs(audio));
%% 參數設定
T = length(s_t) / Fs;
t = (0:length(s_t)-1)' / Fs;
f0 = 80;         
M = 100;
m = 10;
M_max = floor((Fs/2) / f0);
M = min(M, M_max);
N = length(s_t); 
f = (-N/2:N/2-1) * (Fs/N); 
%% 計算 (傳統FFT)
S_f_shifted = fftshift(fft(s_t, N)); 
a_k = abs(S_f_shifted); 
%% 計算 a_m
a_m = (1/2) * sum(abs(s_t)) / N; 
%% 非傳統零點遷移
s_t_shifted = s_t + a_m * M * f0 * exp(1j * 2 * pi * M * f0 * t); 
%% 設計濾波器
f_low = f0 / (Fs/2);       % 最低頻率 (正規化)
f_high = ((M-m)*f0) / (Fs/2); % 最高頻率 (正規化)
filter_order = 100;  % 設定濾波器階數
b = fir1(filter_order, [f_low, f_high], 'bandpass'); % 設計帶通濾波器
s_t_filtered = filter(b, 1, real(s_t_shifted));  % 濾波處理
%% 計算 FFT (零點遷移)
S_f_prime = fft(s_t_filtered, N); 
S_f_prime_shifted = fftshift(S_f_prime);
a_k_prime = abs(S_f_prime_shifted); 
%% 播放濾波後的時域波形
filtered_audio = real(s_t_filtered);             % 取實部
filtered_audio = filtered_audio / max(abs(filtered_audio)); % 適度正規化
gain = 80;  % 增益因子，可以根據需要調整
amplified_audio = filtered_audio * gain; % 將超出 [-1, 1] 的值裁剪
amplified_audio(amplified_audio > 1) = 1;
amplified_audio(amplified_audio < -1) = -1;
soundsc(amplified_audio, Fs);     % 播放音訊
pause(length(filtered_audio)/Fs + 1);    % 暫停，確保音訊能播放完
%% 畫圖
figure;
subplot(4,1,1);
plot(t, s_t);
title('原始音訊');
xlabel('時間 (秒)');
ylabel('振幅');

subplot(4,1,2);
plot(t, real(s_t_shifted));
title('零點遷移的時域信號 s''(t) (實部)');
xlabel('時間 (秒)');
ylabel('振幅');

subplot(4,1,3);
plot(f, a_k, 'b');
title('原始 FFT 雙邊頻譜');
xlabel('頻率 (Hz)');
ylabel('正規化振幅');
xlim([min(f) max(f)]);
ylim([0 max(a_k)*1.1]);

subplot(4,1,4);
plot(f, a_k_prime, 'r');
title('零點遷移 雙邊頻譜');
xlabel('頻率 (Hz)');
ylabel('正規化振幅');
xlim([min(f) max(f)]);
ylim([0 max(a_k_prime)*1.1]);
toc;