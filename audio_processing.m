clc; clear; close all;
tic;
% 讀取音訊文件
file_path = 'C:\Users\F112112136\Downloads\321.wav'; 
[audio, Fs] = audioread(file_path); 
if size(audio, 2) > 1
    audio = mean(audio, 2);
end
s_t = audio / max(abs(audio));
%% 參數設定
T = length(s_t) / Fs;
t = (0:length(s_t)-1)' / Fs;
f0 = 80;         % 基頻設定為 80 Hz
M = 100;         % 最大諧波數
m = 10;          % m 用來調整濾波器上截止頻率
M_max = floor((Fs/2) / f0);
M = min(M, M_max);
N = length(s_t); 
f = (-N/2:N/2-1) * (Fs/N);  % 雙邊頻譜頻率軸
%% 計算 (傳統 FFT)
S_f_shifted = fftshift(fft(s_t, N)); 
a_k = abs(S_f_shifted);  % 原始 FFT 振幅
%% 計算 a_m（使用平均能量）
a_m = (1/2) * sum(abs(s_t)) / N; 
%% 非傳統零點遷移
% 這裡我們使用複數形式，但取實部後再進行濾波處理
s_t_shifted = s_t + a_m * M * f0 * exp(1j * 2 * pi * M * f0 * t); 
%% 設計濾波器
% 濾波器通帶設定為 [f0, ((M-m)*f0)]
f_low = f0 / (Fs/2);       
f_high = ((M-m)*f0) / (Fs/2); 
filter_order = 100;
b = fir1(filter_order, [f_low, f_high], 'bandpass'); 
% 為了使用 FIR 濾波器處理實數訊號，取 s_t_shifted 的實部
% 濾波前也可先正規化 s_t_shifted 避免數值過大
s_t_shifted = s_t_shifted / max(abs(s_t_shifted));
s_t_filtered = filter(b, 1, real(s_t_shifted));  

%% 使用 STFT 分析濾波後的訊號
% 設定窗長、重疊長度與 FFT 長度
win_duration = 0.04;                   % 40 毫秒
win_length = round(win_duration * Fs);   % 轉換成樣本數
win = hamming(win_length, 'periodic');   % 使用 Hamming 窗
overlap_length = round(0.75 * win_length); % 75% 重疊
nfft = 1024;                           % FFT 長度

% 使用 MATLAB 內建 stft 函數 (MATLAB R2019b 以上版本)
[S, F_stft, T_stft] = stft(s_t_filtered, Fs, ...
    'Window', win, 'OverlapLength', overlap_length, 'FFTLength', nfft);

%% 利用 STFT 諧波提取與合成
% 1. 設定諧波參數：我們只考慮 f0, 2f0, ...,(M-m)*f0 的諧波
num_tones = M - m;
tone_freqs = (1:num_tones)' * f0;  % 固定諧波頻率

% 2. STFT 的時間幀數
num_frames = length(T_stft);
hop = win_length - overlap_length;   % 幀移
len_recon = (num_frames-1)*hop + win_length;  % 重組後訊號長度
s_tone_recon = zeros(len_recon, 1);

% 建立矩陣儲存每幀提取的諧波幅度與相位（可用於觀察）
tone_amplitudes_STFT = zeros(num_tones, num_frames);
tone_phases_STFT = zeros(num_tones, num_frames);

% 3. 逐幀提取與合成
for m_idx = 1:num_frames
    % 取出第 m_idx 幀的 STFT結果（頻域資訊）
    S_frame = S(:, m_idx);
    % 針對每個諧波 tone 提取幅度與相位
    for n = 1:num_tones
         % 找出 STFT 頻率向量中最接近 n*f0 的頻率 bin
         [~, bin_idx] = min(abs(F_stft - tone_freqs(n)));
         tone_amplitudes_STFT(n, m_idx) = abs(S_frame(bin_idx));
         tone_phases_STFT(n, m_idx) = angle(S_frame(bin_idx));
    end
    % 合成第 m_idx 幀的時域信號
    % 產生一個單幀的時間向量（相對於幀內）
    tau = (0:win_length-1)'/Fs;
    s_frame = zeros(win_length, 1);
    for n = 1:num_tones
         % 以 2 倍幅度補償正負頻率（假設 STFT 只取正頻率）
         s_frame = s_frame + 2 * tone_amplitudes_STFT(n, m_idx) * ...
             cos(2*pi*tone_freqs(n)*tau + tone_phases_STFT(n, m_idx));
    end
    % 可選：乘上與分析窗相同的合成窗（此處用 win）
    s_frame = s_frame .* win;
    
    % 將合成的單幀訊號進行 Overlap-Add 重組
    start_idx = (m_idx-1)*hop + 1;
    end_idx = start_idx + win_length - 1;
    s_tone_recon(start_idx:end_idx) = s_tone_recon(start_idx:end_idx) + s_frame;
end

% 4. 對重組後的諧波合成訊號做正規化、增益與夾幅處理
s_tone_recon_norm = s_tone_recon / max(abs(s_tone_recon));
gain_synth = 5;  % 可根據需求調整增益值
s_tone_recon_amp = s_tone_recon_norm * gain_synth;
s_tone_recon_amp(s_tone_recon_amp > 1) = 1;
s_tone_recon_amp(s_tone_recon_amp < -1) = -1;

%% 計算平均頻譜 (ak')
S_avg = mean(abs(S), 2);

%% 設定 STFT 為 Centered，為了獲得負頻率部分
[S, F_stft, T_stft] = stft(s_t_filtered, Fs, 'Window', win, ...
    'OverlapLength', overlap_length, 'FFTLength', nfft, 'Centered', true);

%% 正規化 FFT 結果（以相同基準比較）
norm_factor = max([max(a_k) max(S_avg)]);
a_k = a_k / norm_factor;
S_avg = S_avg / norm_factor;

%% 計算 iFFT (istft)
s_reconstructed_STFT = istft(S, Fs, 'Window', win, 'OverlapLength', overlap_length, 'FFTLength', nfft);
s_recon_STFT_norm = s_reconstructed_STFT / max(abs(s_reconstructed_STFT));
gain_synth = 80;  % 可根據需求調整增益值
s_recon_STFT_amp = s_recon_STFT_norm * gain_synth;
s_recon_STFT_amp(s_recon_STFT_amp > 1) = 1;
s_recon_STFT_amp(s_recon_STFT_amp < -1) = -1;

t_recon = (0:length(s_recon_STFT_amp)-1)'/Fs;

%% 播放音訊
soundsc(s_tone_recon_amp, Fs);
pause(length(s_tone_recon)/Fs + 1);

%% 將諧波合成後的訊號 s_tone_recon_amp 輸出為 WAV 檔
outputFile = 'C:\Users\F112112136\Downloads\s_tone_recon_amp.wav';
audiowrite(outputFile, s_tone_recon_amp, Fs);

disp(['已將合成後的音訊輸出至： ' outputFile]);

%% 畫圖
figure;
subplot(5,1,1);
plot(t, s_t);
title('原始音訊 s(t)');
xlabel('時間 (秒)');
ylabel('振幅');

subplot(5,1,2);
plot(t_recon, real(s_recon_STFT_amp));
title('零點遷移後iFFT的時域信號');
xlabel('時間 (秒)');
ylabel('振幅');

subplot(5,1,3);
plot(t_recon, real(s_tone_recon_amp));
title('諧波合成還原的時域信號');
xlabel('時間 (秒)');
ylabel('振幅');

subplot(5,1,4);
plot(f, a_k, 'b');
title('直接對原始音訊做 FFT 雙邊頻譜');
xlabel('頻率 (Hz)');
ylabel('正規化振幅');
xlim([min(f) max(f)]);
ylim([0 max(a_k)*1.1]);

subplot(5,1,5);
plot(F_stft, S_avg, 'r');
title('零點遷移後 STFT 雙邊頻譜');
xlabel('頻率 (Hz)');
ylabel('正規化振幅');
xlim([min(F_stft) max(F_stft)]);
ylim([0 max(S_avg)*1.1]);

figure;
subplot(2,1,1);
stem(tone_freqs, tone_amplitudes_STFT(:, m), 'filled');
title('提取的諧波頻率 振幅 (Tone Frequencies)');
xlabel('頻率 (Hz)');
ylabel('正規化振幅');
xlim([0 max(tone_freqs)*1.1]);
subplot(2,1,2);
stem(tone_freqs, tone_phases_STFT(:, m), 'filled');
title('提取的諧波頻率 相位 (Tone Frequencies)');
xlabel('頻率 (Hz)');
ylabel('正規化振幅');
xlim([0 max(tone_freqs)*1.1]);
toc;