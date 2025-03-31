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
f0 = 80;         % 基頻設定為 80 Hz
M = 100;         % 最大諧波數
m = 90;          % m 用來調整濾波器上截止頻率
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
s_t_shifted = s_t + a_m * M * f0 * exp(1j * 2 * pi * M * f0 * t); 

%% 設計濾波器
f_low = f0 / (Fs/2);       
f_high = ((M-m)*f0) / (Fs/2); 
filter_order = 100;
b = fir1(filter_order, [f_low, f_high], 'bandpass'); 
s_t_shifted = s_t_shifted / max(abs(s_t_shifted));
s_t_filtered = filter(b, 1, real(s_t_shifted));  

%% 使用 STFT 分析濾波後的訊號
win_duration = 0.04;                   % 40 毫秒
win_length = round(win_duration * Fs);   % 窗長（樣本數）
win = hamming(win_length, 'periodic');   % 使用 Hamming 窗
overlap_length = round(0.75 * win_length); % 75% 重疊
nfft = 1024;                           % FFT 長度

% 此處先用非Centered方式進行 tone 提取
[S_nonCentered, F_stft_nonCentered, T_stft] = stft(s_t_filtered, Fs, ...
    'Window', win, 'OverlapLength', overlap_length, 'FFTLength', nfft);

%% 利用 STFT 諧波提取與合成（時域合成法）
% 1. 設定諧波參數：僅考慮 f0, 2f0, ...,(M-m)*f0 的諧波
num_tones = M - m;
tone_freqs = (1:num_tones)' * f0;  % 固定諧波頻率

% 2. STFT 的時間幀數（取自 T_stft）
num_frames = length(T_stft);
hop = win_length - overlap_length;   % 幀移
len_recon = (num_frames-1)*hop + win_length;  % 重組後訊號長度
s_tone_recon = zeros(len_recon, 1);

% 建立矩陣儲存每幀提取的諧波幅度與相位
tone_amplitudes_STFT = zeros(num_tones, num_frames);
tone_phases_STFT = zeros(num_tones, num_frames);

for m_idx = 1:num_frames
    S_frame = S_nonCentered(:, m_idx);
    for n = 1:num_tones
         [~, bin_idx] = min(abs(F_stft_nonCentered - tone_freqs(n)));
         tone_amplitudes_STFT(n, m_idx) = abs(S_frame(bin_idx));
         tone_phases_STFT(n, m_idx) = angle(S_frame(bin_idx));
    end
    tau = (0:win_length-1)'/Fs;
    s_frame = zeros(win_length, 1);
    for n = 1:num_tones
         s_frame = s_frame + 2 * tone_amplitudes_STFT(n, m_idx) * ...
             cos(2*pi*tone_freqs(n)*tau + tone_phases_STFT(n, m_idx));
    end
    s_frame = s_frame .* win;
    start_idx = (m_idx-1)*hop + 1;
    end_idx = start_idx + win_length - 1;
    s_tone_recon(start_idx:end_idx) = s_tone_recon(start_idx:end_idx) + s_frame;
end

% 正規化與增益處理（時域諧波合成結果）
s_tone_recon_norm = s_tone_recon / max(abs(s_tone_recon));
gain_synth = 1;  
s_tone_recon_amp = s_tone_recon_norm * gain_synth;
s_tone_recon_amp(s_tone_recon_amp > 1) = 1;
s_tone_recon_amp(s_tone_recon_amp < -1) = -1;

%% 利用提取的 tone 資訊構造新的 STFT 矩陣並進行 iSTFT
% 使用 'Centered' 的 STFT 參數，以獲得負頻率部分
[S_centered, F_stft_centered, T_stft_centered] = stft(s_t_filtered, Fs, ...
    'Window', win, 'OverlapLength', overlap_length, 'FFTLength', nfft, 'Centered', true);

% 新的 STFT 矩陣初始化 (與 S_centered 同尺寸)
[nFreq, numFrames] = size(S_centered);
S_tones_extracted = zeros(nFreq, numFrames);

% 對於每個提取的 tone，將其填入對應的正負頻率 bin
for n = 1:num_tones
    % 找正頻率 bin
    [~, pos_idx] = min(abs(F_stft_centered - tone_freqs(n)));
    % 找負頻率 bin (對於實信號，應取共軛)
    [~, neg_idx] = min(abs(F_stft_centered + tone_freqs(n)));
    
    for m_idx = 1:numFrames
        tone_val = tone_amplitudes_STFT(n, m_idx) * exp(1j * tone_phases_STFT(n, m_idx));
        S_tones_extracted(pos_idx, m_idx) = tone_val;
        S_tones_extracted(neg_idx, m_idx) = conj(tone_val);
    end
end

% 從新的 STFT 矩陣做 iSTFT還原時域訊號
s_recon_tones = istft(S_tones_extracted, Fs, 'Window', win, 'OverlapLength', overlap_length, 'FFTLength', nfft);
s_recon_tones_norm = s_recon_tones / max(abs(s_recon_tones));
gain_tones = 10;  
s_recon_tones_amp = s_recon_tones_norm * gain_tones;
s_recon_tones_amp(s_recon_tones_amp > 1) = 1;
s_recon_tones_amp(s_recon_tones_amp < -1) = -1;
t_recon_tones = (0:length(s_recon_tones_amp)-1)'/Fs;

%% 計算平均頻譜 (ak')
S_avg = mean(abs(S_centered), 2);

%% 正規化 FFT 結果（以相同基準比較）
norm_factor = max([max(a_k) max(S_avg)]);
a_k = a_k / norm_factor;
S_avg = S_avg / norm_factor;

%% 播放音訊（你可以分別播放比較）
% disp('播放 時域 諧波合成結果 ...');
% soundsc(s_tone_recon_amp, Fs);
% pause(length(s_tone_recon)/Fs + 1);
% 
% disp('播放 提取 tone 後 iSTFT 還原結果 ...');
% soundsc(s_recon_tones_amp, Fs);
% pause(length(s_recon_tones_amp)/Fs + 1);

%% 將諧波合成後的訊號 s_tone_recon_amp 輸出為 WAV 檔
outputFile = 'C:\Users\F112112136\Downloads\s_tone_recon_amp.wav';
audiowrite(outputFile, s_tone_recon_amp, Fs);
disp(['已將合成後的音訊輸出至： ' outputFile]);

%% (過零點檢測) 利用「包絡」判斷波包，然後找包絡低谷當邊界
s_wave = s_tone_recon_amp;
env = abs(hilbert(s_wave)); 
smoothEnv = movmean(env, 100);   % 平滑

% 閾值判斷
threshold = 0.02;
lowEnvIdx = find(smoothEnv < threshold);  % 包絡低於閾值的樣本索引

% 全域過零點
zcs = find(diff(sign(s_wave)) ~= 0);

% 分類
zcs_inGap = [];   % 落在包絡低區
zcs_inBurst = []; % 落在包絡高區
for z = 1:length(zcs)
    zcPos_1 = zcs(z);
    if smoothEnv(zcPos_1) < threshold
        zcs_inGap(end+1) = zcPos_1;
    else
        zcs_inBurst(end+1) = zcPos_1;
    end
end

% 產生與 s_wave 對應的時間向量
t_wave = (0:length(s_wave)-1)'/Fs;

% 繪圖：在原始訊號上標出選取的過零點
figure;

subplot(2,1,1);
plot(t_wave, s_wave, 'b'); hold on;
plot(t_wave(zcs_inBurst), s_wave(zcs_inBurst), 'ro', 'MarkerSize',8, 'LineWidth',1.5);
xlabel('時間 (秒)');
ylabel('振幅');
title('波包內部的過零點');
legend('訊號','波包內部');

subplot(2,1,2);
plot(t_wave, s_wave, 'b'); hold on;
plot(t_wave(zcs_inGap), s_wave(zcs_inGap), 'gs', 'MarkerSize',8, 'LineWidth',1.5);
xlabel('時間 (秒)');
ylabel('振幅');
title('波包間的過零點');
legend('訊號','波包間');

% 這樣 zcs_inGap 就是「波包間」的過零點，zcs_inBurst 就是「波包內部」的過零點


%% 畫圖比較
figure;
subplot(5,1,1);
plot(t, s_t);
title('原始音訊 s(t)');
xlabel('時間 (秒)');
ylabel('振幅');

subplot(5,1,2);
plot(t_recon_tones, real(s_recon_tones_amp));
title('提取tone後 iSTFT 還原的訊號');
xlabel('時間 (秒)');
ylabel('振幅');

subplot(5,1,3);
plot(t_recon_tones, real(s_tone_recon_amp));
title('諧波合成還原的訊號');
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
plot(F_stft_centered, S_avg, 'r');
title('零點遷移後 STFT 雙邊頻譜');
xlabel('頻率 (Hz)');
ylabel('正規化振幅');
xlim([min(F_stft_centered) max(F_stft_centered)]);
ylim([0 max(S_avg)*1.1]);

figure;
subplot(2,1,1);
stem(tone_freqs, tone_amplitudes_STFT(:, m), 'filled');
title('第 m 幀：提取的諧波頻率 振幅');
xlabel('頻率 (Hz)');
ylabel('正規化振幅');
xlim([0 max(tone_freqs)*1.1]);

subplot(2,1,2);
stem(tone_freqs, tone_phases_STFT(:, m), 'filled');
title('第 m 幀：提取的諧波頻率 相位');
xlabel('頻率 (Hz)');
ylabel('弧度');
xlim([0 max(tone_freqs)*1.1]);
toc;