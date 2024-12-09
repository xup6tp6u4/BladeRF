% 1. 讀取音訊檔案
[audio, Fs] = audioread('C:\Users\IVAN\Downloads\rick astley - never gonna give you up_5sec.wav'); % 替換成你的音訊檔案
audio = audio(:, 1); % 如果是多聲道，取第一聲道

% 播放原始音訊
disp('播放原始音訊...');
sound(audio, Fs);
pause(length(audio)/Fs + 1);

% 自定義FFT取樣長度和子載波間隔
N = 240e3; % 設定FFT的取樣長度（可以根據需要調整）
subcarrier_spacing = Fs / N; % 自定義子載波間隔（頻率解析度）

% 在音訊後加入零向量，使其長度為 N 的倍數
audio_padded = [audio; zeros(N - mod(length(audio), N), 1)]; % 填充零，確保長度為 N 的倍數
time_padded = (0:length(audio_padded)-1) / Fs; % 時間向量對應填充後的音訊

% 時間向量只基於原始音訊長度
time_original = (0:length(audio)-1) / Fs; % 原始音訊的時間向量

% 2. 對音訊進行 FFT
AudioFFT = fft(audio_padded, N); % 使用自定義的 N 進行 FFT

% 雙邊頻率範圍
freq = (-N/2:N/2-1) * (Fs / N);

% 使用 fftshift 將頻譜中心化
AudioFFT_shifted = fftshift(AudioFFT);

disp('AudioFFT（帶編號的輸出）：');
for k = 1:length(AudioFFT)
    fprintf('%d: %.6f + %.6fi\n', k, real(AudioFFT(k)), imag(AudioFFT(k)));
end

% 3. 零點遷移
am_final = (sum(AudioFFT)) / 2 + 1; % 計算新增點
AudioFFT(end) = am_final; % 將最後一個值替換為 am_final

% 更新處理後頻譜
ProcessedFFT = fftshift(AudioFFT);

% 4. 使用 iFFT 還原信號
ProcessedAudio = ifft(AudioFFT, 'symmetric'); % 還原信號
ProcessedAudio = ProcessedAudio / (max(abs(ProcessedAudio))); % 正規化
ProcessedAudio = ProcessedAudio * 0.45; % 降低音量至 80%

% 播放處理後的音訊
disp('播放處理後的音訊...');
sound(ProcessedAudio, Fs);
pause(length(ProcessedAudio)/Fs + 1);

% 5. 繪製時域和頻域波形
figure(1);
subplot(2, 2, 1);
plot(time_original, audio);
title('原始音訊 - 時域');
xlabel('時間 (秒)');
ylabel('振幅');

% 原始音訊頻譜（雙邊頻譜，中心化）
subplot(2, 2, 2);
plot(freq, abs(AudioFFT_shifted));
title('原始音訊 - 雙邊頻域');
xlabel('頻率 (Hz)');
ylabel('幅度');

subplot(2, 2, 3);
plot(time_original, ProcessedAudio(1:length(audio))); % 使用原始音訊長度繪圖
title('處理後音訊 - 時域');
xlabel('時間 (秒)');
ylabel('振幅');

% 處理後音訊頻譜（雙邊頻譜，中心化）
subplot(2, 2, 4);
plot(freq, abs(ProcessedFFT));
title('處理後音訊 - 雙邊頻域');
xlabel('頻率 (Hz)');
ylabel('幅度');

% 6. 儲存處理後的音訊（可選）
audiowrite('C:\Users\IVAN\Downloads\processed_audio.wav', ProcessedAudio, Fs);
disp('處理後的音訊已儲存為 processed_audio.wav');
