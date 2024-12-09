% 1. 讀取音訊檔案
[audio, Fs] = audioread('C:\Users\IVAN\Downloads\rick astley - never gonna give you up_5sec.wav'); % 替換成你的音訊檔案
audio = audio(:, 1); % 如果是多聲道，取第一聲道

% 在音訊後加入一個零向量
audio = [audio; 0]; % 在末尾增加一個零樣本
N = length(audio); % 更新樣本數
time = (0:N-1) / Fs; % 時間向量

% 播放原始音訊
disp('播放原始音訊...');
sound(audio, Fs);
pause(length(audio)/Fs + 1);

% 2. 對音訊進行 FFT
AudioFFT = fft(audio);

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
plot(time, audio);
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
plot(time, ProcessedAudio);
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