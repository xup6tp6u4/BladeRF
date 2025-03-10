clear;
clc;

%% 參數設定
fc = 100e3;     % 降低載波頻率到 100 kHz，便於時間域觀察
fm = 1e3;       % 調變訊號頻率 1 kHz
delta_f = 20e3;  % 頻率偏移量 (±20 kHz)
fs = 40e6;       % 取樣率：40 MHz (相對 100 kHz 載波已足夠高)
T = 3/fm;       % 總模擬時間：這裡設 3 個調變週期 (3 ms)
t = 0:1/fs:T;   % 時間向量

%% 產生調變訊號 (1 kHz 正弦波)
m = sin(2*pi*fm*t);

%% 計算調變訊號的積分 (近似積分)
% 使用 cumtrapz 進行數值積分
int_m = cumtrapz(t, m);

%% 產生 FM 訊號
% FM 訊號公式：s(t) = cos(2πfc t + 2πΔf * ∫m(t)dt)
s_fm = cos(2*pi*fc*t + 2*pi*delta_f*int_m);

%% 繪圖檢視
figure;
subplot(2,1,1);
plot(t, m);
title('調變訊號 m(t) (1 kHz 正弦波)');
xlabel('時間 (s)');
ylabel('幅值');


subplot(2,1,2);
plot(t, s_fm);
title('FM 訊號 s_{FM}(t)');
xlabel('時間 (s)');
ylabel('幅值');
grid on;


% 輸出波形資料成 CSV 檔案
% 將時間、波形實部與虛部組合成一個矩陣，每一列分別為一個時間點的數值
waveform_data = [real(s_fm(:)), imag(s_fm(:))];
csv_filename = 'D:\FM_1kHz_Generate.csv';
writematrix(waveform_data, csv_filename); % MATLAB R2019a 以上版本建議使用 writematrix
disp(['波形資料已輸出至 ', csv_filename]);