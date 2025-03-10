% 指定CSV檔案的路徑 
filename = 'D:\WaveData_SA_black.csv';  % 替換為你的實際路徑
filename = 'D:\WaveData_SA_white.csv';  % 攝影機全白
filename = 'D:\WaveData_SP.csv';  % 攝影機條紋圖案(水平)
filename = 'D:\WaveData_TZ.csv';  % 攝影機條紋圖案(垂直)
filename = 'D:\WaveData_SA_QR.csv';  % 攝影機QR Code
filename = 'C:\Program Files\bladeRF\2.43G_baseband.csv';  

% 讀取整個CSV檔案
data = readtable(filename, 'ReadVariableNames', false);

% 提取Wave Data (Delta Hz)
wave_data_neg = data{22:end, 1:2};  % 假設從第22行開始是Wave Data，並且取前兩列數據

% 提取時間點
time_points = (1:length(wave_data_neg))';  % 使用行號作為時間點

% 提取頻率偏移數據
delta_freq1 = wave_data_neg(:, 1); % 第一列頻率偏移數據
delta_freq2 = wave_data_neg(:, 2); % 第二列頻率偏移數據

% 創建一個新圖表窗口並設置兩個子圖
figure;  % 創建新圖表窗口

% 繪製第一個子圖（Delta Frequency 1）
subplot(2, 1, 1);  % 2行1列的子圖佈局，這是第一個子圖
plot(time_points, delta_freq1, '-o', 'DisplayName', 'Delta Frequency 1');
title('Delta Frequency 1 vs Time');
xlabel('Time Points');
ylabel('Frequency Offset (Hz)');
legend('show');
grid on;

% 繪製第二個子圖（Delta Frequency 2）
subplot(2, 1, 2);  % 2行1列的子圖佈局，這是第二個子圖
plot(time_points, delta_freq2, '-x','Color', 'r', 'DisplayName', 'Delta Frequency 2');
title('Delta Frequency 2 vs Time');
xlabel('Time Points');
ylabel('Frequency Offset (Hz)');
legend('show');
grid on;
