%基本參數設定
fc = 1.8e9;
fspace = 15e3;
dt = 27.65 / fc;
t = 0:dt:1/fspace;
m = 2;  % 每個符碼 m-bit
M = 2^m;
nCarriers = M + 1; % M+1 個子載波
N = 4;  % 假設 N = 4 維空間
p = 2;  % 設定範數

% 設定內切圓的半徑及其他參數
syms x;
vol_WN = (2^(N-1) * pi^((N-1)/2)) / gamma((N+1)/2);  % W_N 體積
inner_radius = sqrt(vol_WN / pi);  % 根據 W_N 體積調整內切圓的半徑

% 隨機生成符合條件的符碼點
n_points = 1;  % 生成的點數量，這裡只需要一個點
a1 = zeros(1, n_points);
a2 = zeros(1, n_points);
a3 = zeros(1, n_points);
for i = 1:n_points
    % 隨機生成極座標
    r = inner_radius * (rand^(1/N)); % 隨機半徑 (均勻分布)
    theta1 = rand * pi;              % 隨機極角 (範圍: [0, π])
    theta2 = rand * 2 * pi;          % 隨機方位角 (範圍: [0, 2π])
    % 將球座標轉換為笛卡兒座標
    a1(i) = r * sin(theta1) * cos(theta2);
    a2(i) = r * sin(theta1) * sin(theta2);
    a3(i) = r * cos(theta1);
end

% 檢查生成的點是否正確位於內切圓內 (可選)
assert(all(sqrt(a1.^2 + a2.^2 + a3.^2) <= inner_radius), '生成的點不在內切圓內');
% 生成對應的多項式
a_poly = x^4 + sqrt(2)/2 * (a1 + a3*1i) * x^3 + a2 * x^2 + sqrt(2)/2 * (a1 - a3*1i) * x + 1;
disp('生成的多項式:');
disp(a_poly);
% 提取多項式係數作為符碼
coeffs_poly = [1, sqrt(2)/2 * (a1 + a3*1i), a2, sqrt(2)/2 * (a1 - a3*1i), 1];
symbolmap = coeffs_poly;  % 將多項式係數存入 symbolmap

% 生成 OFDM 波形
OFDM = zeros(nCarriers, length(t));
for n = 1:nCarriers
    OFDM(n, :) = symbolmap(n) .* exp(1i * 2 * pi * fspace * n * t); % 將符碼載入個不同的子載波
end

% 合併 OFDM 波形
OFDM_combined = sum(OFDM, 1);
% 計算每個子載波的 RF 時域信號
RF = zeros(nCarriers, length(t));
for n = 1:nCarriers
    RF(n, :) = cos(2 * pi * fc * t) .* real(OFDM(n, :));  % 每個子載波的 RF 波形
end

% 僅選擇計算子載波範圍 2 ~ (nCarriers-1) 的信號
OFDM_combined_partial = sum(OFDM(2:nCarriers, :), 1);
% OFDM 合併後的 RF 波形 (僅限選定子載波範圍)
RF_combined = cos(2 * pi * fc * t) .* OFDM_combined_partial;
% 繪製 OFDM 波形
figure(1);
% 設定 Y 軸的範圍
y_limits = [-6, 6];  % 根據需要調整範圍
for n = 1:nCarriers
    subplot(nCarriers + 1, 1, n);  % 每個子載波一個子圖
    plot(t, real(OFDM(n, :)));  % 只繪製 OFDM 的實部
    title(['第 ', num2str(n), ' 個子載波的 OFDM 波形']);
    xlabel('Time (s)');
    ylabel('Amplitude');
    grid on;
    ylim(y_limits);  % 設定統一的 Y 軸範圍
end
subplot(nCarriers + 1, 1, nCarriers + 1);  % 最後一個子圖
plot(t, real(OFDM_combined));
title('所有子載波的 OFDM 波形合併結果');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;
ylim(y_limits);  % 設定統一的 Y 軸範圍
sgtitle('OFDM 波形');  % 整個圖的標題

% 繪製 RF 波形及合併結果
figure(2);
% 設定 Y 軸的範圍
y_limits_2 = [-6, 6];  % 根據需要調整範圍
% 每個子載波的 RF 波形
for n = 1:nCarriers
    subplot(nCarriers + 1, 1, n);  % nCarriers+1 行，將合併波形放在最後
    xlabel('Time (s)');
    ylabel('Amplitude');
    grid on;
    ylim(y_limits_2);  % 設定統一的 Y 軸範圍
end
subplot(nCarriers + 1, 1, nCarriers + 1);  % 最後一個子圖
plot(t, real(RF_combined));
title('所有子載波的 RF 波形合併結果');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;
ylim(y_limits_2);  % 設定統一的 Y 軸範圍
sgtitle('RF 波形及合併結果');  % 添加整體標題

% 符碼點的 3D 投影
figure(3);
scatter3(a1, a2, a3, 50, 'filled', 'MarkerFaceColor', 'b');  % 投影點 (a1, a2, a3)
hold on;

% 產生內切圓的 3D 球體
[theta, phi] = meshgrid(linspace(0, pi, 50), linspace(0, 2 * pi, 50));  % 球面座標
x_inner = inner_radius * sin(theta) .* cos(phi);
y_inner = inner_radius * sin(theta) .* sin(phi);
z_inner = inner_radius * cos(theta);
surf(x_inner, y_inner, z_inner, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor', [1 0 0]);

% 添加圖形標籤和格式
title('取點(a1, a2, a3) + WN的三維空間投影');
xlabel('w1');
ylabel('w2');
zlabel('w3');