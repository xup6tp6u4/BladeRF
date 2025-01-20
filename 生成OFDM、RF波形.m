% 基本參數設定
fc = 1.8e9;
fspace = 15e3;
dt = 27.65 / fc;
t = 0:dt:1/fspace;
m = 2;  % 每個符碼 m-bit
M = 2^m;
nCarriers = M+1; %M+1個子載波
co = sqrt(M / 4); % 每象限的點間距
N = 4;  % 假設 N = 4 維空間
p = 2;  % 設定範數

% 設定內切圓的半徑及其他參數
syms x;
vol_WN = (2^(N-1) * pi^((N-1)/2)) / gamma((N+1)/2);  % WN 體積
inner_radius = sqrt(vol_WN / pi);  % 根據 WN 體積調整內切圓的半徑

% 隨機生成符合條件的符碼點
symbolmap = zeros(1, nCarriers);  % 存放對應子載波的符碼
while true
    % 隨機生成 [a1, a2, a3] 在 [-inner_radius, inner_radius] 範圍內
    a1 = (rand * 2 - 1) * inner_radius;

    % 檢查點是否位於內切圓內
    if norm([a1, a2, a3], p) <= inner_radius
        % 若點符合條件，生成對應的多項式
        a_poly = x^4 + sqrt(2)/2 * (a1 + a3*1i) * x^3 + a2 * x^2 + sqrt(2)/2 * (a1 - a3*1i) * x + 1;
        disp('生成的多項式:');
        disp(a_poly);
        % 提取多項式係數作為符碼
        coeffs_poly = [1, sqrt(2)/2 * (a1 + a3*1i), a2, sqrt(2)/2 * (a1 - a3*1i), 1];
        symbolmap = coeffs_poly;  % 將多項式係數存入 symbolmap
        break;  % 跳出迴圈
    end
end

% 確認符碼點
disp('符碼 (symbolmap):');
disp(symbolmap);

% 生成 OFDM 波形
OFDM = zeros(nCarriers, length(t));
for n = 1:nCarriers
    OFDM(n, :) = symbolmap(n) .* exp(1i * 2 * pi * fspace * n * t); % 將符碼載入個不同的子載波
end

% 計算每個子載波的 RF 時域信號
RF = zeros(nCarriers, length(t));
for n = 1:nCarriers
    RF(n, :) = cos(2 * pi * fc * t) .* real(OFDM(n, :));  % 每個子載波的 RF 波形
end

% 繪製 OFDM 波形
figure(1);
for n = 1:nCarriers
    subplot(nCarriers, 1, n);  % 每個子載波一個子圖
    plot(t, real(OFDM(n, :)));  % 只繪製 OFDM 的實部
    title(['第 ', num2str(n), ' 個子載波的 OFDM 波形']);
    xlabel('Time (s)');
    ylabel('Amplitude');
    grid on;
end
sgtitle('OFDM 波形');  % 整個圖的標題

% 繪製 RF 波形及合併結果
figure(2);
% 每個子載波的 RF 波形
for n = 1:nCarriers
    subplot(nCarriers + 1, 1, n);  % nCarriers+1 行，將合併波形放在最後
    plot(t, real(RF(n, :)));  % 只繪製 RF 的實部
    title(['第 ', num2str(n), ' 個子載波的 RF 波形']);
    xlabel('Time (s)');
    ylabel('Amplitude');
    grid on;
end
% 合併 RF 波形
RF_combined = sum(RF, 1);  % 合併所有 RF 波形
subplot(nCarriers + 1, 1, nCarriers + 1);  % 最後一個子圖
plot(t, real(RF_combined));
title('所有子載波的 RF 波形合併結果');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% 添加整體標題
sgtitle('RF 波形及合併結果');


% 符碼點的 3D 投影
figure(3);
scatter3(a1, a2, a3, 50, 'filled', 'MarkerFaceColor', 'b');  % 投影點 (a1, a2, a3)
hold on;

% 產生內切圓的 3D 球體
[theta, phi] = meshgrid(linspace(0, pi, 50), linspace(0, 2*pi, 50));  % 球面座標
x_inner = inner_radius * sin(theta) .* cos(phi);
y_inner = inner_radius * sin(theta) .* sin(phi);
z_inner = inner_radius * cos(theta);

% 繪製內切圓的 3D 球體
surf(x_inner, y_inner, z_inner, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor', [1 0 0]);

% 添加圖形標籤和格式
title('取點(a1, a2, a3) 在三維空間內的投影');
xlabel('w1');
ylabel('w2');
zlabel('w3');
axis equal;
grid on;
hold off;