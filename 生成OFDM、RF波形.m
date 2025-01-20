% 基本參數設定
fc = 1.8e9;
fspace = 15e3;
dt = 27.65 / fc;
t = 0:dt:1/fspace;
m = 2;  % 每個符碼 m-bit
M = 2^m;
nCarriers = 4;
nc = nCarriers + 1;
co = sqrt(M / 4); % 每象限的點間距
N = 4;  % 假設 N = 4 維空間
p = 2;  % 設定範數

% 設定內切圓的半徑及其他參數
syms x;  % 定義符號變量 x

vol_WN = (2^(N-1) * pi^((N-1)/2)) / gamma((N+1)/2);  % WN 體積
inner_radius = sqrt(vol_WN / pi);  % 根據 WN 體積調整內切圓的半徑

% 隨機生成符合條件的符碼點
symbolmap = zeros(1, N); 
while true
    % 隨機生成 [a1, a2, a3] 在 [-inner_radius, inner_radius] 範圍內
    a1 = (rand * 2 - 1) * inner_radius;
    a2 = (rand * 2 - 1) * inner_radius;
    a3 = (rand * 2 - 1) * inner_radius;

    % 檢查點是否位於內切圓內
    if norm([a1, a2, a3], p) <= inner_radius
        symbolmap(1, :) = [a1, a2, a3, 0]; % 4D 空間，最後一維設為 0
        break;  % 符碼點符合條件，跳出迴圈
    end
    % 如果不在內切圓內，會重新生成點
end

% 使用選定的點生成多項式
a = symbolmap(1, 1:3); % 提取 [a1, a2, a3]
a_poly = x^4 + sqrt(2)/2 * (a(1) + a(3)*1i) * x^3 + a(2) * x^2 + sqrt(2)/2 * (a(1) - a(3)*1i) * x + 1;
disp('生成的多項式:');
disp(a_poly);

symbolmap_final = symbolmap;

% 生成 OFDM 波形
symbols = symbolmap_final(:);
zzz = zeros(nCarriers, length(t)); % 初始化 zzz 為一個二維數組
xxx = zeros(nCarriers, length(t));

for n = 1:nCarriers
    nn = n + 1;
    zzz(n, :) = symbols(n) .* exp(1i * 2 * pi * fspace * n * t); % 第 n 個子載波
    xxx(n, :) = symbols(1) .* exp(1i * 2 * pi * fspace * n * t); % 第 n 個子載波
end

zzxx_mcl = zzz .* xxx;

% 生成第 nc 個子載波
zzz_nc = symbolmap(1,4) .* exp(1i * 2 * pi * fspace * nc * t); % 第 nc 個子載波
zzz_nc = repmat(zzz_nc, 1, 2); % 重複 zzz_nc 使其長度與 zzxx_add 匹配

% 合併

rfff_add = cos(2 * pi * fc * t) .* sum(zzxx_add_final(:, 1:length(t)), 1); % 修改：只使用前半部分
rfff_mcl = cos(2 * pi * fc * t) .* sum(zzxx_mcl_final, 1);

% 確認數值用
disp('symbolmap_final[]:');
for idx = 1:length(symbolmap_final)
    fprintf('%d: %.4f + %.4fi\n', idx, real(symbolmap_final(idx)), imag(symbolmap_final(idx)));

subplot(4, 1, 1), plot(t, real(rfff_zzz)), title('rfff_ zzz'), xlabel('Time'), ylabel('Amplitude');
subplot(4, 1, 2), plot(t, real(rfff_xxx)), title('rfff_ xxx'), xlabel('Time'), ylabel('Amplitude');
subplot(4, 1, 3), plot(t, real(rfff_add)), title('rfff_ add'), xlabel('Time'), ylabel('Amplitude');
subplot(4, 1, 4), plot(t, real(rfff_mcl)), title('rfff_ mcl'), xlabel('Time'), ylabel('Amplitude');
% 使用 PCA 將 4 維符碼點降到 3 維
symbolmap_final_3D = symbolmap_final(:, 1:3); % 取前三維進行三維投影
% 繪製符碼點在三維空間的投影
figure;
scatter3(symbolmap_final_3D(:,1), symbolmap_final_3D(:,2), symbolmap_final_3D(:,3), 50, 'filled');
hold on;
% 產生內切圓的 3D 球體
[theta, phi] = meshgrid(linspace(0, pi, 50), linspace(0, 2*pi, 50));  % 球面座標
x_inner = inner_radius * sin(theta) .* cos(phi);  % X 座標
y_inner = inner_radius * sin(theta) .* sin(phi);  % Y 座標
z_inner = inner_radius * cos(theta);  % Z 座標
% 繪製內切圓的 3D 球體
surf(x_inner, y_inner, z_inner, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor', [1 0 0]);
% 添加圖形標籤和格式
title('四維符碼點的三維投影，含 W4 體積的內切圓');
xlabel('X');
ylabel('Y');
zlabel('Z');
axis equal;
grid on;
hold off;