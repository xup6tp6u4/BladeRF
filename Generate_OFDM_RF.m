fc = 1.8e9;
fspace = 15e3;
dt = 27.65 / fc;
t = 0:dt:1/fspace;
m = 4;
nCarriers = 4;
nc = nCarriers + 1;
co = sqrt(m / 4); % 每象限的點間距
% 定義 N，這裡 N 是空間的維度
N = 4;  % 假設 N = 4 維空間
p = 2;  % 設定範數

% % 所有 I 和 Q 的組合
% I = [1:-1/co:1/co, -1/co:-1/co:-1];
% Q = [1:-1/co:1/co, -1/co:-1/co:-1] * 1j; % 使用 j 代表虛數

% % 使用 meshgrid 生成星座圖
% [Igrid, Qgrid] = meshgrid(I, Q);
% symbolmap = Igrid + Qgrid; % 組合 I 和 Q

% WN 的公式
vol_WN = (2^(N-1) * pi^((N-1)/2)) / gamma((N+1)/2);

% 定義星座點數量和符號生成
symbolmap = zeros(m, N); % 初始化符號矩陣

% 生成符號位於高維球面上
for i = 1:m
    % 計算均勻的角度分佈
    % 使用不同的角度，根據符號的索引 i
    theta = (0:N-1) * (2 * pi * i / m);  % 為每個符號使用不同的角度
    % 計算每個維度
    r = 1;  % 半徑為 1，確保在單位球面上
    am = zeros(1, N);  % 初始化符碼
    for j = 1:N-1
        am(j) = r * cos(theta(j));  % 實部
    end
    am(N) = r * sin(theta(1));  % 將最後一維單獨計算以保持維度
    % 將符號存儲到符號映射中
    symbolmap(i, :) = am;
end
% 確認符號是否不同
disp('symbolmap_final[]:');
for idx = 1:length(symbolmap)
    fprintf('%d: %.4f + %.4fi\n', idx, real(symbolmap(idx)), imag(symbolmap(idx)));
end
%零點遷移
%am_final = (sum(am)) / 2 + 1;
%am_final = 2 + ((2^p)/(m-1)^(p-1));
am_final = 0;%用以測試無零點遷移的波形
symbolmap_final = [symbolmap; am_final * ones(1, N)];  % 添加 am_final 作為一個符號

% 生成 OFDM 波形
symbols = symbolmap_final(:);
zzz = zeros(nCarriers, length(t)); % 初始化 zzz 為一個二維數組
xxx = zeros(nCarriers, length(t));

rfff_add = cos(2 * pi * fc * t) .* sum(zzxx_add_final(:, 1:length(t)), 1); % 修改：只使用前半部分
rfff_mcl = cos(2 * pi * fc * t) .* sum(zzxx_mcl_final, 1);

% MIMO 系統參數
h11 = 1; h12 = 0.5;
h21 = 0.3; h22 = 1.2;
% 選擇 OFDM 的 zzz 子載波作為 MIMO 的輸入
x1 = zzz(1, :); % 使用 zzz 的第 1 子載波作為 x1
x2 = zzz(2, :); % 使用 zzz 的第 2 子載波作為 x2
% 雜訊向量，這裡假設為隨機的雜訊
n1 = randn(size(x1)) * 0.01; 
n2 = randn(size(x2)) * 0.01;
% MIMO的矩陣公式
y1 = h11 * x1 + h12 * x2 + n1;
y2 = h21 * x1 + h22 * x2 + n2;
% 繪製MIMO接收訊號
% figure(1);
% subplot(2, 1, 1);
% plot(real(y1)), title('y1 (接收訊號)');
% xlabel('時間');
% ylabel('幅值');
% subplot(2, 1, 2);
% plot(real(y2)), title('y2 (接收訊號)');
% xlabel('時間');
% ylabel('幅值');
% 確認數值用
%disp(['am_final: ', num2str(am_final), '(零點遷移)']);
disp('symbolmap_final[]:');
for idx = 1:length(symbolmap_final)
    fprintf('%d: %.4f + %.4fi\n', idx, real(symbolmap_final(idx)), imag(symbolmap_final(idx)));
end

% 繪製波形
figure(nCarriers + 1);
subplot(4, 1, 1), plot(t, real(rfff_zzz)), title('rfff_ zzz'), xlabel('Time'), ylabel('Amplitude');
subplot(4, 1, 2), plot(t, real(rfff_xxx)), title('rfff_ xxx'), xlabel('Time'), ylabel('Amplitude');
subplot(4, 1, 3), plot(t, real(rfff_add)), title('rfff_ add'), xlabel('Time'), ylabel('Amplitude');
subplot(4, 1, 4), plot(t, real(rfff_mcl)), title('rfff_ mcl'), xlabel('Time'), ylabel('Amplitude');