fc = 1.8e9;
fspace = 15e3;
dt = 27.65 / fc;
t = 0:dt:1/fspace;
m = 4;
nCarriers = 4;
nc = nCarriers + 1;
co = sqrt(m / 4); % 每象限的點間距
% 所有 I 和 Q 的組合
I = [1:-1/co:1/co, -1/co:-1/co:-1];
Q = [1:-1/co:1/co, -1/co:-1/co:-1] * 1j; % 使用 j 代表虛數
% 使用 meshgrid 生成星座圖
[Igrid, Qgrid] = meshgrid(I, Q);
symbolmap = Igrid + Qgrid; % 組合 I 和 Q
% 計算 am 符碼
am = abs(symbolmap(:)); % 將 symbolmap 展平成一維數組
am_final = (sum(am)) / 2 + 2;
% 最終的 symbolmap
symbolmap_final = [symbolmap(:); am_final];
% 生成 OFDM 波形
symbols = symbolmap(:);
zzz = zeros(nCarriers, length(t)); % 初始化 zzz 為一個二維數組
xxx = zeros(nCarriers, length(t));
for n = 1:nCarriers
    zzz(n, :) = symbols(n) .* exp(1i * 2 * pi * fspace * n * t); % 第 n 個子載波
    xxx(n, :) = symbols(1) .* exp(1i * 2 * pi * fspace * n * t); % 第 n 個子載波
end
zzxx_add = [zzz, xxx];
zzxx_mcl = zzz .* xxx;
% 生成第 nc 個子載波
zzz_nc = am_final .* exp(1i * 2 * pi * fspace * nc * t); % 第 nc 個子載波
zzz_nc = repmat(zzz_nc, 1, 2); % 重複 zzz_nc 使其長度與 zzxx_add 匹配
% 合併
zzxx_add_final = [zzxx_add; zzz_nc]; % 修改：直接添加 zzz_nc，無需轉置
zzxx_mcl_final = [zzxx_mcl; zzz_nc(:, 1:size(zzxx_mcl, 2))]; % 修改：確保尺寸匹配
% 計算 rfff 和 rffff
rfff_zzz = cos(2 * pi * fc * t) .* sum(zzz, 1); % 合併所有子載波
rfff_xxx = cos(2 * pi * fc * t) .* sum(xxx, 1);
rfff_add = cos(2 * pi * fc * t) .* sum(zzxx_add_final(:, 1:length(t)), 1); % 修改：只使用前半部分
rfff_mcl = cos(2 * pi * fc * t) .* sum(zzxx_mcl_final, 1);
% 確認數值用
disp(['am_final: ', num2str(am_final), '(零點遷移)']);
disp('symbolmap_final[]:');
for idx = 1:length(symbolmap_final)
    fprintf('%d: %.4f + %.4fi\n', idx, real(symbolmap_final(idx)), imag(symbolmap_final(idx)));
end
% 繪製波形
figure(nCarriers + 1);
subplot(4, 1, 1), plot(t, real(rfff_zzz)), title('rfff_zzz'), xlabel('Time'), ylabel('Amplitude');
subplot(4, 1, 2), plot(t, real(rfff_xxx)), title('rfff_xxx'), xlabel('Time'), ylabel('Amplitude');
subplot(4, 1, 3), plot(t, real(rfff_add)), title('rfff_add'), xlabel('Time'), ylabel('Amplitude');
subplot(4, 1, 4), plot(t, real(rfff_mcl)), title('rfff_mcl'), xlabel('Time'), ylabel('Amplitude');