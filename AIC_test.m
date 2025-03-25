% 示例数据
x = TEM_summer_list(:);
y = GCAS_ER_summer_list(:); % 线性关系加噪声


% 线性拟合
p1 = polyfit(x, y, 1);
y1 = polyval(p1, x);
residuals1 = y - y1;
SSE1 = sum(residuals1.^2);

% 二次多项式拟合
p2 = polyfit(x, y, 2);
y2 = polyval(p2, x);
residuals2 = y - y2;
SSE2 = sum(residuals2.^2);

% 计算 AIC
n = length(y);
k1 = length(p1);
k2 = length(p2);
AIC1 = n * log(SSE1/n) + 2 * k1;
AIC2 = n * log(SSE2/n) + 2 * k2;

% 输出结果
fprintf('线性拟合的 AIC: %.4f\n', AIC1);
fprintf('二次多项式拟合的 AIC: %.4f\n', AIC2);

% 比较 AIC
if AIC1 < AIC2
    disp('线性拟合更好');
else
    disp('二次多项式拟合更好');
end