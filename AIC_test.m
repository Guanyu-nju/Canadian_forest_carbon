% Sample data
x = TEM_summer_list(:);
y = GCAS_ER_summer_list(:); % Linear relationship with noise

% Linear fitting (1st degree polynomial)
p1 = polyfit(x, y, 1);
y1 = polyval(p1, x);
residuals1 = y - y1;
SSE1 = sum(residuals1.^2);  % Sum of squared errors

% Quadratic polynomial fitting (2nd degree polynomial)
p2 = polyfit(x, y, 2);
y2 = polyval(p2, x);
residuals2 = y - y2;
SSE2 = sum(residuals2.^2);  % Sum of squared errors

% Calculate AIC (Akaike Information Criterion)
n = length(y);  % Number of data points
k1 = length(p1);  % Number of parameters in linear model (2 parameters: slope and intercept)
k2 = length(p2);  % Number of parameters in quadratic model (3 parameters)
AIC1 = n * log(SSE1/n) + 2 * k1;  % AIC for linear model
AIC2 = n * log(SSE2/n) + 2 * k2;  % AIC for quadratic model

% Output results
fprintf('AIC for linear fitting: %.4f\n', AIC1);
fprintf('AIC for quadratic polynomial fitting: %.4f\n', AIC2);

% Compare AIC values (lower AIC indicates better model)
if AIC1 < AIC2
    disp('Linear fitting is better');
else
    disp('Quadratic polynomial fitting is better');
end
