function [p_value] = F_test(x, y, num)


% 拟合线性模型（1次多项式）
p = polyfit(x, y, num); % p 是拟合参数
y_pred = polyval(p, x); % 计算预测值

% 计算 R²
SS_total = sum((y - mean(y)).^2); % 总平方和
SS_residual = sum((y - y_pred).^2); % 残差平方和
R2 = 1 - SS_residual / SS_total; % 决定系数 R²
% fprintf('R²: %.4f\n', R2);

% 计算 F 统计量和 p 值
n = length(y); % 样本数量
p_params = length(p) - 1; % 模型参数数量（多项式阶数）

F = ((SS_total - SS_residual) / p_params) / (SS_residual / (n - p_params - 1)); % F 统计量
p_value = 1 - fcdf(F, p_params, n - p_params - 1); % p 值
%
% fprintf('F 统计量: %.4f\n', F);
% fprintf('p 值: %.4f\n', p_value);
%
% % 判断显著性
% alpha = 0.05; % 显著性水平
% if p_value < alpha
%     disp('模型显著');
% else
%     disp('模型不显著');
end