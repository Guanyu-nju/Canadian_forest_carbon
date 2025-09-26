function [p_value] = F_test(x, y, num)
% F_TEST Performs F-test for polynomial regression model significance
%   [p_value] = F_test(x, y, num)
%
%   Inputs:
%   x - Independent variable data
%   y - Dependent variable data  
%   num - Degree of polynomial to fit (1 for linear, 2 for quadratic, etc.)
%
%   Output:
%   p_value - p-value from F-test for model significance

    % Fit polynomial model with specified degree
    p = polyfit(x, y, num); % p contains polynomial coefficients
    y_pred = polyval(p, x); % Calculate predicted values
    
    % Calculate R-squared (coefficient of determination)
    SS_total = sum((y - mean(y)).^2); % Total sum of squares
    SS_residual = sum((y - y_pred).^2); % Residual sum of squares
    R2 = 1 - SS_residual / SS_total; % R-squared value
    % fprintf('R²: %.4f\n', R2);
    
    % Calculate F-statistic and p-value
    n = length(y); % Number of data points
    p_params = length(p) - 1; % Number of model parameters (polynomial degree)
    
    % F-statistic calculation
    F = ((SS_total - SS_residual) / p_params) / (SS_residual / (n - p_params - 1));
    
    % p-value calculation using F-distribution
    p_value = 1 - fcdf(F, p_params, n - p_params - 1);
    
    %
    % fprintf('F statistic: %.4f\n', F);
    % fprintf('p-value: %.4f\n', p_value);
    %
    % % Significance test
    % alpha = 0.05; % Significance level
    % if p_value < alpha
    %     disp('Model is significant');
    % else
    %     disp('Model is not significant');
    % end
    
end
