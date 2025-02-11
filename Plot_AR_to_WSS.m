clear;
close all;

% Load the .mat files for type1, type2, and type3
data1 = load('type1.mat');
data2 = load('type2.mat');
data3 = load('type3.mat');
data4 = load('type4.mat');
Type1 = data1.Type1;
Type2 = data2.Type2;
Type3 = data3.Type3;
Type4 = data4.Type4;
%Column 1: Aspect Ratio (AR)
%Column 2: Low WSS Region (Low Wall Shear Stress Region)
%Column 3: High WSS Region (High Wall Shear Stress Region)
%Column 4: Average WSS (Average Wall Shear Stress)
%Column 5: High OSI Region (High Oscillatory Shear Index Region)
%Column 6: Average OSI (Average Oscillatory Shear Index)
%Column 7: High RRT Region (High Recirculation Residence Time Region)
%Column 8: Average RRT (Average Recirculation Residence Time)
%Column 9: Area Ratio (New x-axis variable for your latest plot)
%Column 10: Volume Normalized Vorticity
%Column 11: Volume
%Column 12: Surface Area
%Column 13: qcrit vortex

% Define the fit type: 'linear' or 'power'
fit_type = 'linear';  % Change this to 'power' to use the power fit

% --- Low WSS Region vs Aspect Ratio ---
figure;

% Extract, filter, and sort relevant columns for low WSS
[AR_filtered1, low_WSS_region_filtered1] = process_data(Type1, 1, 2);
[AR_filtered2, low_WSS_region_filtered2] = process_data(Type2, 1, 2);
[AR_filtered3, low_WSS_region_filtered3] = process_data(Type3, 1, 2);

% Combine data for fitting
AR_combined = [AR_filtered1; AR_filtered2; AR_filtered3];
low_WSS_combined = [low_WSS_region_filtered1; low_WSS_region_filtered2; low_WSS_region_filtered3];

% Plot data points
hold on;
plot(AR_filtered1, low_WSS_region_filtered1, 'o', 'Color', 'b', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Type 1');
plot(AR_filtered2, low_WSS_region_filtered2, '^', 'Color', 'r', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Type 2');
plot(AR_filtered3, low_WSS_region_filtered3, 's', 'Color', 'm', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Type 3');

% Generate x-values for fit line
AR_fit = linspace(min(AR_combined), max(AR_combined), 100);

% Fit and plot based on fit_type
if strcmp(fit_type, 'linear')
    % Linear fit
    p_linear = polyfit(AR_combined, low_WSS_combined, 1);
    yfit = polyval(p_linear, AR_fit);
    plot(AR_fit, yfit, '--k', 'LineWidth', 2, 'DisplayName', 'Linear Fit');

    % Calculate R-squared for linear fit
    yfit_combined = polyval(p_linear, AR_combined);
    SS_res = sum((low_WSS_combined - yfit_combined).^2);
    SS_tot = sum((low_WSS_combined - mean(low_WSS_combined)).^2);
    R_squared = 1 - (SS_res / SS_tot);

elseif strcmp(fit_type, 'power')
    % Power fit
    log_AR_combined = log(AR_combined);
    log_low_WSS_combined = log(low_WSS_combined);
    p_power = polyfit(log_AR_combined, log_low_WSS_combined, 1);
    a_power = exp(p_power(2));
    b_power = p_power(1);

    yfit = a_power * AR_fit.^b_power;
    plot(AR_fit, yfit, '--k', 'LineWidth', 2, 'DisplayName', 'Power Fit');

    % Calculate R-squared for power fit
    yfit_combined = a_power * AR_combined.^b_power;
    SS_res = sum((low_WSS_combined - yfit_combined).^2);
    SS_tot = sum((low_WSS_combined - mean(low_WSS_combined)).^2);
    R_squared = 1 - (SS_res / SS_tot);
else
    error('Invalid fit_type. Choose ''linear'' or ''power''.');
end

% Display R-squared value
text(mean(AR_combined), max(low_WSS_combined), ['R^2 = ', num2str(R_squared, 2)], 'FontSize', 12, 'Color', 'k');

xlabel('Aspect Ratio');
ylabel('Low WSS Region');
title(['Low WSS Region vs Aspect Ratio (' fit_type ' fit)']);
legend('Location', 'best');
grid off;
AR_filtered4_row1 = Type4(1, 1);  % Aspect Ratio for row 1
low_WSS_filtered4_row1 = Type4(1, 2);  % Low WSS Region for row 1

AR_filtered4_row2 = Type4(2, 1);  % Aspect Ratio for row 2
low_WSS_filtered4_row2 = Type4(2, 2);  % Low WSS Region for row 2

% Plot the first row point with a unique marker and color
plot(AR_filtered4_row1, low_WSS_filtered4_row1, 'h', 'Color', 'c', ...
     'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Type 4 (ANY_007)');

% Plot the second row point with a unique marker and color
plot(AR_filtered4_row2, low_WSS_filtered4_row2, '*', 'Color', 'k', ...
     'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Type 4 (ANY_111)');
% --- Distance Calculation for Type4 Points ---
% Fit line parameters
m = p_linear(1); % Slope of the fit line
c = p_linear(2); % Intercept of the fit line

% Type4 Data Points
x4_row1 = AR_filtered4_row1; % Surface Area for ANY_007
y4_row1 = low_WSS_filtered4_row1; % Low WSS Region for ANY_007

x4_row2 = AR_filtered4_row2; % Surface Area for ANY_111
y4_row2 = low_WSS_filtered4_row2; % Low WSS Region for ANY_111

% Predicted y-values on the fit line for Type4 points
yfit_row1 = m * x4_row1 + c;
yfit_row2 = m * x4_row2 + c;

% Calculate the vertical distances (residuals)
distance_row1 = abs(y4_row1 - yfit_row1);
distance_row2 = abs(y4_row2 - yfit_row2);

% Display the distances
fprintf('low_WSS vs aspect ratio Distance of ANY_007 to the fit line: %.4f\n', distance_row1);
fprintf('low_WSS vs aspect ratio Distance of ANY_111 to the fit line: %.4f\n', distance_row2);

% Optionally annotate the plot with the distances
text(x4_row1, y4_row1, sprintf('%.4f', distance_row1), 'FontSize', 10, 'Color', 'c', 'VerticalAlignment', 'bottom');
text(x4_row2, y4_row2, sprintf('%.4f', distance_row2), 'FontSize', 10, 'Color', 'k', 'VerticalAlignment', 'bottom');
% --- High WSS Region vs Aspect Ratio ---
figure;

% Extract, filter, and sort relevant columns for high WSS
[AR_filtered1, high_WSS_region_filtered1] = process_data(Type1, 1, 3);
[AR_filtered2, high_WSS_region_filtered2] = process_data(Type2, 1, 3);
[AR_filtered3, high_WSS_region_filtered3] = process_data(Type3, 1, 3);

% Combine data for fitting
AR_combined = [AR_filtered1; AR_filtered2; AR_filtered3];
high_WSS_combined = [high_WSS_region_filtered1; high_WSS_region_filtered2; high_WSS_region_filtered3];

% Plot data points
hold on;
plot(AR_filtered1, high_WSS_region_filtered1, 'o', 'Color', 'b', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Type 1');
plot(AR_filtered2, high_WSS_region_filtered2, '^', 'Color', 'r', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Type 2');
plot(AR_filtered3, high_WSS_region_filtered3, 's', 'Color', 'm', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Type 3');

% Generate x-values for fit line
AR_fit = linspace(min(AR_combined), max(AR_combined), 100);

% Fit and plot based on fit_type
if strcmp(fit_type, 'linear')
    % Linear fit
    p_linear = polyfit(AR_combined, high_WSS_combined, 1);
    yfit = polyval(p_linear, AR_fit);
    plot(AR_fit, yfit, '--k', 'LineWidth', 2, 'DisplayName', 'Linear Fit');

    % Calculate R-squared for linear fit
    yfit_combined = polyval(p_linear, AR_combined);
    SS_res = sum((high_WSS_combined - yfit_combined).^2);
    SS_tot = sum((high_WSS_combined - mean(high_WSS_combined)).^2);
    R_squared = 1 - (SS_res / SS_tot);

elseif strcmp(fit_type, 'power')
    % Power fit
    log_AR_combined = log(AR_combined);
    log_high_WSS_combined = log(high_WSS_combined);
    p_power = polyfit(log_AR_combined, log_high_WSS_combined, 1);
    a_power = exp(p_power(2));
    b_power = p_power(1);

    yfit = a_power * AR_fit.^b_power;
    plot(AR_fit, yfit, '--k', 'LineWidth', 2, 'DisplayName', 'Power Fit');

    % Calculate R-squared for power fit
    yfit_combined = a_power * AR_combined.^b_power;
    SS_res = sum((high_WSS_combined - yfit_combined).^2);
    SS_tot = sum((high_WSS_combined - mean(high_WSS_combined)).^2);
    R_squared = 1 - (SS_res / SS_tot);
else
    error('Invalid fit_type. Choose ''linear'' or ''power''.');
end

% Display R-squared value
text(mean(AR_combined), max(high_WSS_combined), ['R^2 = ', num2str(R_squared, 2)], 'FontSize', 12, 'Color', 'k');

xlabel('Aspect Ratio');
ylabel('High WSS Region');
title(['High WSS Region vs Aspect Ratio (' fit_type ' fit)']);
legend('Location', 'best');
grid off;
% Plot the Type4 points for High WSS Region
AR_filtered4_row1 = Type4(1, 1);  % Aspect Ratio for row 1
high_WSS_filtered4_row1 = Type4(1, 3);  % High WSS Region for row 1

AR_filtered4_row2 = Type4(2, 1);  % Aspect Ratio for row 2
high_WSS_filtered4_row2 = Type4(2, 3);  % High WSS Region for row 2

% Plot the first row point
plot(AR_filtered4_row1, high_WSS_filtered4_row1, 'h', 'Color', 'c', ...
     'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Type 4 (ANY_007)');

% Plot the second row point
plot(AR_filtered4_row2, high_WSS_filtered4_row2, '*', 'Color', 'k', ...
     'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Type 4 (ANY_111)');

legend('Location', 'best');
% --- Distance Calculation for Type4 Points ---
% Fit line parameters
m = p_linear(1); % Slope of the fit line
c = p_linear(2); % Intercept of the fit line

% Type4 Data Points
x4_row1 = AR_filtered4_row1; % Surface Area for ANY_007
y4_row1 = high_WSS_filtered4_row1; % Low WSS Region for ANY_007

x4_row2 = AR_filtered4_row2; % Surface Area for ANY_111
y4_row2 = high_WSS_filtered4_row2; % Low WSS Region for ANY_111

% Predicted y-values on the fit line for Type4 points
yfit_row1 = m * x4_row1 + c;
yfit_row2 = m * x4_row2 + c;

% Calculate the vertical distances (residuals)
distance_row1 = abs(y4_row1 - yfit_row1);
distance_row2 = abs(y4_row2 - yfit_row2);

% Display the distances
fprintf('high_WSS vs aspect ratio Distance of ANY_007 to the fit line: %.4f\n', distance_row1);
fprintf('high_WSS vs aspect ratio Distance of ANY_111 to the fit line: %.4f\n', distance_row2);

% Optionally annotate the plot with the distances
text(x4_row1, y4_row1, sprintf('%.4f', distance_row1), 'FontSize', 10, 'Color', 'c', 'VerticalAlignment', 'bottom');
text(x4_row2, y4_row2, sprintf('%.4f', distance_row2), 'FontSize', 10, 'Color', 'k', 'VerticalAlignment', 'bottom');
% --- High OSI Region vs Aspect Ratio ---
figure;

% Extract data for High OSI Region
[AR_filtered1, high_OSI_region_filtered1] = process_data(Type1, 1, 5);
[AR_filtered2, high_OSI_region_filtered2] = process_data(Type2, 1, 5);
[AR_filtered3, high_OSI_region_filtered3] = process_data(Type3, 1, 5);

% Combine data for fitting
AR_combined = [AR_filtered1; AR_filtered2; AR_filtered3];
high_OSI_combined = [high_OSI_region_filtered1; high_OSI_region_filtered2; high_OSI_region_filtered3];

% Plot data points
hold on;
plot(AR_filtered1, high_OSI_region_filtered1, 'o', 'Color', 'b', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Type 1');
plot(AR_filtered2, high_OSI_region_filtered2, '^', 'Color', 'r', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Type 2');
plot(AR_filtered3, high_OSI_region_filtered3, 's', 'Color', 'm', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Type 3');

% Apply linear fit
p = polyfit(AR_combined, high_OSI_combined, 1);
AR_fit = linspace(min(AR_combined), max(AR_combined), 100);
yfit = polyval(p, AR_fit);
plot(AR_fit, yfit, '--k', 'LineWidth', 2, 'DisplayName', 'Linear Fit');

% Calculate R-squared
SS_res = sum((high_OSI_combined - polyval(p, AR_combined)).^2);
SS_tot = sum((high_OSI_combined - mean(high_OSI_combined)).^2);
R_squared = 1 - (SS_res / SS_tot);
%text(mean(AR_combined), max(high_OSI_combined), ['R^2 = ', num2str(R_squared, 2)], 'FontSize', 12, 'Color', 'k');

xlabel('Aspect Ratio');
ylabel('High OSI Region');
title('High OSI Region vs Aspect Ratio');
legend('Location', 'best');
grid off;
% Plot the Type4 points for High OSI Region
AR_filtered4_row1 = Type4(1, 1);  % Aspect Ratio for row 1
high_OSI_filtered4_row1 = Type4(1, 5);  % High OSI Region for row 1

AR_filtered4_row2 = Type4(2, 1);  % Aspect Ratio for row 2
high_OSI_filtered4_row2 = Type4(2, 5);  % High OSI Region for row 2

% Plot the first row point
plot(AR_filtered4_row1, high_OSI_filtered4_row1, 'h', 'Color', 'c', ...
     'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Type 4 (ANY_007)');

% Plot the second row point
plot(AR_filtered4_row2, high_OSI_filtered4_row2, '*', 'Color', 'k', ...
     'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Type 4 (ANY_111)');

legend('Location', 'best');
% --- Distance Calculation for Type4 Points ---
% Fit line parameters
m = p(1); % Slope of the fit line
c = p(2); % Intercept of the fit line

% Type4 Data Points
x4_row1 = AR_filtered4_row1; % Surface Area for ANY_007
y4_row1 = high_OSI_filtered4_row1; % Low WSS Region for ANY_007

x4_row2 = AR_filtered4_row2; % Surface Area for ANY_111
y4_row2 = high_OSI_filtered4_row2; % Low WSS Region for ANY_111

% Predicted y-values on the fit line for Type4 points
yfit_row1 = m * x4_row1 + c;
yfit_row2 = m * x4_row2 + c;

% Calculate the vertical distances (residuals)
distance_row1 = abs(y4_row1 - yfit_row1);
distance_row2 = abs(y4_row2 - yfit_row2);

% Display the distances
fprintf('high_OSI vs aspect ratio Distance of ANY_007 to the fit line: %.4f\n', distance_row1);
fprintf('high_OSI vs aspect ratio Distance of ANY_111 to the fit line: %.4f\n', distance_row2);

% Optionally annotate the plot with the distances
text(x4_row1, y4_row1, sprintf('%.4f', distance_row1), 'FontSize', 10, 'Color', 'c', 'VerticalAlignment', 'bottom');
text(x4_row2, y4_row2, sprintf('%.4f', distance_row2), 'FontSize', 10, 'Color', 'k', 'VerticalAlignment', 'bottom');
% --- High RRT Region vs Aspect Ratio ---
figure;

% Extract data for High RRT Region
[AR_filtered1, high_RRT_region_filtered1] = process_data(Type1, 1, 7);
[AR_filtered2, high_RRT_region_filtered2] = process_data(Type2, 1, 7);
[AR_filtered3, high_RRT_region_filtered3] = process_data(Type3, 1, 7);

% Combine data for fitting
AR_combined = [AR_filtered1; AR_filtered2; AR_filtered3];
high_RRT_combined = [high_RRT_region_filtered1; high_RRT_region_filtered2; high_RRT_region_filtered3];

% Plot data points
hold on;
plot(AR_filtered1, high_RRT_region_filtered1, 'o', 'Color', 'b', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Type 1');
plot(AR_filtered2, high_RRT_region_filtered2, '^', 'Color', 'r', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Type 2');
plot(AR_filtered3, high_RRT_region_filtered3, 's', 'Color', 'm', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Type 3');

% Apply linear fit
p = polyfit(AR_combined, high_RRT_combined, 1);
AR_fit = linspace(min(AR_combined), max(AR_combined), 100);
yfit = polyval(p, AR_fit);
plot(AR_fit, yfit, '--k', 'LineWidth', 2, 'DisplayName', 'Linear Fit');

% Calculate R-squared
SS_res = sum((high_RRT_combined - polyval(p, AR_combined)).^2);
SS_tot = sum((high_RRT_combined - mean(high_RRT_combined)).^2);
R_squared = 1 - (SS_res / SS_tot);
text(mean(AR_combined), max(high_RRT_combined), ['R^2 = ', num2str(R_squared, 2)], 'FontSize', 12, 'Color', 'k');

xlabel('Aspect Ratio');
ylabel('High RRT Region');
title('High RRT Region vs Aspect Ratio');
legend('Location', 'best');
grid off;
% Plot the Type4 points for High RRT Region
AR_filtered4_row1 = Type4(1, 1);  % Aspect Ratio for row 1
high_RRT_filtered4_row1 = Type4(1, 7);  % High RRT Region for row 1

AR_filtered4_row2 = Type4(2, 1);  % Aspect Ratio for row 2
high_RRT_filtered4_row2 = Type4(2, 7);  % High RRT Region for row 2

% Plot the first row point
plot(AR_filtered4_row1, high_RRT_filtered4_row1, 'h', 'Color', 'c', ...
     'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Type 4 (ANY_007)');

% Plot the second row point
plot(AR_filtered4_row2, high_RRT_filtered4_row2, '*', 'Color', 'k', ...
     'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Type 4 (ANY_111)');

legend('Location', 'best');
% --- Distance Calculation for Type4 Points ---
% Fit line parameters
m = p(1); % Slope of the fit line
c = p(2); % Intercept of the fit line

% Type4 Data Points
x4_row1 = AR_filtered4_row1; % Surface Area for ANY_007
y4_row1 = high_RRT_filtered4_row1; % Low WSS Region for ANY_007

x4_row2 = AR_filtered4_row2; % Surface Area for ANY_111
y4_row2 = high_RRT_filtered4_row2; % Low WSS Region for ANY_111

% Predicted y-values on the fit line for Type4 points
yfit_row1 = m * x4_row1 + c;
yfit_row2 = m * x4_row2 + c;

% Calculate the vertical distances (residuals)
distance_row1 = abs(y4_row1 - yfit_row1);
distance_row2 = abs(y4_row2 - yfit_row2);

% Display the distances
fprintf('high_RRT vs aspect ratio Distance of ANY_007 to the fit line: %.4f\n', distance_row1);
fprintf('high_RRT vs aspect ratio Distance of ANY_111 to the fit line: %.4f\n', distance_row2);

% Optionally annotate the plot with the distances
text(x4_row1, y4_row1, sprintf('%.4f', distance_row1), 'FontSize', 10, 'Color', 'c', 'VerticalAlignment', 'bottom');
text(x4_row2, y4_row2, sprintf('%.4f', distance_row2), 'FontSize', 10, 'Color', 'k', 'VerticalAlignment', 'bottom');
% --- Vorticity Plot ---
figure;

% Extract, filter, and sort relevant columns for Volume Normalized Vorticity
[AR_filtered1, vorticity_filtered1] = process_data(Type1, 1, 10);
[AR_filtered2, vorticity_filtered2] = process_data(Type2, 1, 10);
[AR_filtered3, vorticity_filtered3] = process_data(Type3, 1, 10);

% Combine data for fitting
AR_combined = [AR_filtered1; AR_filtered2; AR_filtered3];
vorticity_combined = [vorticity_filtered1; vorticity_filtered2; vorticity_filtered3];

% Plot for Type1, Type2, and Type3 for Vorticity as dots
hold on;
plot(AR_filtered1, vorticity_filtered1, 'o', 'Color', 'b', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Type 1');
plot(AR_filtered2, vorticity_filtered2, '^', 'Color', 'r', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Type 2');
plot(AR_filtered3, vorticity_filtered3, 's', 'Color', 'm', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Type 3');

% Apply linear fit for Vorticity
p = polyfit(AR_combined, vorticity_combined, 1);
AR_fit = linspace(min(AR_combined), max(AR_combined), 100);
yfit = polyval(p, AR_fit);
plot(AR_fit, yfit, '--k', 'LineWidth', 2, 'DisplayName', 'Linear Fit');

% Calculate R-squared
SS_res = sum((vorticity_combined - polyval(p, AR_combined)).^2);
SS_tot = sum((vorticity_combined - mean(vorticity_combined)).^2);
R_squared = 1 - (SS_res / SS_tot);
text(mean(AR_combined), max(vorticity_combined), ['R^2 = ', num2str(R_squared, 2)], 'FontSize', 12, 'Color', 'k');

xlabel('Aspect Ratio');
ylabel('Volume Normalized Vorticity');
title('Volume Normalized Vorticity vs Aspect Ratio');
legend('Location', 'best');
ylim([100, 300]);

grid off;
% Plot the Type4 points for Vorticity
AR_filtered4_row1 = Type4(1, 1);  % Aspect Ratio for row 1
vorticity_filtered4_row1 = Type4(1, 10);  % Vorticity for row 1

AR_filtered4_row2 = Type4(2, 1);  % Aspect Ratio for row 2
vorticity_filtered4_row2 = Type4(2, 10);  % Vorticity for row 2

% Plot the first row point
plot(AR_filtered4_row1, vorticity_filtered4_row1, 'h', 'Color', 'c', ...
     'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Type 4 (ANY_007)');

% Plot the second row point
plot(AR_filtered4_row2, vorticity_filtered4_row2, '*', 'Color', 'k', ...
     'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Type 4 (ANY_111)');

% Update the legend
legend('Location', 'best');
% --- Distance Calculation for Type4 Points ---
% Fit line parameters
m = p(1); % Slope of the fit line
c = p(2); % Intercept of the fit line

% Type4 Data Points
x4_row1 = AR_filtered4_row1; % Surface Area for ANY_007
y4_row1 = vorticity_filtered4_row1; % Low WSS Region for ANY_007

x4_row2 = AR_filtered4_row2; % Surface Area for ANY_111
y4_row2 = vorticity_filtered4_row2; % Low WSS Region for ANY_111

% Predicted y-values on the fit line for Type4 points
yfit_row1 = m * x4_row1 + c;
yfit_row2 = m * x4_row2 + c;

% Calculate the vertical distances (residuals)
distance_row1 = abs(y4_row1 - yfit_row1);
distance_row2 = abs(y4_row2 - yfit_row2);

% Display the distances
fprintf('vorticity vs aspect ratio Distance of ANY_007 to the fit line: %.4f\n', distance_row1);
fprintf('vorticity vs aspect ratio Distance of ANY_111 to the fit line: %.4f\n', distance_row2);

% Optionally annotate the plot with the distances
text(x4_row1, y4_row1, sprintf('%.4f', distance_row1), 'FontSize', 10, 'Color', 'c', 'VerticalAlignment', 'bottom');
text(x4_row2, y4_row2, sprintf('%.4f', distance_row2), 'FontSize', 10, 'Color', 'k', 'VerticalAlignment', 'bottom');
% --- qcrit Vortex Plot ---
figure;

% Extract, filter, and sort relevant columns for qcrit vortex
[AR_filtered1, vortex_filtered1] = process_data(Type1, 1, 13);
[AR_filtered2, vortex_filtered2] = process_data(Type2, 1, 13);
[AR_filtered3, vortex_filtered3] = process_data(Type3, 1, 13);

% Combine data for fitting
AR_combined = [AR_filtered1; AR_filtered2; AR_filtered3];
vortex_combined = [vortex_filtered1; vortex_filtered2; vortex_filtered3];

% Plot for Type1, Type2, and Type3 for qcrit vortex as dots
hold on;
plot(AR_filtered1, vortex_filtered1, 'o', 'Color', 'b', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Type 1');
plot(AR_filtered2, vortex_filtered2, '^', 'Color', 'r', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Type 2');
plot(AR_filtered3, vortex_filtered3, 's', 'Color', 'm', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Type 3');

% Apply power fit for qcrit vortex
log_AR_combined = log(AR_combined);
log_vortex_combined = log(vortex_combined);
p = polyfit(log_AR_combined, log_vortex_combined, 1);
a = exp(p(2));
b = p(1);

% Generate power fit line
AR_fit = linspace(min(AR_combined), max(AR_combined), 100);
yfit = a * AR_fit.^b;
plot(AR_fit, yfit, '--k', 'LineWidth', 2, 'DisplayName', 'Power Fit');

% Calculate R-squared for power fit
vortex_fit_combined = a * AR_combined.^b;
SS_res = sum((vortex_combined - vortex_fit_combined).^2);
SS_tot = sum((vortex_combined - mean(vortex_combined)).^2);
R_squared = 1 - (SS_res / SS_tot);
text(mean(AR_combined), max(vortex_combined), ['R^2 = ', num2str(R_squared, 2)], 'FontSize', 12, 'Color', 'k');

xlabel('Aspect Ratio');
ylabel('qcrit Vortex');
title('qcrit Vortex vs Aspect Ratio');
legend('Location', 'best');
ylim([2, 18]);

grid off;
% Plot the Type4 points for Vortex (qcrit)
AR_filtered4_row1 = Type4(1, 1);  % Aspect Ratio for row 1
vortex_filtered4_row1 = Type4(1, 13);  % Vortex (qcrit) for row 1

AR_filtered4_row2 = Type4(2, 1);  % Aspect Ratio for row 2
vortex_filtered4_row2 = Type4(2, 13);  % Vortex (qcrit) for row 2

% Plot the first row point
plot(AR_filtered4_row1, vortex_filtered4_row1, 'h', 'Color', 'c', ...
     'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Type 4 (ANY_007)');

% Plot the second row point
plot(AR_filtered4_row2, vortex_filtered4_row2, '*', 'Color', 'k', ...
     'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Type 4 (ANY_111)');

% Update the legend
legend('Location', 'best');
% --- Distance Calculation for Type4 Points ---

% Power fit parameters
a = exp(p(2)); % Scale factor (intercept in log-log fit)
b = p(1);      % Exponent (slope in log-log fit)

% Type 4 Data Points
x4_row1 = AR_filtered4_row1; % Independent variable (e.g., Surface Area for ANY_007)
y4_row1 = vortex_filtered4_row1; % Dependent variable (e.g., Low WSS Region for ANY_007)

x4_row2 = AR_filtered4_row2; % Independent variable (e.g., Surface Area for ANY_111)
y4_row2 = vortex_filtered4_row2; % Dependent variable (e.g., Low WSS Region for ANY_111)

% Predicted y-values on the power fit line for Type 4 points
yfit_row1 = a * (x4_row1^b);
yfit_row2 = a * (x4_row2^b);

% Calculate the vertical distances (residuals)
distance_row1 = abs(y4_row1 - yfit_row1);
distance_row2 = abs(y4_row2 - yfit_row2);

% Display the distances
fprintf('vortex vs area ratio Distance of ANY_007 to the power fit line: %.4f\n', distance_row1);
fprintf('vortex vs area ratio Distance of ANY_111 to the power fit line: %.4f\n', distance_row2);

% Optionally annotate the plot with the distances
text(x4_row1, y4_row1, sprintf('%.4f', distance_row1), 'FontSize', 10, 'Color', 'c', 'VerticalAlignment', 'bottom');
text(x4_row2, y4_row2, sprintf('%.4f', distance_row2), 'FontSize', 10, 'Color', 'k', 'VerticalAlignment', 'bottom');
% --- Helper function for filtering and sorting data ---
function [x_filtered, y_filtered] = process_data(Type, x_col, y_col)
    x = Type(:, x_col);
    y = Type(:, y_col);
    
    % Filter out zero or NaN values
    valid_indices = ~isnan(y) & y ~= 0;
    x_filtered = x(valid_indices);
    y_filtered = y(valid_indices);
    
    % Sort by x-axis variable
    [x_filtered, sort_idx] = sort(x_filtered);
    y_filtered = y_filtered(sort_idx);
end
