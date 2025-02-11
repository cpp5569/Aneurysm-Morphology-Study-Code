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

% --- Low WSS Region vs Area Ratio ---
figure;

% Extract, filter, and sort relevant columns for Low WSS with Area Ratio
[area_ratio_filtered1, low_WSS_region_filtered1] = process_data(Type1, 9, 2);
[area_ratio_filtered2, low_WSS_region_filtered2] = process_data(Type2, 9, 2);
[area_ratio_filtered3, low_WSS_region_filtered3] = process_data(Type3, 9, 2);

% Combine data for fitting
area_combined = [area_ratio_filtered1; area_ratio_filtered2; area_ratio_filtered3];
low_WSS_combined = [low_WSS_region_filtered1; low_WSS_region_filtered2; low_WSS_region_filtered3];

% Plot data points
hold on;
plot(area_ratio_filtered1, low_WSS_region_filtered1, 'o', 'Color', 'b', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Type 1');
plot(area_ratio_filtered2, low_WSS_region_filtered2, '^', 'Color', 'r', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Type 2');
plot(area_ratio_filtered3, low_WSS_region_filtered3, 's', 'Color', 'm', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Type 3');

% Generate x-values for fit line
area_fit = linspace(min(area_combined), max(area_combined), 100);

% Fit and plot based on fit_type
p = polyfit(area_combined, low_WSS_combined, 1);
yfit = polyval(p, area_fit);
plot(area_fit, yfit, '--k', 'LineWidth', 2, 'DisplayName', 'Linear Fit');

% Calculate R-squared for linear fit
yfit_combined = polyval(p, area_combined);
SS_res = sum((low_WSS_combined - yfit_combined).^2);
SS_tot = sum((low_WSS_combined - mean(low_WSS_combined)).^2);
R_squared = 1 - (SS_res / SS_tot);

% Display R-squared value
text(mean(area_combined), max(low_WSS_combined), ['R^2 = ', num2str(R_squared, 2)], 'FontSize', 12, 'Color', 'k');

xlabel('Area Ratio');
ylabel('Low WSS Region');
title('Low WSS Region vs Area Ratio');
legend('Location', 'best');
grid off;

% Plot the Type4 points for Low WSS Region vs Area Ratio
area_ratio_filtered4_row1 = Type4(1, 9);  % Area Ratio for row 1
low_WSS_filtered4_row1 = Type4(1, 2);    % Low WSS Region for row 1

area_ratio_filtered4_row2 = Type4(2, 9);  % Area Ratio for row 2
low_WSS_filtered4_row2 = Type4(2, 2);    % Low WSS Region for row 2

% Plot the first row point
plot(area_ratio_filtered4_row1, low_WSS_filtered4_row1, 'h', 'Color', 'c', ...
     'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Type 4 (ANY_007)');

% Plot the second row point
plot(area_ratio_filtered4_row2, low_WSS_filtered4_row2, '*', 'Color', 'k', ...
     'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Type 4 (ANY_111)');

% Update the legend
legend('Location', 'best');

% --- Distance Calculation for Type4 Points ---
% Fit line parameters
m = p(1); % Slope of the fit line
c = p(2); % Intercept of the fit line

% Type4 Data Points
x4_row1 = area_ratio_filtered4_row1; % Surface Area for ANY_007
y4_row1 = low_WSS_filtered4_row1; % Low WSS Region for ANY_007

x4_row2 = area_ratio_filtered4_row2; % Surface Area for ANY_111
y4_row2 = low_WSS_filtered4_row2; % Low WSS Region for ANY_111

% Predicted y-values on the fit line for Type4 points
yfit_row1 = m * x4_row1 + c;
yfit_row2 = m * x4_row2 + c;

% Calculate the vertical distances (residuals)
distance_row1 = abs(y4_row1 - yfit_row1);
distance_row2 = abs(y4_row2 - yfit_row2);

% Display the distances
fprintf('low_WSS vs area ratio Distance of ANY_007 to the fit line: %.4f\n', distance_row1);
fprintf('low_WSS vs area ratio Distance of ANY_111 to the fit line: %.4f\n', distance_row2);

% Optionally annotate the plot with the distances
text(x4_row1, y4_row1, sprintf('%.4f', distance_row1), 'FontSize', 10, 'Color', 'c', 'VerticalAlignment', 'bottom');
text(x4_row2, y4_row2, sprintf('%.4f', distance_row2), 'FontSize', 10, 'Color', 'k', 'VerticalAlignment', 'bottom');
% --- High WSS Region vs Area Ratio ---
figure;

% Extract, filter, and sort relevant columns for High WSS with Area Ratio
[area_ratio_filtered1, high_WSS_region_filtered1] = process_data(Type1, 9, 3);
[area_ratio_filtered2, high_WSS_region_filtered2] = process_data(Type2, 9, 3);
[area_ratio_filtered3, high_WSS_region_filtered3] = process_data(Type3, 9, 3);

% Combine data for fitting
area_combined = [area_ratio_filtered1; area_ratio_filtered2; area_ratio_filtered3];
high_WSS_combined = [high_WSS_region_filtered1; high_WSS_region_filtered2; high_WSS_region_filtered3];

% Plot data points
hold on;
plot(area_ratio_filtered1, high_WSS_region_filtered1, 'o', 'Color', 'b', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Type 1');
plot(area_ratio_filtered2, high_WSS_region_filtered2, '^', 'Color', 'r', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Type 2');
plot(area_ratio_filtered3, high_WSS_region_filtered3, 's', 'Color', 'm', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Type 3');

% Generate x-values for fit line
area_fit = linspace(min(area_combined), max(area_combined), 100);

% Linear fit across all points
p = polyfit(area_combined, high_WSS_combined, 1);
yfit = polyval(p, area_fit);
plot(area_fit, yfit, '--k', 'LineWidth', 2, 'DisplayName', 'Linear Fit');

% Calculate R-squared
SS_res = sum((high_WSS_combined - polyval(p, area_combined)).^2);
SS_tot = sum((high_WSS_combined - mean(high_WSS_combined)).^2);
R_squared = 1 - (SS_res / SS_tot);
text(mean(area_combined), max(high_WSS_combined), ['R^2 = ', num2str(R_squared, 2)], 'FontSize', 12, 'Color', 'k');

xlabel('Area Ratio');
ylabel('High WSS Region');
title('High WSS Region vs Area Ratio');
legend('Location', 'best');
grid off;

% Plot the Type4 points for High WSS Region vs Area Ratio
area_ratio_filtered4_row1 = Type4(1, 9);  % Area Ratio for row 1
high_WSS_filtered4_row1 = Type4(1, 3);   % High WSS Region for row 1

area_ratio_filtered4_row2 = Type4(2, 9);  % Area Ratio for row 2
high_WSS_filtered4_row2 = Type4(2, 3);   % High WSS Region for row 2

% Plot the first row point
plot(area_ratio_filtered4_row1, high_WSS_filtered4_row1, 'h', 'Color', 'c', ...
     'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Type 4 (ANY_007)');

% Plot the second row point
plot(area_ratio_filtered4_row2, high_WSS_filtered4_row2, '*', 'Color', 'k', ...
     'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Type 4 (ANY_111)');

% Update the legend
legend('Location', 'best');
% --- Distance Calculation for Type4 Points ---
% Fit line parameters
m = p(1); % Slope of the fit line
c = p(2); % Intercept of the fit line

% Type4 Data Points
x4_row1 = area_ratio_filtered4_row1; % Surface Area for ANY_007
y4_row1 = high_WSS_filtered4_row1; % Low WSS Region for ANY_007

x4_row2 = area_ratio_filtered4_row2; % Surface Area for ANY_111
y4_row2 = high_WSS_filtered4_row2; % Low WSS Region for ANY_111

% Predicted y-values on the fit line for Type4 points
yfit_row1 = m * x4_row1 + c;
yfit_row2 = m * x4_row2 + c;

% Calculate the vertical distances (residuals)
distance_row1 = abs(y4_row1 - yfit_row1);
distance_row2 = abs(y4_row2 - yfit_row2);

% Display the distances
fprintf('high_WSS vs area ratio Distance of ANY_007 to the fit line: %.4f\n', distance_row1);
fprintf('high_WSS vs area ratio Distance of ANY_111 to the fit line: %.4f\n', distance_row2);

% Optionally annotate the plot with the distances
text(x4_row1, y4_row1, sprintf('%.4f', distance_row1), 'FontSize', 10, 'Color', 'c', 'VerticalAlignment', 'bottom');
text(x4_row2, y4_row2, sprintf('%.4f', distance_row2), 'FontSize', 10, 'Color', 'k', 'VerticalAlignment', 'bottom');
% --- Vortex (qcrit) vs Area Ratio ---
figure;

% Extract, filter, and sort relevant columns for Vortex with Area Ratio
[area_ratio_filtered1, vortex_filtered1] = process_data(Type1, 9, 13);
[area_ratio_filtered2, vortex_filtered2] = process_data(Type2, 9, 13);
[area_ratio_filtered3, vortex_filtered3] = process_data(Type3, 9, 13);

% Combine data for fitting
area_combined = [area_ratio_filtered1; area_ratio_filtered2; area_ratio_filtered3];
vortex_combined = [vortex_filtered1; vortex_filtered2; vortex_filtered3];

% Plot data points
hold on;
plot(area_ratio_filtered1, vortex_filtered1, 'o', 'Color', 'b', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Type 1');
plot(area_ratio_filtered2, vortex_filtered2, '^', 'Color', 'r', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Type 2');
plot(area_ratio_filtered3, vortex_filtered3, 's', 'Color', 'm', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Type 3');

% Generate x-values for fit line
area_fit = linspace(min(area_combined), max(area_combined), 100);

% Power fit (y = a * x^b)
log_area_combined = log(area_combined);
log_vortex_combined = log(vortex_combined);
p_power = polyfit(log_area_combined, log_vortex_combined, 1);
a_power = exp(p_power(2)); % Convert intercept back from log scale
b_power = p_power(1);      % Slope is the exponent

% Generate power fit line
yfit = a_power * area_fit.^b_power;
plot(area_fit, yfit, '--k', 'LineWidth', 2, 'DisplayName', 'Power Fit');

% Calculate R-squared for power fit
vortex_fit = a_power * area_combined.^b_power;
SS_res = sum((vortex_combined - vortex_fit).^2);
SS_tot = sum((vortex_combined - mean(vortex_combined)).^2);
R_squared = 1 - (SS_res / SS_tot);
text(mean(area_combined), max(vortex_combined), ['R^2 = ', num2str(R_squared, 2)], 'FontSize', 12, 'Color', 'k');

xlabel('Area Ratio');
ylabel('qcrit Vortex');
title('qcrit Vortex vs Area Ratio');
legend('Location', 'best');
grid off;

% Plot the Type4 points for Vortex vs Area Ratio
area_ratio_filtered4_row1 = Type4(1, 9);  % Area Ratio for row 1
vortex_filtered4_row1 = Type4(1, 13);    % Vortex for row 1

area_ratio_filtered4_row2 = Type4(2, 9);  % Area Ratio for row 2
vortex_filtered4_row2 = Type4(2, 13);    % Vortex for row 2

% Plot the first row point
plot(area_ratio_filtered4_row1, vortex_filtered4_row1, 'h', 'Color', 'c', ...
     'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Type 4 (ANY_007)');

% Plot the second row point
plot(area_ratio_filtered4_row2, vortex_filtered4_row2, '*', 'Color', 'k', ...
     'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Type 4 (ANY_111)');

% Update the legend
legend('Location', 'best');
% --- Distance Calculation for Type4 Points ---

% Power fit parameters
a = exp(p_power(2)); % Scale factor (intercept in log-log fit)
b = p_power(1);      % Exponent (slope in log-log fit)

% Type 4 Data Points
x4_row1 = area_ratio_filtered4_row1; % Independent variable (e.g., Surface Area for ANY_007)
y4_row1 = vortex_filtered4_row1; % Dependent variable (e.g., Low WSS Region for ANY_007)

x4_row2 = area_ratio_filtered4_row2; % Independent variable (e.g., Surface Area for ANY_111)
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

% --- High OSI Region vs Area Ratio ---
figure;

% Extract, filter, and sort relevant columns for High OSI with Area Ratio
[area_ratio_filtered1, high_OSI_region_filtered1] = process_data(Type1, 9, 5);
[area_ratio_filtered2, high_OSI_region_filtered2] = process_data(Type2, 9, 5);
[area_ratio_filtered3, high_OSI_region_filtered3] = process_data(Type3, 9, 5);

% Combine data for fitting
area_combined = [area_ratio_filtered1; area_ratio_filtered2; area_ratio_filtered3];
high_OSI_combined = [high_OSI_region_filtered1; high_OSI_region_filtered2; high_OSI_region_filtered3];

% Plot data points
hold on;
plot(area_ratio_filtered1, high_OSI_region_filtered1, 'o', 'Color', 'b', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Type 1');
plot(area_ratio_filtered2, high_OSI_region_filtered2, '^', 'Color', 'r', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Type 2');
plot(area_ratio_filtered3, high_OSI_region_filtered3, 's', 'Color', 'm', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Type 3');

% Generate x-values for fit line
area_fit = linspace(min(area_combined), max(area_combined), 100);

% Linear fit across all points
p = polyfit(area_combined, high_OSI_combined, 1);
yfit = polyval(p, area_fit);
plot(area_fit, yfit, '--k', 'LineWidth', 2, 'DisplayName', 'Linear Fit');

% Calculate R-squared
SS_res = sum((high_OSI_combined - polyval(p, area_combined)).^2);
SS_tot = sum((high_OSI_combined - mean(high_OSI_combined)).^2);
R_squared = 1 - (SS_res / SS_tot);
%text(mean(area_combined), max(high_OSI_combined), ['R^2 = ', num2str(R_squared, 2)], 'FontSize', 12, 'Color', 'k');

xlabel('Area Ratio');
ylabel('High OSI Region');
title('High OSI Region vs Area Ratio');
legend('Location', 'best');
grid off;

% Plot the Type4 points for High OSI Region vs Area Ratio
area_ratio_filtered4_row1 = Type4(1, 9);  % Area Ratio for row 1
high_OSI_filtered4_row1 = Type4(1, 5);   % High OSI Region for row 1

area_ratio_filtered4_row2 = Type4(2, 9);  % Area Ratio for row 2
high_OSI_filtered4_row2 = Type4(2, 5);   % High OSI Region for row 2

% Plot the first row point
plot(area_ratio_filtered4_row1, high_OSI_filtered4_row1, 'h', 'Color', 'c', ...
     'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Type 4 (ANY_007)');

% Plot the second row point
plot(area_ratio_filtered4_row2, high_OSI_filtered4_row2, '*', 'Color', 'k', ...
     'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Type 4 (ANY_111)');

% Update the legend
legend('Location', 'best');
% --- Distance Calculation for Type4 Points ---
% Fit line parameters
m = p(1); % Slope of the fit line
c = p(2); % Intercept of the fit line

% Type4 Data Points
x4_row1 = area_ratio_filtered4_row1; % Surface Area for ANY_007
y4_row1 = high_OSI_filtered4_row1; % Low WSS Region for ANY_007

x4_row2 = area_ratio_filtered4_row2; % Surface Area for ANY_111
y4_row2 = high_OSI_filtered4_row2; % Low WSS Region for ANY_111

% Predicted y-values on the fit line for Type4 points
yfit_row1 = m * x4_row1 + c;
yfit_row2 = m * x4_row2 + c;

% Calculate the vertical distances (residuals)
distance_row1 = abs(y4_row1 - yfit_row1);
distance_row2 = abs(y4_row2 - yfit_row2);

% Display the distances
fprintf('OSI vs area ratio Distance of ANY_007 to the fit line: %.4f\n', distance_row1);
fprintf('OSI vs area ratio Distance of ANY_111 to the fit line: %.4f\n', distance_row2);

% Optionally annotate the plot with the distances
text(x4_row1, y4_row1, sprintf('%.4f', distance_row1), 'FontSize', 10, 'Color', 'c', 'VerticalAlignment', 'bottom');
text(x4_row2, y4_row2, sprintf('%.4f', distance_row2), 'FontSize', 10, 'Color', 'k', 'VerticalAlignment', 'bottom');
% --- High RRT Region vs Area Ratio ---
figure;

% Extract, filter, and sort relevant columns for High RRT with Area Ratio
[area_ratio_filtered1, high_RRT_region_filtered1] = process_data(Type1, 9, 7);
[area_ratio_filtered2, high_RRT_region_filtered2] = process_data(Type2, 9, 7);
[area_ratio_filtered3, high_RRT_region_filtered3] = process_data(Type3, 9, 7);

% Combine data for fitting
area_combined = [area_ratio_filtered1; area_ratio_filtered2; area_ratio_filtered3];
high_RRT_combined = [high_RRT_region_filtered1; high_RRT_region_filtered2; high_RRT_region_filtered3];

% Plot data points
hold on;
plot(area_ratio_filtered1, high_RRT_region_filtered1, 'o', 'Color', 'b', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Type 1');
plot(area_ratio_filtered2, high_RRT_region_filtered2, '^', 'Color', 'r', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Type 2');
plot(area_ratio_filtered3, high_RRT_region_filtered3, 's', 'Color', 'm', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Type 3');

% Generate x-values for fit line
area_fit = linspace(min(area_combined), max(area_combined), 100);

% Linear fit across all points
p = polyfit(area_combined, high_RRT_combined, 1);
yfit = polyval(p, area_fit);
plot(area_fit, yfit, '--k', 'LineWidth', 2, 'DisplayName', 'Linear Fit');

% Calculate R-squared
SS_res = sum((high_RRT_combined - polyval(p, area_combined)).^2);
SS_tot = sum((high_RRT_combined - mean(high_RRT_combined)).^2);
R_squared = 1 - (SS_res / SS_tot);
text(mean(area_combined), max(high_RRT_combined), ['R^2 = ', num2str(R_squared, 2)], 'FontSize', 12, 'Color', 'k');

xlabel('Area Ratio');
ylabel('High RRT Region');
title('High RRT Region vs Area Ratio');
legend('Location', 'best');
grid off;

% Plot the Type4 points for High RRT Region vs Area Ratio
area_ratio_filtered4_row1 = Type4(1, 9);  % Area Ratio for row 1
high_RRT_filtered4_row1 = Type4(1, 7);   % High RRT Region for row 1

area_ratio_filtered4_row2 = Type4(2, 9);  % Area Ratio for row 2
high_RRT_filtered4_row2 = Type4(2, 7);   % High RRT Region for row 2

% Plot the first row point
plot(area_ratio_filtered4_row1, high_RRT_filtered4_row1, 'h', 'Color', 'c', ...
     'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Type 4 (ANY_007)');

% Plot the second row point
plot(area_ratio_filtered4_row2, high_RRT_filtered4_row2, '*', 'Color', 'k', ...
     'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Type 4 (ANY_111)');

% Update the legend
legend('Location', 'best');
% --- Distance Calculation for Type4 Points ---
% Fit line parameters
m = p(1); % Slope of the fit line
c = p(2); % Intercept of the fit line

% Type4 Data Points
x4_row1 = area_ratio_filtered4_row1; % Surface Area for ANY_007
y4_row1 = high_RRT_filtered4_row1; % Low WSS Region for ANY_007

x4_row2 = area_ratio_filtered4_row2; % Surface Area for ANY_111
y4_row2 = high_RRT_filtered4_row2; % Low WSS Region for ANY_111

% Predicted y-values on the fit line for Type4 points
yfit_row1 = m * x4_row1 + c;
yfit_row2 = m * x4_row2 + c;

% Calculate the vertical distances (residuals)
distance_row1 = abs(y4_row1 - yfit_row1);
distance_row2 = abs(y4_row2 - yfit_row2);

% Display the distances
fprintf('RRT vs area ratio Distance of ANY_007 to the fit line: %.4f\n', distance_row1);
fprintf('RRT vs area ratio Distance of ANY_111 to the fit line: %.4f\n', distance_row2);

% Optionally annotate the plot with the distances
text(x4_row1, y4_row1, sprintf('%.4f', distance_row1), 'FontSize', 10, 'Color', 'c', 'VerticalAlignment', 'bottom');
text(x4_row2, y4_row2, sprintf('%.4f', distance_row2), 'FontSize', 10, 'Color', 'k', 'VerticalAlignment', 'bottom');
% --- Vorticity vs Area Ratio ---
figure;

% Extract, filter, and sort relevant columns for Vorticity with Area Ratio
[area_ratio_filtered1, vorticity_filtered1] = process_data(Type1, 9, 10); % Area Ratio, Volume Normalized Vorticity
[area_ratio_filtered2, vorticity_filtered2] = process_data(Type2, 9, 10);
[area_ratio_filtered3, vorticity_filtered3] = process_data(Type3, 9, 10);

% Combine the data
area_combined = [area_ratio_filtered1; area_ratio_filtered2; area_ratio_filtered3];
vorticity_combined = [vorticity_filtered1; vorticity_filtered2; vorticity_filtered3];

% Sort the combined data
[area_combined, sort_idx] = sort(area_combined);
vorticity_combined = vorticity_combined(sort_idx);

% Plot dots for Type 1, Type 2, and Type 3
hold on;
plot(area_ratio_filtered1, vorticity_filtered1, 'o', 'Color', 'b', 'MarkerSize', 8, 'LineWidth', 2);
plot(area_ratio_filtered2, vorticity_filtered2, '^', 'Color', 'r', 'MarkerSize', 8, 'LineWidth', 2);
plot(area_ratio_filtered3, vorticity_filtered3, 's', 'Color', 'm', 'MarkerSize', 8, 'LineWidth', 2);

% Linear fit across all points
p = polyfit(area_combined, vorticity_combined, 1);
yfit = polyval(p, area_combined);
plot(area_combined, yfit, '--k', 'LineWidth', 2);

% Calculate R-squared
SS_res = sum((vorticity_combined - yfit).^2);
SS_tot = sum((vorticity_combined - mean(vorticity_combined)).^2);
R_squared = 1 - (SS_res / SS_tot);
text(mean(area_combined), max(vorticity_combined), ['R^2 = ', num2str(R_squared, 2)], 'FontSize', 12, 'Color', 'k');

xlabel('Area Ratio');
ylabel('Volume Normalized Vorticity');
title('Vorticity vs Area Ratio');
ylim([100, 300]);

%legend({'Type 1', 'Type 2', 'Type 3', 'Linear Fit'});
grid off;
% Plot the Type4 points for Vorticity vs Area Ratio
area_ratio_filtered4_row1 = Type4(1, 9);  % Area Ratio for row 1
vorticity_filtered4_row1 = Type4(1, 10);  % Vorticity for row 1

area_ratio_filtered4_row2 = Type4(2, 9);  % Area Ratio for row 2
vorticity_filtered4_row2 = Type4(2, 10);  % Vorticity for row 2

% Plot the first row point
plot(area_ratio_filtered4_row1, vorticity_filtered4_row1, 'h', 'Color', 'c', ...
     'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Type 4 (ANY_007)');

% Plot the second row point
plot(area_ratio_filtered4_row2, vorticity_filtered4_row2, '*', 'Color', 'k', ...
     'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Type 4 (ANY_111)');

% Update the legend
legend('Location', 'best');
% --- Distance Calculation for Type4 Points ---
% Fit line parameters
m = p(1); % Slope of the fit line
c = p(2); % Intercept of the fit line

% Type4 Data Points
x4_row1 = area_ratio_filtered4_row1; % Surface Area for ANY_007
y4_row1 = vorticity_filtered4_row1; % Low WSS Region for ANY_007

x4_row2 = area_ratio_filtered4_row2; % Surface Area for ANY_111
y4_row2 = vorticity_filtered4_row2; % Low WSS Region for ANY_111

% Predicted y-values on the fit line for Type4 points
yfit_row1 = m * x4_row1 + c;
yfit_row2 = m * x4_row2 + c;

% Calculate the vertical distances (residuals)
distance_row1 = abs(y4_row1 - yfit_row1);
distance_row2 = abs(y4_row2 - yfit_row2);

% Display the distances
fprintf('vorticity vs area ratio Distance of ANY_007 to the fit line: %.4f\n', distance_row1);
fprintf('vorticity vs area ratio Distance of ANY_111 to the fit line: %.4f\n', distance_row2);

% Optionally annotate the plot with the distances
text(x4_row1, y4_row1, sprintf('%.4f', distance_row1), 'FontSize', 10, 'Color', 'c', 'VerticalAlignment', 'bottom');
text(x4_row2, y4_row2, sprintf('%.4f', distance_row2), 'FontSize', 10, 'Color', 'k', 'VerticalAlignment', 'bottom');
% --- Helper function for filtering and sorting data ---
function [area_ratio_filtered, region_filtered] = process_data(Type, area_col, region_col)
    area_ratio = Type(:, area_col);  % Area Ratio (x-axis)
    region = Type(:, region_col);    % Region data (y-axis)

    % Filter out zero or NaN values in the region
    valid_indices = ~isnan(region) & region ~= 0 & area_ratio > 0;
    area_ratio_filtered = area_ratio(valid_indices);
    region_filtered = region(valid_indices);

    % Sort by area ratio
    [area_ratio_filtered, sort_idx] = sort(area_ratio_filtered);
    region_filtered = region_filtered(sort_idx);
end
