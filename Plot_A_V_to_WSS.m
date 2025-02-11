close all;
clear;
data1 = load('type1.mat');
data2 = load('type2.mat');
data3 = load('type3.mat');
data4 = load('type4.mat');
Type1 = data1.Type1;
Type2 = data2.Type2;
Type3 = data3.Type3;
Type4 = data4.Type4;

% Column Definitions:
% Column 1: Aspect Ratio (AR)
% Column 2: Low WSS Region (Low Wall Shear Stress Region)
% Column 3: High WSS Region (High Wall Shear Stress Region)
% Column 4: Average WSS (Average Wall Shear Stress)
% Column 5: High OSI Region (High Oscillatory Shear Index Region)
% Column 6: Average OSI (Average Oscillatory Shear Index)
% Column 7: High RRT Region (High Recirculation Residence Time Region)
% Column 8: Average RRT (Average Recirculation Residence Time)
% Column 9: Area Ratio (New x-axis variable for your latest plot)
% Column 10: Volume Normalized Vorticity
% Column 11: Volume
% Column 12: Surface Area
% Column 13: qcrit vortex
% ANY_007 = Basilar
% ANY_111 = ICA
% --- Low WSS Region vs Area ---
figure;

% Extract, filter, and sort relevant columns for Low WSS with Area
[area_filtered1, low_WSS_region_filtered1] = process_data(Type1, 12, 2); % Column 12 is Surface Area
[area_filtered2, low_WSS_region_filtered2] = process_data(Type2, 12, 2);
[area_filtered3, low_WSS_region_filtered3] = process_data(Type3, 12, 2);

hold on;

% Plot each type separately with different markers
plot(area_filtered1, low_WSS_region_filtered1, 'o', 'Color', 'b', 'MarkerSize', 8, 'LineWidth', 2);
plot(area_filtered2, low_WSS_region_filtered2, '^', 'Color', 'r', 'MarkerSize', 8, 'LineWidth', 2);
plot(area_filtered3, low_WSS_region_filtered3, 's', 'Color', 'm', 'MarkerSize', 8, 'LineWidth', 2);

% Combine all the data for fitting
area_combined = [area_filtered1; area_filtered2; area_filtered3];
low_WSS_combined = [low_WSS_region_filtered1; low_WSS_region_filtered2; low_WSS_region_filtered3];

% Sort combined data for fitting and plotting
[area_combined, sort_idx] = sort(area_combined);
low_WSS_combined = low_WSS_combined(sort_idx);

% Linear fit across all points
p = polyfit(area_combined, low_WSS_combined, 1);
yfit = polyval(p, area_combined);
plot(area_combined, yfit, '--k', 'LineWidth', 2);

% Calculate R-squared
SS_res = sum((low_WSS_combined - yfit).^2);
SS_tot = sum((low_WSS_combined - mean(low_WSS_combined)).^2);
R_squared = 1 - (SS_res / SS_tot);
text(mean(area_combined), max(low_WSS_combined), ['R^2 = ', num2str(R_squared, 2)], 'FontSize', 12, 'Color', 'k');

xlabel('Surface Area');
ylabel('Low WSS Region');
title('Low WSS Region vs Surface Area');
% legend({'Type 1', 'Type 2', 'Type 3', 'Fit Line'});
grid off;

area_filtered4_row1 = Type4(1, 12);
low_WSS_region_filtered4_row1 = Type4(1, 2);

area_filtered4_row2 = Type4(2, 12);
low_WSS_region_filtered4_row2 = Type4(2, 2);

% Plot Type4 points
plot(area_filtered4_row1, low_WSS_region_filtered4_row1, 'h', 'Color', 'c', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Type 4 (ANY_007)');
plot(area_filtered4_row2, low_WSS_region_filtered4_row2, '*', 'Color', 'k', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Type 4 (ANY_111)');
% --- Distance Calculation for Type4 Points ---
% Fit line parameters
m = p(1); % Slope of the fit line
c = p(2); % Intercept of the fit line

% Type4 Data Points
x4_row1 = area_filtered4_row1; % Surface Area for ANY_007
y4_row1 = low_WSS_region_filtered4_row1; % Low WSS Region for ANY_007

x4_row2 = area_filtered4_row2; % Surface Area for ANY_111
y4_row2 = low_WSS_region_filtered4_row2; % Low WSS Region for ANY_111

% Predicted y-values on the fit line for Type4 points
yfit_row1 = m * x4_row1 + c;
yfit_row2 = m * x4_row2 + c;

% Calculate the vertical distances (residuals)
distance_row1 = abs(y4_row1 - yfit_row1);
distance_row2 = abs(y4_row2 - yfit_row2);

% Display the distances
fprintf('low_WSS vs area Distance of ANY_007 to the fit line: %.4f\n', distance_row1);
fprintf('low_WSS vs area Distance of ANY_111 to the fit line: %.4f\n', distance_row2);

% Optionally annotate the plot with the distances
text(x4_row1, y4_row1, sprintf('%.4f', distance_row1), 'FontSize', 10, 'Color', 'c', 'VerticalAlignment', 'bottom');
text(x4_row2, y4_row2, sprintf('%.4f', distance_row2), 'FontSize', 10, 'Color', 'k', 'VerticalAlignment', 'bottom');

% --- Low WSS Region vs Volume ---
figure;

% Extract, filter, and sort relevant columns for Low WSS with Volume
[volume_filtered1, low_WSS_region_filtered1] = process_data(Type1, 11, 2); % Column 11 is Volume
[volume_filtered2, low_WSS_region_filtered2] = process_data(Type2, 11, 2);
[volume_filtered3, low_WSS_region_filtered3] = process_data(Type3, 11, 2);

hold on;

% Plot each type separately with different markers
plot(volume_filtered1, low_WSS_region_filtered1, 'o', 'Color', 'b', 'MarkerSize', 8, 'LineWidth', 2);
plot(volume_filtered2, low_WSS_region_filtered2, '^', 'Color', 'r', 'MarkerSize', 8, 'LineWidth', 2);
plot(volume_filtered3, low_WSS_region_filtered3, 's', 'Color', 'm', 'MarkerSize', 8, 'LineWidth', 2);

% Combine all the data for fitting
volume_combined = [volume_filtered1; volume_filtered2; volume_filtered3];
low_WSS_combined = [low_WSS_region_filtered1; low_WSS_region_filtered2; low_WSS_region_filtered3];

% Sort combined data for fitting and plotting
[volume_combined, sort_idx] = sort(volume_combined);
low_WSS_combined = low_WSS_combined(sort_idx);

% Linear fit across all points
p = polyfit(volume_combined, low_WSS_combined, 1);
yfit = polyval(p, volume_combined);
plot(volume_combined, yfit, '--k', 'LineWidth', 2);

% Calculate R-squared
SS_res = sum((low_WSS_combined - yfit).^2);
SS_tot = sum((low_WSS_combined - mean(low_WSS_combined)).^2);
R_squared = 1 - (SS_res / SS_tot);
%text(mean(volume_combined), max(low_WSS_combined), ['R^2 = ', num2str(R_squared, 2)], 'FontSize', 12, 'Color', 'k');

xlabel('Volume');
ylabel('Low WSS Region');
title('Low WSS Region vs Volume');
% legend({'Type 1', 'Type 2', 'Type 3', 'Fit Line'});
grid off;
ylim([0 45]);
% Load and add Type4 data
volume_filtered4_row1 = Type4(1, 11);
low_WSS_region_filtered4_row1 = Type4(1, 2);

volume_filtered4_row2 = Type4(2, 11);
low_WSS_region_filtered4_row2 = Type4(2, 2);

% Plot Type4 points
plot(volume_filtered4_row1, low_WSS_region_filtered4_row1, 'h', 'Color', 'c', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Type 4 (ANY_007)');
plot(volume_filtered4_row2, low_WSS_region_filtered4_row2, '*', 'Color', 'k', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Type 4 (ANY_111)');
% --- Distance Calculation for Type4 Points ---
% Fit line parameters
m = p(1); % Slope of the fit line
c = p(2); % Intercept of the fit line

% Type4 Data Points
x4_row1 = volume_filtered4_row1; % Surface Area for ANY_007
y4_row1 = low_WSS_region_filtered4_row1; % Low WSS Region for ANY_007

x4_row2 = volume_filtered4_row2; % Surface Area for ANY_111
y4_row2 = low_WSS_region_filtered4_row2; % Low WSS Region for ANY_111

% Predicted y-values on the fit line for Type4 points
yfit_row1 = m * x4_row1 + c;
yfit_row2 = m * x4_row2 + c;

% Calculate the vertical distances (residuals)
distance_row1 = abs(y4_row1 - yfit_row1);
distance_row2 = abs(y4_row2 - yfit_row2);

% Display the distances
fprintf('low_WSS vs volume Distance of ANY_007 to the fit line: %.4f\n', distance_row1);
fprintf('low_WSS vs volume Distance of ANY_111 to the fit line: %.4f\n', distance_row2);

% Optionally annotate the plot with the distances
text(x4_row1, y4_row1, sprintf('%.4f', distance_row1), 'FontSize', 10, 'Color', 'c', 'VerticalAlignment', 'bottom');
text(x4_row2, y4_row2, sprintf('%.4f', distance_row2), 'FontSize', 10, 'Color', 'k', 'VerticalAlignment', 'bottom');

% --- High WSS Region vs Area ---
figure;

% Extract, filter, and sort relevant columns for High WSS with Area
[area_filtered1, high_WSS_region_filtered1] = process_data(Type1, 12, 3); % Column 12 is Surface Area
[area_filtered2, high_WSS_region_filtered2] = process_data(Type2, 12, 3);
[area_filtered3, high_WSS_region_filtered3] = process_data(Type3, 12, 3);

hold on;

% Plot each type separately with different markers
plot(area_filtered1, high_WSS_region_filtered1, 'o', 'Color', 'b', 'MarkerSize', 8, 'LineWidth', 2);
plot(area_filtered2, high_WSS_region_filtered2, '^', 'Color', 'r', 'MarkerSize', 8, 'LineWidth', 2);
plot(area_filtered3, high_WSS_region_filtered3, 's', 'Color', 'm', 'MarkerSize', 8, 'LineWidth', 2);

% Combine all the data for fitting
area_combined = [area_filtered1; area_filtered2; area_filtered3];
high_WSS_combined = [high_WSS_region_filtered1; high_WSS_region_filtered2; high_WSS_region_filtered3];

% Sort combined data for fitting and plotting
[area_combined, sort_idx] = sort(area_combined);
high_WSS_combined = high_WSS_combined(sort_idx);

% Linear fit across all points
p = polyfit(area_combined, high_WSS_combined, 1);
yfit = polyval(p, area_combined);
plot(area_combined, yfit, '--k', 'LineWidth', 2);

% Calculate R-squared
SS_res = sum((high_WSS_combined - yfit).^2);
SS_tot = sum((high_WSS_combined - mean(high_WSS_combined)).^2);
R_squared = 1 - (SS_res / SS_tot);
text(mean(area_combined), max(high_WSS_combined), ['R^2 = ', num2str(R_squared, 2)], 'FontSize', 12, 'Color', 'k');

xlabel('Surface Area');
ylabel('High WSS Region');
title('High WSS Region vs Surface Area');
% legend({'Type 1', 'Type 2', 'Type 3', 'Fit Line'});
grid off;
% Load and add Type4 data
area_filtered4_row1 = Type4(1, 12);
high_WSS_region_filtered4_row1 = Type4(1, 3);

area_filtered4_row2 = Type4(2, 12);
high_WSS_region_filtered4_row2 = Type4(2, 3);

% Plot Type4 points
plot(area_filtered4_row1, high_WSS_region_filtered4_row1, 'h', 'Color', 'c', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Type 4 (ANY_007)');
plot(area_filtered4_row2, high_WSS_region_filtered4_row2, '*', 'Color', 'k', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Type 4 (ANY_111)');
% --- Distance Calculation for Type4 Points ---
% Fit line parameters
m = p(1); % Slope of the fit line
c = p(2); % Intercept of the fit line

% Type4 Data Points
x4_row1 = area_filtered4_row1; % Surface Area for ANY_007
y4_row1 = high_WSS_region_filtered4_row1; % Low WSS Region for ANY_007

x4_row2 = area_filtered4_row2; % Surface Area for ANY_111
y4_row2 = high_WSS_region_filtered4_row2; % Low WSS Region for ANY_111

% Predicted y-values on the fit line for Type4 points
yfit_row1 = m * x4_row1 + c;
yfit_row2 = m * x4_row2 + c;

% Calculate the vertical distances (residuals)
distance_row1 = abs(y4_row1 - yfit_row1);
distance_row2 = abs(y4_row2 - yfit_row2);

% Display the distances
fprintf('high_WSS vs area Distance of ANY_007 to the fit line: %.4f\n', distance_row1);
fprintf('high_WSS vs area Distance of ANY_111 to the fit line: %.4f\n', distance_row2);

% Optionally annotate the plot with the distances
text(x4_row1, y4_row1, sprintf('%.4f', distance_row1), 'FontSize', 10, 'Color', 'c', 'VerticalAlignment', 'bottom');
text(x4_row2, y4_row2, sprintf('%.4f', distance_row2), 'FontSize', 10, 'Color', 'k', 'VerticalAlignment', 'bottom');

% --- High WSS Region vs Volume ---
figure;

% Extract, filter, and sort relevant columns for High WSS with Volume
[volume_filtered1, high_WSS_region_filtered1] = process_data(Type1, 11, 3); % Column 11 is Volume
[volume_filtered2, high_WSS_region_filtered2] = process_data(Type2, 11, 3);
[volume_filtered3, high_WSS_region_filtered3] = process_data(Type3, 11, 3);

hold on;

% Plot each type separately with different markers
plot(volume_filtered1, high_WSS_region_filtered1, 'o', 'Color', 'b', 'MarkerSize', 8, 'LineWidth', 2);
plot(volume_filtered2, high_WSS_region_filtered2, '^', 'Color', 'r', 'MarkerSize', 8, 'LineWidth', 2);
plot(volume_filtered3, high_WSS_region_filtered3, 's', 'Color', 'm', 'MarkerSize', 8, 'LineWidth', 2);

% Combine all the data for fitting
volume_combined = [volume_filtered1; volume_filtered2; volume_filtered3];
high_WSS_combined = [high_WSS_region_filtered1; high_WSS_region_filtered2; high_WSS_region_filtered3];

% Sort combined data for fitting and plotting
[volume_combined, sort_idx] = sort(volume_combined);
high_WSS_combined = high_WSS_combined(sort_idx);

% Linear fit across all points
p = polyfit(volume_combined, high_WSS_combined, 1);
yfit = polyval(p, volume_combined);
plot(volume_combined, yfit, '--k', 'LineWidth', 2);

% Calculate R-squared
SS_res = sum((high_WSS_combined - yfit).^2);
SS_tot = sum((high_WSS_combined - mean(high_WSS_combined)).^2);
R_squared = 1 - (SS_res / SS_tot);
%text(mean(volume_combined), max(high_WSS_combined), ['R^2 = ', num2str(R_squared, 2)], 'FontSize', 12, 'Color', 'k');

xlabel('Volume');
ylabel('High WSS Region');
title('High WSS Region vs Volume');
% legend({'Type 1', 'Type 2', 'Type 3', 'Fit Line'});
grid off;
% Load and add Type4 data
volume_filtered4_row1 = Type4(1, 11);
high_WSS_region_filtered4_row1 = Type4(1, 3);

volume_filtered4_row2 = Type4(2, 11);
high_WSS_region_filtered4_row2 = Type4(2, 3);

% Plot Type4 points
plot(volume_filtered4_row1, high_WSS_region_filtered4_row1, 'h', 'Color', 'c', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Type 4 (ANY_007)');
plot(volume_filtered4_row2, high_WSS_region_filtered4_row2, '*', 'Color', 'k', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Type 4 (ANY_111)');
% --- Distance Calculation for Type4 Points ---
% Fit line parameters
m = p(1); % Slope of the fit line
c = p(2); % Intercept of the fit line

% Type4 Data Points
x4_row1 = volume_filtered4_row1; % Surface Area for ANY_007
y4_row1 = high_WSS_region_filtered4_row1; % Low WSS Region for ANY_007

x4_row2 = volume_filtered4_row2; % Surface Area for ANY_111
y4_row2 = high_WSS_region_filtered4_row2; % Low WSS Region for ANY_111

% Predicted y-values on the fit line for Type4 points
yfit_row1 = m * x4_row1 + c;
yfit_row2 = m * x4_row2 + c;

% Calculate the vertical distances (residuals)
distance_row1 = abs(y4_row1 - yfit_row1);
distance_row2 = abs(y4_row2 - yfit_row2);

% Display the distances
fprintf('high_WSS vs volume Distance of ANY_007 to the fit line: %.4f\n', distance_row1);
fprintf('high_WSS vs volume Distance of ANY_111 to the fit line: %.4f\n', distance_row2);

% Optionally annotate the plot with the distances
text(x4_row1, y4_row1, sprintf('%.4f', distance_row1), 'FontSize', 10, 'Color', 'c', 'VerticalAlignment', 'bottom');
text(x4_row2, y4_row2, sprintf('%.4f', distance_row2), 'FontSize', 10, 'Color', 'k', 'VerticalAlignment', 'bottom');

% --- Vortex vs Surface Area ---
figure;

% Extract, filter, and sort relevant columns for vortex with Surface Area
[area_filtered1, vortex_filtered1] = process_data(Type1, 12, 13); % Column 12 is Surface Area, Column 13 is qcrit vortex
[area_filtered2, vortex_filtered2] = process_data(Type2, 12, 13);
[area_filtered3, vortex_filtered3] = process_data(Type3, 12, 13);

% Plot the dots with different shapes
hold on;

% Plot Type 1 dots as circles (blue)
plot(area_filtered1, vortex_filtered1, 'o', 'Color', 'b', 'MarkerSize', 8, 'LineWidth', 2);

% Plot Type 2 dots as triangles (red)
plot(area_filtered2, vortex_filtered2, '^', 'Color', 'r', 'MarkerSize', 8, 'LineWidth', 2);

% Plot Type 3 dots as squares (magenta)
plot(area_filtered3, vortex_filtered3, 's', 'Color', 'm', 'MarkerSize', 8, 'LineWidth', 2);

% Combine the areas and vortex regions for fitting
area_combined = [area_filtered1; area_filtered2; area_filtered3];
vortex_combined = [vortex_filtered1; vortex_filtered2; vortex_filtered3];

% Sort combined data for fitting and plotting
[area_combined, sort_idx] = sort(area_combined);
vortex_combined = vortex_combined(sort_idx);

% Power fit (y = ax^b)
log_x = log(area_combined);
log_y = log(vortex_combined);

% Perform linear fit on log-transformed data
p = polyfit(log_x, log_y, 1);

% Extract the parameters for the power law: y = a * x^b
a = exp(p(2)); % Convert intercept back from log scale
b = p(1); % Slope is the power

% Generate the power fit curve
yfit = a * area_combined.^b;

% Plot the fit line
plot(area_combined, yfit, '--k', 'LineWidth', 2); % Black line for the power fit

% Calculate R-squared for the power fit
SS_res = sum((vortex_combined - yfit).^2); % Residual sum of squares
SS_tot = sum((vortex_combined - mean(vortex_combined)).^2); % Total sum of squares
R_squared = 1 - (SS_res / SS_tot); % R-squared value

% Display the R-squared value on the plot
text(mean(area_combined), max(vortex_combined), ['R^2 = ', num2str(R_squared, 2)], 'FontSize', 12, 'Color', 'k');

xlabel('Surface Area');
ylabel('qcrit Vortex');
title('qcrit Vortex vs Surface Area');
% legend({'Type 1', 'Type 2', 'Type 3', 'Power Fit'});
grid off;
% Load and add Type4 data
area_filtered4_row1 = Type4(1, 12);
vortex_filtered4_row1 = Type4(1, 13);

area_filtered4_row2 = Type4(2, 12);
vortex_filtered4_row2 = Type4(2, 13);

% Plot Type4 points
plot(area_filtered4_row1, vortex_filtered4_row1, 'h', 'Color', 'c', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Type 4 (ANY_007)');
plot(area_filtered4_row2, vortex_filtered4_row2, '*', 'Color', 'k', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Type 4 (ANY_111)');
% --- Distance Calculation for Type4 Points ---

% Power fit parameters
a = exp(p(2)); % Scale factor (intercept in log-log fit)
b = p(1);      % Exponent (slope in log-log fit)

% Type 4 Data Points
x4_row1 = area_filtered4_row1; % Independent variable (e.g., Surface Area for ANY_007)
y4_row1 = vortex_filtered4_row1; % Dependent variable (e.g., Low WSS Region for ANY_007)

x4_row2 = area_filtered4_row2; % Independent variable (e.g., Surface Area for ANY_111)
y4_row2 = vortex_filtered4_row2; % Dependent variable (e.g., Low WSS Region for ANY_111)

% Predicted y-values on the power fit line for Type 4 points
yfit_row1 = a * (x4_row1^b);
yfit_row2 = a * (x4_row2^b);

% Calculate the vertical distances (residuals)
distance_row1 = abs(y4_row1 - yfit_row1);
distance_row2 = abs(y4_row2 - yfit_row2);

% Display the distances
fprintf('vortex vs area Distance of ANY_007 to the power fit line: %.4f\n', distance_row1);
fprintf('vortex vs area Distance of ANY_111 to the power fit line: %.4f\n', distance_row2);

% Optionally annotate the plot with the distances
text(x4_row1, y4_row1, sprintf('%.4f', distance_row1), 'FontSize', 10, 'Color', 'c', 'VerticalAlignment', 'bottom');
text(x4_row2, y4_row2, sprintf('%.4f', distance_row2), 'FontSize', 10, 'Color', 'k', 'VerticalAlignment', 'bottom');


% --- Vortex vs Volume ---
figure;

% Extract, filter, and sort relevant columns for vortex with Volume
[volume_filtered1, vortex_filtered1] = process_data(Type1, 11, 13); % Column 11 is Volume, Column 13 is qcrit vortex
[volume_filtered2, vortex_filtered2] = process_data(Type2, 11, 13);
[volume_filtered3, vortex_filtered3] = process_data(Type3, 11, 13);

% Plot the dots with different shapes
hold on;

% Plot Type 1 dots as circles (blue)
plot(volume_filtered1, vortex_filtered1, 'o', 'Color', 'b', 'MarkerSize', 8, 'LineWidth', 2);

% Plot Type 2 dots as triangles (red)
plot(volume_filtered2, vortex_filtered2, '^', 'Color', 'r', 'MarkerSize', 8, 'LineWidth', 2);

% Plot Type 3 dots as squares (magenta)
plot(volume_filtered3, vortex_filtered3, 's', 'Color', 'm', 'MarkerSize', 8, 'LineWidth', 2);

% Combine the volumes and vortex regions for fitting
volume_combined = [volume_filtered1; volume_filtered2; volume_filtered3];
vortex_combined = [vortex_filtered1; vortex_filtered2; vortex_filtered3];

% Sort combined data for fitting and plotting
[volume_combined, sort_idx] = sort(volume_combined);
vortex_combined = vortex_combined(sort_idx);

% Power fit (y = ax^b)
log_x = log(volume_combined);
log_y = log(vortex_combined);

% Perform linear fit on log-transformed data
p = polyfit(log_x, log_y, 1);

% Extract the parameters for the power law: y = a * x^b
a = exp(p(2)); % Convert intercept back from log scale
b = p(1); % Slope is the power

% Generate the power fit curve
yfit = a * volume_combined.^b;

% Plot the fit line
plot(volume_combined, yfit, '--k', 'LineWidth', 2); % Black line for the power fit

% Calculate R-squared for the power fit
SS_res = sum((vortex_combined - yfit).^2); % Residual sum of squares
SS_tot = sum((vortex_combined - mean(vortex_combined)).^2); % Total sum of squares
R_squared = 1 - (SS_res / SS_tot); % R-squared value

% Display the R-squared value on the plot
text(mean(volume_combined), max(vortex_combined), ['R^2 = ', num2str(R_squared, 2)], 'FontSize', 12, 'Color', 'k');

xlabel('Volume');
ylabel('qcrit Vortex');
title('qcrit Vortex vs Volume');
% legend({'Type 1', 'Type 2', 'Type 3', 'Power Fit'});
grid off;
% Load and add Type4 data
volume_filtered4_row1 = Type4(1, 11);
vortex_filtered4_row1 = Type4(1, 13);

volume_filtered4_row2 = Type4(2, 11);
vortex_filtered4_row2 = Type4(2, 13);

% Plot Type4 points
plot(volume_filtered4_row1, vortex_filtered4_row1, 'h', 'Color', 'c', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Type 4 (ANY_007)');
plot(volume_filtered4_row2, vortex_filtered4_row2, '*', 'Color', 'k', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Type 4 (ANY_111)');

% --- Distance Calculation for Type4 Points ---

% Power fit parameters
a = exp(p(2)); % Scale factor (intercept in log-log fit)
b = p(1);      % Exponent (slope in log-log fit)

% Type 4 Data Points
x4_row1 = volume_filtered4_row1; % Independent variable (e.g., Surface Area for ANY_007)
y4_row1 = vortex_filtered4_row1; % Dependent variable (e.g., Low WSS Region for ANY_007)

x4_row2 = volume_filtered4_row2; % Independent variable (e.g., Surface Area for ANY_111)
y4_row2 = vortex_filtered4_row2; % Dependent variable (e.g., Low WSS Region for ANY_111)

% Predicted y-values on the power fit line for Type 4 points
yfit_row1 = a * (x4_row1^b);
yfit_row2 = a * (x4_row2^b);

% Calculate the vertical distances (residuals)
distance_row1 = abs(y4_row1 - yfit_row1);
distance_row2 = abs(y4_row2 - yfit_row2);

% Display the distances
fprintf('vortex vs volume Distance of ANY_007 to the power fit line: %.4f\n', distance_row1);
fprintf('vortex vs volume Distance of ANY_111 to the power fit line: %.4f\n', distance_row2);

% Optionally annotate the plot with the distances
text(x4_row1, y4_row1, sprintf('%.4f', distance_row1), 'FontSize', 10, 'Color', 'c', 'VerticalAlignment', 'bottom');
text(x4_row2, y4_row2, sprintf('%.4f', distance_row2), 'FontSize', 10, 'Color', 'k', 'VerticalAlignment', 'bottom');


% --- High OSI Region vs Surface Area ---
figure;

% Extract, filter, and sort relevant columns for High OSI with Surface Area
[area_filtered1, high_OSI_region_filtered1] = process_data(Type1, 12, 5); % Surface Area, High OSI Region
[area_filtered2, high_OSI_region_filtered2] = process_data(Type2, 12, 5);
[area_filtered3, high_OSI_region_filtered3] = process_data(Type3, 12, 5);

% Plot the dots with different shapes
hold on;

% Plot Type 1 dots as circles (blue)
plot(area_filtered1, high_OSI_region_filtered1, 'o', 'Color', 'b', 'MarkerSize', 8, 'LineWidth', 2);

% Plot Type 2 dots as triangles (red)
plot(area_filtered2, high_OSI_region_filtered2, '^', 'Color', 'r', 'MarkerSize', 8, 'LineWidth', 2);

% Plot Type 3 dots as squares (magenta)
plot(area_filtered3, high_OSI_region_filtered3, 's', 'Color', 'm', 'MarkerSize', 8, 'LineWidth', 2);

% Combine the areas and high OSI regions for fitting
area_combined = [area_filtered1; area_filtered2; area_filtered3];
high_OSI_combined = [high_OSI_region_filtered1; high_OSI_region_filtered2; high_OSI_region_filtered3];

% Sort combined data for fitting and plotting
[area_combined, sort_idx] = sort(area_combined);
high_OSI_combined = high_OSI_combined(sort_idx);

% Linear fit
p = polyfit(area_combined, high_OSI_combined, 1);
yfit = polyval(p, area_combined);
% Plot the fit line
plot(area_combined, yfit, '--k', 'LineWidth', 2); % Black dashed line for the fit

% Calculate R-squared
SS_res = sum((high_OSI_combined - yfit).^2); % Residual sum of squares
SS_tot = sum((high_OSI_combined - mean(high_OSI_combined)).^2); % Total sum of squares
R_squared = 1 - (SS_res / SS_tot);

% Display the R-squared value on the plot
%text(mean(area_combined), max(high_OSI_combined), ['R^2 = ', num2str(R_squared, 2)], 'FontSize', 12, 'Color', 'k');

xlabel('Surface Area');
ylabel('High OSI Region');
title('High OSI Region vs Surface Area');
% legend({'Type 1', 'Type 2', 'Type 3', 'Fit Line'});
grid off;
% Load and add Type4 data
area_filtered4_row1 = Type4(1, 12);
high_OSI_region_filtered4_row1 = Type4(1, 5);

area_filtered4_row2 = Type4(2, 12);
high_OSI_region_filtered4_row2 = Type4(2, 5);

% Plot Type4 points
plot(area_filtered4_row1, high_OSI_region_filtered4_row1, 'h', 'Color', 'c', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Type 4 (ANY_007)');
plot(area_filtered4_row2, high_OSI_region_filtered4_row2, '*', 'Color', 'k', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Type 4 (ANY_111)');

% --- Distance Calculation for Type4 Points ---
% Fit line parameters
m = p(1); % Slope of the fit line
c = p(2); % Intercept of the fit line

% Type4 Data Points
x4_row1 = area_filtered4_row1; % Surface Area for ANY_007
y4_row1 = high_OSI_region_filtered4_row1; % Low WSS Region for ANY_007

x4_row2 = area_filtered4_row2; % Surface Area for ANY_111
y4_row2 = high_OSI_region_filtered4_row2; % Low WSS Region for ANY_111

% Predicted y-values on the fit line for Type4 points
yfit_row1 = m * x4_row1 + c;
yfit_row2 = m * x4_row2 + c;

% Calculate the vertical distances (residuals)
distance_row1 = abs(y4_row1 - yfit_row1);
distance_row2 = abs(y4_row2 - yfit_row2);

% Display the distances
fprintf('high_OSI vs area Distance of ANY_007 to the fit line: %.4f\n', distance_row1);
fprintf('high_OSI vs area Distance of ANY_111 to the fit line: %.4f\n', distance_row2);

% Optionally annotate the plot with the distances
text(x4_row1, y4_row1, sprintf('%.4f', distance_row1), 'FontSize', 10, 'Color', 'c', 'VerticalAlignment', 'bottom');
text(x4_row2, y4_row2, sprintf('%.4f', distance_row2), 'FontSize', 10, 'Color', 'k', 'VerticalAlignment', 'bottom');
% --- High OSI Region vs Volume ---
figure;

% Extract, filter, and sort relevant columns for High OSI with Volume
[volume_filtered1, high_OSI_region_filtered1] = process_data(Type1, 11, 5); % Volume, High OSI Region
[volume_filtered2, high_OSI_region_filtered2] = process_data(Type2, 11, 5);
[volume_filtered3, high_OSI_region_filtered3] = process_data(Type3, 11, 5);

% Plot the dots with different shapes
hold on;

% Plot Type 1 dots as circles (blue)
plot(volume_filtered1, high_OSI_region_filtered1, 'o', 'Color', 'b', 'MarkerSize', 8, 'LineWidth', 2);

% Plot Type 2 dots as triangles (red)
plot(volume_filtered2, high_OSI_region_filtered2, '^', 'Color', 'r', 'MarkerSize', 8, 'LineWidth', 2);

% Plot Type 3 dots as squares (magenta)
plot(volume_filtered3, high_OSI_region_filtered3, 's', 'Color', 'm', 'MarkerSize', 8, 'LineWidth', 2);

% Combine the volumes and high OSI regions for fitting
volume_combined = [volume_filtered1; volume_filtered2; volume_filtered3];
high_OSI_combined = [high_OSI_region_filtered1; high_OSI_region_filtered2; high_OSI_region_filtered3];

% Sort combined data for fitting and plotting
[volume_combined, sort_idx] = sort(volume_combined);
high_OSI_combined = high_OSI_combined(sort_idx);

% Linear fit
p = polyfit(volume_combined, high_OSI_combined, 1);
yfit = polyval(p, volume_combined);

% Plot the fit line
plot(volume_combined, yfit, '--k', 'LineWidth', 2); % Black dashed line for the fit

% Calculate R-squared
SS_res = sum((high_OSI_combined - yfit).^2); % Residual sum of squares
SS_tot = sum((high_OSI_combined - mean(high_OSI_combined)).^2); % Total sum of squares
R_squared = 1 - (SS_res / SS_tot);

% Display the R-squared value on the plot
%text(mean(volume_combined), max(high_OSI_combined), ['R^2 = ', num2str(R_squared, 2)], 'FontSize', 12, 'Color', 'k');

xlabel('Volume');
ylabel('High OSI Region');
title('High OSI Region vs Volume');
% legend({'Type 1', 'Type 2', 'Type 3', 'Fit Line'});
grid off;
% Load and add Type4 data
volume_filtered4_row1 = Type4(1, 11);
high_OSI_region_filtered4_row1 = Type4(1, 5);

volume_filtered4_row2 = Type4(2, 11);
high_OSI_region_filtered4_row2 = Type4(2, 5);

% Plot Type4 points
plot(volume_filtered4_row1, high_OSI_region_filtered4_row1, 'h', 'Color', 'c', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Type 4 (ANY_007)');
plot(volume_filtered4_row2, high_OSI_region_filtered4_row2, '*', 'Color', 'k', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Type 4 (ANY_111)');

% --- Distance Calculation for Type4 Points ---
% Fit line parameters
m = p(1); % Slope of the fit line
c = p(2); % Intercept of the fit line

% Type4 Data Points
x4_row1 = volume_filtered4_row1; % Surface Area for ANY_007
y4_row1 = high_OSI_region_filtered4_row1; % Low WSS Region for ANY_007

x4_row2 = volume_filtered4_row2; % Surface Area for ANY_111
y4_row2 = high_OSI_region_filtered4_row2; % Low WSS Region for ANY_111

% Predicted y-values on the fit line for Type4 points
yfit_row1 = m * x4_row1 + c;
yfit_row2 = m * x4_row2 + c;

% Calculate the vertical distances (residuals)
distance_row1 = abs(y4_row1 - yfit_row1);
distance_row2 = abs(y4_row2 - yfit_row2);

% Display the distances
fprintf('high_OSI vs volume Distance of ANY_007 to the fit line: %.4f\n', distance_row1);
fprintf('high_OSI vs volume Distance of ANY_111 to the fit line: %.4f\n', distance_row2);

% Optionally annotate the plot with the distances
text(x4_row1, y4_row1, sprintf('%.4f', distance_row1), 'FontSize', 10, 'Color', 'c', 'VerticalAlignment', 'bottom');
text(x4_row2, y4_row2, sprintf('%.4f', distance_row2), 'FontSize', 10, 'Color', 'k', 'VerticalAlignment', 'bottom');
% --- High RRT Region vs Surface Area ---
figure;

% Extract, filter, and sort relevant columns for High RRT with Surface Area
[area_filtered1, high_RRT_region_filtered1] = process_data(Type1, 12, 7); % Surface Area, High RRT Region
[area_filtered2, high_RRT_region_filtered2] = process_data(Type2, 12, 7);
[area_filtered3, high_RRT_region_filtered3] = process_data(Type3, 12, 7);

% Plot the dots with different shapes
hold on;

% Plot Type 1 dots as circles (blue)
plot(area_filtered1, high_RRT_region_filtered1, 'o', 'Color', 'b', 'MarkerSize', 8, 'LineWidth', 2);

% Plot Type 2 dots as triangles (red)
plot(area_filtered2, high_RRT_region_filtered2, '^', 'Color', 'r', 'MarkerSize', 8, 'LineWidth', 2);

% Plot Type 3 dots as squares (magenta)
plot(area_filtered3, high_RRT_region_filtered3, 's', 'Color', 'm', 'MarkerSize', 8, 'LineWidth', 2);

% Combine the areas and high RRT regions for fitting
area_combined = [area_filtered1; area_filtered2; area_filtered3];
high_RRT_combined = [high_RRT_region_filtered1; high_RRT_region_filtered2; high_RRT_region_filtered3];

% Sort combined data for fitting and plotting
[area_combined, sort_idx] = sort(area_combined);
high_RRT_combined = high_RRT_combined(sort_idx);

% Linear fit
p = polyfit(area_combined, high_RRT_combined, 1);
yfit = polyval(p, area_combined);

% Plot the fit line
plot(area_combined, yfit, '--k', 'LineWidth', 2); % Black dashed line for the fit

% Calculate R-squared
SS_res = sum((high_RRT_combined - yfit).^2); % Residual sum of squares
SS_tot = sum((high_RRT_combined - mean(high_RRT_combined)).^2); % Total sum of squares
R_squared = 1 - (SS_res / SS_tot);

% Display the R-squared value on the plot
text(mean(area_combined), max(high_RRT_combined), ['R^2 = ', num2str(R_squared, 2)], 'FontSize', 12, 'Color', 'k');

xlabel('Surface Area');
ylabel('High RRT Region');
title('High RRT Region vs Surface Area');
% legend({'Type 1', 'Type 2', 'Type 3', 'Fit Line'});
grid off;
% Load and add Type4 data
area_filtered4_row1 = Type4(1, 12);
RRT_region_filtered4_row1 = Type4(1, 7);

area_filtered4_row2 = Type4(2, 12);
RRT_region_filtered4_row2 = Type4(2, 7);

% Plot Type4 points
plot(area_filtered4_row1, RRT_region_filtered4_row1, 'h', 'Color', 'c', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Type 4 (ANY_007)');
plot(area_filtered4_row2, RRT_region_filtered4_row2, '*', 'Color', 'k', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Type 4 (ANY_111)');

% --- Distance Calculation for Type4 Points ---
% Fit line parameters
m = p(1); % Slope of the fit line
c = p(2); % Intercept of the fit line

% Type4 Data Points
x4_row1 = area_filtered4_row1; % Surface Area for ANY_007
y4_row1 = RRT_region_filtered4_row1; % Low WSS Region for ANY_007

x4_row2 = area_filtered4_row2; % Surface Area for ANY_111
y4_row2 = RRT_region_filtered4_row2; % Low WSS Region for ANY_111

% Predicted y-values on the fit line for Type4 points
yfit_row1 = m * x4_row1 + c;
yfit_row2 = m * x4_row2 + c;

% Calculate the vertical distances (residuals)
distance_row1 = abs(y4_row1 - yfit_row1);
distance_row2 = abs(y4_row2 - yfit_row2);

% Display the distances
fprintf('high_RRT vs area Distance of ANY_007 to the fit line: %.4f\n', distance_row1);
fprintf('high_RRT vs area Distance of ANY_111 to the fit line: %.4f\n', distance_row2);

% Optionally annotate the plot with the distances
text(x4_row1, y4_row1, sprintf('%.4f', distance_row1), 'FontSize', 10, 'Color', 'c', 'VerticalAlignment', 'bottom');
text(x4_row2, y4_row2, sprintf('%.4f', distance_row2), 'FontSize', 10, 'Color', 'k', 'VerticalAlignment', 'bottom');
% --- High RRT Region vs Volume ---
figure;

% Extract, filter, and sort relevant columns for High RRT with Volume
[volume_filtered1, high_RRT_region_filtered1] = process_data(Type1, 11, 7); % Volume, High RRT Region
[volume_filtered2, high_RRT_region_filtered2] = process_data(Type2, 11, 7);
[volume_filtered3, high_RRT_region_filtered3] = process_data(Type3, 11, 7);

% Plot the dots with different shapes
hold on;

% Plot Type 1 dots as circles (blue)
plot(volume_filtered1, high_RRT_region_filtered1, 'o', 'Color', 'b', 'MarkerSize', 8, 'LineWidth', 2);

% Plot Type 2 dots as triangles (red)
plot(volume_filtered2, high_RRT_region_filtered2, '^', 'Color', 'r', 'MarkerSize', 8, 'LineWidth', 2);

% Plot Type 3 dots as squares (magenta)
plot(volume_filtered3, high_RRT_region_filtered3, 's', 'Color', 'm', 'MarkerSize', 8, 'LineWidth', 2);

% Combine the volumes and high RRT regions for fitting
volume_combined = [volume_filtered1; volume_filtered2; volume_filtered3];
high_RRT_combined = [high_RRT_region_filtered1; high_RRT_region_filtered2; high_RRT_region_filtered3];

% Sort combined data for fitting and plotting
[volume_combined, sort_idx] = sort(volume_combined);
high_RRT_combined = high_RRT_combined(sort_idx);

% Linear fit
p = polyfit(volume_combined, high_RRT_combined, 1);
yfit = polyval(p, volume_combined);

% Plot the fit line
plot(volume_combined, yfit, '--k', 'LineWidth', 2); % Black dashed line for the fit

% Calculate R-squared
SS_res = sum((high_RRT_combined - yfit).^2); % Residual sum of squares
SS_tot = sum((high_RRT_combined - mean(high_RRT_combined)).^2); % Total sum of squares
R_squared = 1 - (SS_res / SS_tot);

% Display the R-squared value on the plot
text(mean(volume_combined), max(high_RRT_combined), ['R^2 = ', num2str(R_squared, 2)], 'FontSize', 12, 'Color', 'k');

xlabel('Volume');
ylabel('High RRT Region');
title('High RRT Region vs Volume');
% legend({'Type 1', 'Type 2', 'Type 3', 'Fit Line'});
grid off;
% Load and add Type4 data
volume_filtered4_row1 = Type4(1, 11);
RRT_region_filtered4_row1 = Type4(1, 7);

volume_filtered4_row2 = Type4(2, 11);
RRT_region_filtered4_row2 = Type4(2, 7);

% Plot Type4 points
plot(volume_filtered4_row1, RRT_region_filtered4_row1, 'h', 'Color', 'c', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Type 4 (ANY_007)');
plot(volume_filtered4_row2, RRT_region_filtered4_row2, '*', 'Color', 'k', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Type 4 (ANY_111)');

% --- Distance Calculation for Type4 Points ---
% Fit line parameters
m = p(1); % Slope of the fit line
c = p(2); % Intercept of the fit line

% Type4 Data Points
x4_row1 = volume_filtered4_row1; % Surface Area for ANY_007
y4_row1 = RRT_region_filtered4_row1; % Low WSS Region for ANY_007

x4_row2 = volume_filtered4_row2; % Surface Area for ANY_111
y4_row2 = RRT_region_filtered4_row2; % Low WSS Region for ANY_111

% Predicted y-values on the fit line for Type4 points
yfit_row1 = m * x4_row1 + c;
yfit_row2 = m * x4_row2 + c;

% Calculate the vertical distances (residuals)
distance_row1 = abs(y4_row1 - yfit_row1);
distance_row2 = abs(y4_row2 - yfit_row2);

% Display the distances
fprintf('high_RRT vs volume Distance of ANY_007 to the fit line: %.4f\n', distance_row1);
fprintf('high_RRT vs volume Distance of ANY_111 to the fit line: %.4f\n', distance_row2);

% Optionally annotate the plot with the distances
text(x4_row1, y4_row1, sprintf('%.4f', distance_row1), 'FontSize', 10, 'Color', 'c', 'VerticalAlignment', 'bottom');
text(x4_row2, y4_row2, sprintf('%.4f', distance_row2), 'FontSize', 10, 'Color', 'k', 'VerticalAlignment', 'bottom');
% --- Vorticity vs Surface Area ---
figure;

% Extract, filter, and sort relevant columns for Vorticity with Surface Area
[area_filtered1, vorticity_filtered1] = process_data(Type1, 12, 10); % Surface Area, Vorticity
[area_filtered2, vorticity_filtered2] = process_data(Type2, 12, 10);
[area_filtered3, vorticity_filtered3] = process_data(Type3, 12, 10);

% Combine the data
area_combined = [area_filtered1; area_filtered2; area_filtered3];
vorticity_combined = [vorticity_filtered1; vorticity_filtered2; vorticity_filtered3];

% Sort the combined data by area
[area_combined, sort_idx] = sort(area_combined);
vorticity_combined = vorticity_combined(sort_idx);

% Plot data points for each type
hold on;

% Plot Type 1 dots as circles (blue)
plot(area_filtered1, vorticity_filtered1, 'o', 'Color', 'b', 'MarkerSize', 8, 'LineWidth', 2);

% Plot Type 2 dots as triangles (red)
plot(area_filtered2, vorticity_filtered2, '^', 'Color', 'r', 'MarkerSize', 8, 'LineWidth', 2);

% Plot Type 3 dots as squares (magenta)
plot(area_filtered3, vorticity_filtered3, 's', 'Color', 'm', 'MarkerSize', 8, 'LineWidth', 2);

% Add linear fit across all points
p = polyfit(area_combined, vorticity_combined, 1); % Linear fit
yfit = polyval(p, area_combined);

% Plot the fit line
plot(area_combined, yfit, '--k', 'LineWidth', 2); % Black dashed line for the fit

% Calculate R-squared
SS_res = sum((vorticity_combined - yfit).^2); % Residual sum of squares
SS_tot = sum((vorticity_combined - mean(vorticity_combined)).^2); % Total sum of squares
R_squared = 1 - (SS_res / SS_tot); % R-squared value

% Display the R-squared value on the plot
text(mean(area_combined), max(vorticity_combined), ['R^2 = ', num2str(R_squared, 2)], 'FontSize', 12, 'Color', 'k');

xlabel('Surface Area');
ylabel('Volume Normalized Vorticity');
title('Vorticity vs Surface Area');
ylim([100 300]);
%legend({'Type 1', 'Type 2', 'Type 3', 'Linear Fit'});
grid off;
% Load and add Type4 data
area_filtered4_row1 = Type4(1, 12);
vorticity_filtered4_row1 = Type4(1, 10);

area_filtered4_row2 = Type4(2, 12);
vorticity_filtered4_row2 = Type4(2, 10);

% Plot Type4 points
plot(area_filtered4_row1, vorticity_filtered4_row1, 'h', 'Color', 'c', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Type 4 (ANY_007)');
plot(area_filtered4_row2, vorticity_filtered4_row2, '*', 'Color', 'k', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Type 4 (ANY_111)');

% --- Distance Calculation for Type4 Points ---
% Fit line parameters
m = p(1); % Slope of the fit line
c = p(2); % Intercept of the fit line

% Type4 Data Points
x4_row1 = area_filtered4_row1; % Surface Area for ANY_007
y4_row1 = vorticity_filtered4_row1; % Low WSS Region for ANY_007

x4_row2 = area_filtered4_row2; % Surface Area for ANY_111
y4_row2 = vorticity_filtered4_row2; % Low WSS Region for ANY_111

% Predicted y-values on the fit line for Type4 points
yfit_row1 = m * x4_row1 + c;
yfit_row2 = m * x4_row2 + c;

% Calculate the vertical distances (residuals)
distance_row1 = abs(y4_row1 - yfit_row1);
distance_row2 = abs(y4_row2 - yfit_row2);

% Display the distances
fprintf('vorticity vs area Distance of ANY_007 to the fit line: %.4f\n', distance_row1);
fprintf('vorticity vs area Distance of ANY_111 to the fit line: %.4f\n', distance_row2);

% Optionally annotate the plot with the distances
text(x4_row1, y4_row1, sprintf('%.4f', distance_row1), 'FontSize', 10, 'Color', 'c', 'VerticalAlignment', 'bottom');
text(x4_row2, y4_row2, sprintf('%.4f', distance_row2), 'FontSize', 10, 'Color', 'k', 'VerticalAlignment', 'bottom');
% --- Vorticity vs Volume ---
figure;

% Extract, filter, and sort relevant columns for Vorticity with Volume
[volume_filtered1, vorticity_filtered1] = process_data(Type1, 11, 10); % Volume, Vorticity
[volume_filtered2, vorticity_filtered2] = process_data(Type2, 11, 10);
[volume_filtered3, vorticity_filtered3] = process_data(Type3, 11, 10);

% Combine the data
volume_combined = [volume_filtered1; volume_filtered2; volume_filtered3];
vorticity_combined = [vorticity_filtered1; vorticity_filtered2; vorticity_filtered3];

% Sort the combined data by volume
[volume_combined, sort_idx] = sort(volume_combined);
vorticity_combined = vorticity_combined(sort_idx);

% Plot data points for each type
hold on;

% Plot Type 1 dots as circles (blue)
plot(volume_filtered1, vorticity_filtered1, 'o', 'Color', 'b', 'MarkerSize', 8, 'LineWidth', 2);

% Plot Type 2 dots as triangles (red)
plot(volume_filtered2, vorticity_filtered2, '^', 'Color', 'r', 'MarkerSize', 8, 'LineWidth', 2);

% Plot Type 3 dots as squares (magenta)
plot(volume_filtered3, vorticity_filtered3, 's', 'Color', 'm', 'MarkerSize', 8, 'LineWidth', 2);

% Add linear fit across all points
p = polyfit(volume_combined, vorticity_combined, 1); % Linear fit
yfit = polyval(p, volume_combined);

% Plot the fit line
plot(volume_combined, yfit, '--k', 'LineWidth', 2); % Black dashed line for the fit

% Calculate R-squared
SS_res = sum((vorticity_combined - yfit).^2); % Residual sum of squares
SS_tot = sum((vorticity_combined - mean(vorticity_combined)).^2); % Total sum of squares
R_squared = 1 - (SS_res / SS_tot); % R-squared value

% Display the R-squared value on the plot
text(mean(volume_combined), max(vorticity_combined), ['R^2 = ', num2str(R_squared, 2)], 'FontSize', 12, 'Color', 'k');

xlabel('Volume');
ylabel('Volume Normalized Vorticity');
title('Vorticity vs Volume');
%legend({'Type 1', 'Type 2', 'Type 3', 'Linear Fit'});
ylim([100, 300]);
grid off;
% Load and add Type4 data
volume_filtered4_row1 = Type4(1, 11);
vorticity_filtered4_row1 = Type4(1, 10);

volume_filtered4_row2 = Type4(2, 11);
vorticity_filtered4_row2 = Type4(2, 10);

% Plot Type4 points
plot(volume_filtered4_row1, vorticity_filtered4_row1, 'h', 'Color', 'c', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Type 4 (ANY_007)');
plot(volume_filtered4_row2, vorticity_filtered4_row2, '*', 'Color', 'k', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Type 4 (ANY_111)');

% --- Distance Calculation for Type4 Points ---
% Fit line parameters
m = p(1); % Slope of the fit line
c = p(2); % Intercept of the fit line

% Type4 Data Points
x4_row1 = volume_filtered4_row1; % Surface Area for ANY_007
y4_row1 = vorticity_filtered4_row1; % Low WSS Region for ANY_007

x4_row2 = volume_filtered4_row2; % Surface Area for ANY_111
y4_row2 = vorticity_filtered4_row2; % Low WSS Region for ANY_111

% Predicted y-values on the fit line for Type4 points
yfit_row1 = m * x4_row1 + c;
yfit_row2 = m * x4_row2 + c;

% Calculate the vertical distances (residuals)
distance_row1 = abs(y4_row1 - yfit_row1);
distance_row2 = abs(y4_row2 - yfit_row2);

% Display the distances
fprintf('vorticity vs volume Distance of ANY_007 to the fit line: %.4f\n', distance_row1);
fprintf('vorticity vs volume Distance of ANY_111 to the fit line: %.4f\n', distance_row2);

% Optionally annotate the plot with the distances
text(x4_row1, y4_row1, sprintf('%.4f', distance_row1), 'FontSize', 10, 'Color', 'c', 'VerticalAlignment', 'bottom');
text(x4_row2, y4_row2, sprintf('%.4f', distance_row2), 'FontSize', 10, 'Color', 'k', 'VerticalAlignment', 'bottom');
% --- Helper function for filtering and sorting data ---
function [x_filtered, region_filtered] = process_data(Type, x_col, region_col)
    x = Type(:, x_col);  % New x-axis (Area or Volume)
    region = Type(:, region_col);    % Region data (y-axis remains the same)

    % Filter out zero values in the region
    valid_indices = region ~= 0;
    x_filtered = x(valid_indices);
    region_filtered = region(valid_indices);

    % Sort by x-axis variable (ascending order)
    [x_filtered, sort_idx] = sort(x_filtered);
    region_filtered = region_filtered(sort_idx);
end
