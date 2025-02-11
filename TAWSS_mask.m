clear;
close all;
office_pc = 'C:\Users\67072\OneDrive - The Pennsylvania State University\Aneurysm CFD project\Post_export\';
% Define the case number dynamically as a string
case_number = '28';
save_dir = fullfile(office_pc, ['case' case_number]);
tawss_mat_file = 'TAWSS.mat'; % TAWSS file
load(fullfile(save_dir, tawss_mat_file));

%{
x = tawss_data(:, 1); % for case 1
y = tawss_data(:, 2);
z = tawss_data(:, 3);
tawss_values = tawss_data(:, 4);
%}

% Extract coordinates and TAWSS values (beside case 1)

x = tawss_data.x;
y = tawss_data.y;
z = tawss_data.z;
tawss_values = tawss_data.tawss_mag;

% Define the mask type (1 = Ball Mask, 2 = Range Mask)
mask_type = 2;  % Set to 1 for ball mask, 2 for range mask

% Define reference point for the ball mask
x_ref = 0;  
y_ref = 0e-3;  
z_ref = 0;  

% Calculate the distance from the arbitrary reference point (for ball mask)
distances = sqrt((x - x_ref).^2 + (y - y_ref).^2 + (z - z_ref).^2);
ball_mask = distances <= 5e-3;

% Define the range mask conditions (updated with y and z logic switched)
z_range_min = -0.0086; z_range_max = 0.0086;
y_range_min = -0.0068;  % Only y above this value is kept
range_mask = (z >= z_range_min & z <= z_range_max) & (y >= y_range_min);

% Apply the selected mask
if mask_type == 1
    % Use ball mask
    mask = ball_mask;
    fprintf('Using ball mask for case %s.\n', case_number);
else
    % Use range mask
    mask = range_mask;
    fprintf('Using range mask for case %s.\n', case_number);
end

% Filter the data using the selected mask
filtered_x = x(mask);
filtered_y = y(mask);
filtered_z = z(mask);
filtered_tawss = tawss_values(mask);

% Compute the average TAWSS for the masked data (median or mean?)
average_tawss = mean(filtered_tawss);

% Set the low and high thresholds as percentages of the average TAWSS
low_threshold_percentage = 60;  % Low threshold (percentage lower than the average)
high_threshold_percentage = 60;  % High threshold (percentage higher than the average)

% Calculate the actual threshold values
low_threshold_value = (1 - low_threshold_percentage / 100) * average_tawss;  % 80% of the average
high_threshold_value = (1 + high_threshold_percentage / 100) * average_tawss;  % 130% of the average

% Filter points that are below the low threshold or above the high threshold
below_threshold_mask = filtered_tawss < low_threshold_value;
above_threshold_mask = filtered_tawss > high_threshold_value;
normal_mask = ~(below_threshold_mask | above_threshold_mask);  % Points in the normal range

% Calculate the percentage of points below the low threshold and above the high threshold
num_points_below_threshold = sum(below_threshold_mask);
num_points_above_threshold = sum(above_threshold_mask);
total_num_points = length(filtered_tawss);
percentage_below_threshold = (num_points_below_threshold / total_num_points) * 100;
percentage_above_threshold = (num_points_above_threshold / total_num_points) * 100;

% Display the average TAWSS, and the percentages of points below and above the thresholds with the case number
fprintf('Average TAWSS for case %s: %.4f\n', case_number, average_tawss);
fprintf('Percentage of points below %.0f%% lower than the average TAWSS for case %s: %.2f%%\n', ...
    low_threshold_percentage, case_number, percentage_below_threshold);
fprintf('Percentage of points above %.0f%% higher than the average TAWSS for case %s: %.2f%%\n', ...
    high_threshold_percentage, case_number, percentage_above_threshold);

% Create a 3D scatter plot of the filtered data
figure;
hold on;

% Plot the points with normal TAWSS (within the thresholds)
scatter3(filtered_x(normal_mask), filtered_z(normal_mask), filtered_y(normal_mask), ...
    100, filtered_tawss(normal_mask), 'filled');

% Plot the points with TAWSS < low threshold in blue
scatter3(filtered_x(below_threshold_mask), filtered_z(below_threshold_mask), filtered_y(below_threshold_mask), ...
    100, filtered_tawss(below_threshold_mask), 'filled', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');

% Plot the points with TAWSS > high threshold in red
scatter3(filtered_x(above_threshold_mask), filtered_z(above_threshold_mask), filtered_y(above_threshold_mask), ...
    100, filtered_tawss(above_threshold_mask), 'filled', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');

% Add labels and colorbar
xlabel('X');
ylabel('Z');
zlabel('Y');
title(sprintf('3D Scatter Plot of TAWSS for Case %s', case_number));
colorbar;
clim([0, 2]);  % Adjust color limits if necessary
view([-90, 0]);
colormap('parula');
rotate3d on;

% Add a dynamic legend with percentages in the labels
legend(sprintf('Normal TAWSS (%.2f%%)', 100 - percentage_below_threshold - percentage_above_threshold), ...
    sprintf('Low TAWSS (< %.0f%%, %.2f%%)', 100 - low_threshold_percentage, percentage_below_threshold), ...
    sprintf('High TAWSS (> %.0f%%, %.2f%%)', 100 + high_threshold_percentage, percentage_above_threshold), ...
    'Location', 'best');

hold off;

figure;
hold on;

% Plot the filtered points (all points after applying the mask)
scatter3(filtered_x, filtered_z, filtered_y, 100, filtered_tawss, 'filled');

% Add labels and colorbar
axis off;
title(sprintf('3D Scatter Plot of Filtered TAWSS for Case %s', case_number));
colorbar;
clim([0, 3]);  
view([75, 30]);
colormap('jet');
rotate3d on;

hold off;
