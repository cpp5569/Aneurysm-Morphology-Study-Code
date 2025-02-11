close all;
clear;
office_pc = 'C:\Users\67072\OneDrive - The Pennsylvania State University\Aneurysm CFD project\Post_export\';
% Define the case number dynamically as a string
case_number = '28';
save_dir = fullfile(office_pc, ['case' case_number]);
rrt_mat_file = 'RRT.mat'; % RRT file
osi_mat_file = 'OSI.mat'; % OSI file
load(fullfile(save_dir, rrt_mat_file));
load(fullfile(save_dir, osi_mat_file));
%{
x = osi_data(:, 1); % X coordinates for case1
y = osi_data(:, 2); % Y coordinates
z = osi_data(:, 3); % Z coordinates
rrt_values = rrt_avg; % RRT values
%}
%
x = osi_data.x; % for others
y = osi_data.y;
z = osi_data.z;
rrt_values = rrt_avg;
%}
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
filtered_rrt = rrt_values(mask);

% Compute the average RRT for the masked data
average_rrt = mean(filtered_rrt);

% Set the high threshold as a percentage of the average RRT
high_threshold_percentage = 60;  % High threshold (percentage higher than the average)

% Calculate the actual high threshold value
high_threshold_value = (1 + high_threshold_percentage / 100) * average_rrt;  

% Filter points that are above the high threshold
above_threshold_mask = filtered_rrt > high_threshold_value;

% Display the average RRT and the percentage of points above the high threshold with the case number
fprintf('Average RRT for case %s: %.4f\n', case_number, average_rrt);
fprintf('Percentage of points above %.0f%% higher than the average RRT for case %s: %.2f%%\n', ...
    high_threshold_percentage, case_number, (sum(above_threshold_mask) / length(filtered_rrt)) * 100);

% Create a 3D scatter plot of the entire RRT data
figure;
hold on;

% Plot all RRT data points using the regular colorbar
scatter3(filtered_x, filtered_z, filtered_y, 100, filtered_rrt, 'filled');

% Overlay the points with RRT > high threshold in red
scatter3(filtered_x(above_threshold_mask), filtered_z(above_threshold_mask), filtered_y(above_threshold_mask), ...
    100, filtered_rrt(above_threshold_mask), 'filled', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');

% Add labels and colorbar
xlabel('X');
ylabel('Z');
zlabel('Y');
title(sprintf('3D Scatter Plot of RRT with High RRT Highlighted for Case %s', case_number));
colorbar;
clim([0, 20]);  % Adjust color limits if necessary
view([-90, 0]);
colormap('parula');
rotate3d on;

% Add a dynamic legend
legend('All RRT values', sprintf('High RRT (> %.0f%%)', 100 + high_threshold_percentage), ...
    'Location', 'best');

hold off;

scatter3(filtered_x, filtered_z, filtered_y, 100, filtered_rrt, 'filled');

% Add labels and colorbar
axis off;
title(sprintf('3D Scatter Plot of Filtered RRT for Case %s', case_number));
colorbar;
clim([0, 30]);  
view([75, 30]);
colormap('jet');
rotate3d on;

hold off;

