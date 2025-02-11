office_pc = 'C:\Users\67072\OneDrive - The Pennsylvania State University\Aneurysm CFD project\Post_export\';
close all;
% Define the case number dynamically as a string
case_number = '28';
save_dir = fullfile(office_pc, ['case' case_number]);
osi_mat_file = 'OSI1.mat'; % OSI file
load(fullfile(save_dir, osi_mat_file));

%{
x = osi_data(:, 1); % for case 1
y = osi_data(:, 2);
z = osi_data(:, 3);
osi_values = osi_data(:, 4);
%}

%
% Extract coordinates and OSI values (adjusting as per OSI data structure)
x = osi_data.x;
y = osi_data.y;
z = osi_data.z;
osi_values = osi_data.osi_mag;
%}
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
filtered_osi = osi_values(mask);

% Compute the average OSI for the masked data
average_osi = mean(filtered_osi);

% Set the high threshold as a percentage of the average OSI
high_threshold_percentage = 80;  % High threshold (percentage higher than the average)

% Calculate the actual high threshold value
high_threshold_value = (1 + high_threshold_percentage / 100) * average_osi;  

% Filter points that are above the high threshold
above_threshold_mask = filtered_osi > high_threshold_value;

% Display the average OSI and the percentage of points above the high threshold with the case number
fprintf('Average OSI for case %s: %.4f\n', case_number, average_osi);
fprintf('Percentage of points above %.0f%% higher than the average OSI for case %s: %.2f%%\n', ...
    high_threshold_percentage, case_number, (sum(above_threshold_mask) / length(filtered_osi)) * 100);

% Create a 3D scatter plot of the entire OSI data
figure;
hold on;

% Plot all OSI data points using the regular colorbar
scatter3(filtered_x, filtered_z, filtered_y, 100, filtered_osi, 'filled');

% Overlay the points with OSI > high threshold in red
scatter3(filtered_x(above_threshold_mask), filtered_z(above_threshold_mask), filtered_y(above_threshold_mask), ...
    100, filtered_osi(above_threshold_mask), 'filled', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');

% Add labels and colorbar
xlabel('X');
ylabel('Z');
zlabel('Y');
title(sprintf('3D Scatter Plot of OSI with High OSI Highlighted for Case %s', case_number));
colorbar;
clim([0, 0.25]);  % Adjust color limits if necessary
view([-90, 0]);
colormap('parula');
rotate3d on;

% Add a dynamic legend
legend('All OSI values', sprintf('High OSI (> %.0f%%)', 100 + high_threshold_percentage), ...
    'Location', 'best');

hold off;

scatter3(filtered_x, filtered_z, filtered_y, 100, filtered_osi, 'filled');

% Add labels and colorbar
axis off;
title(sprintf('3D Scatter Plot of Filtered OSI for Case %s', case_number));
colorbar;
clim([0, 0.5]);  
view([75, 30]);
colormap('jet');
rotate3d on;

hold off;

