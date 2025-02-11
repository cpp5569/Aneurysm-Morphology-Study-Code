%% Load the data
% Replace 'myData.mat' with the actual filename of your MAT file
load(['C:\Users\67072\OneDrive - The Pennsylvania State University\Shared Documents - Brindise Lab\' ...
    'Projects\Cerebral-Aneurysm\4-output_data\postprocessed_data\ANY007\postprocessing_files\any007_STB_full_resolution_TAWSS.mat']);
% Assuming 'tawss_mag' is your 3D TAWSS magnitude data
% First, flatten it into a single column vector if not already done
tawss_mag_col = tawss_mag(:);

% Filter out NaN and zero values
filtered_tawss = tawss_mag_col(~isnan(tawss_mag_col) & tawss_mag_col ~= 0);

% Calculate the mean TAWSS value from the filtered data
average_tawss = mean(filtered_tawss);

% Set the low and high thresholds as percentages relative to the mean TAWSS
low_threshold_percentage = 60;   % Percentage below mean (e.g., mean * (1 - 0.6) = 40% of mean)
high_threshold_percentage = 60;  % Percentage above mean (e.g., mean * (1 + 0.6) = 160% of mean)

% Calculate the actual threshold values
low_threshold_value = (1 - low_threshold_percentage / 100) * average_tawss;  
high_threshold_value = (1 + high_threshold_percentage / 100) * average_tawss;

% Create logical masks for below, above, and normal ranges
below_threshold_mask = filtered_tawss < low_threshold_value;
above_threshold_mask = filtered_tawss > high_threshold_value;
normal_mask = ~(below_threshold_mask | above_threshold_mask);

% Count the number of points in each category
num_points_below_threshold = sum(below_threshold_mask);
num_points_above_threshold = sum(above_threshold_mask);
num_points_normal = sum(normal_mask);

% Calculate total number of points
total_num_points = length(filtered_tawss);

% Calculate percentages
percentage_below_threshold = (num_points_below_threshold / total_num_points) * 100;
percentage_above_threshold = (num_points_above_threshold / total_num_points) * 100;
percentage_normal = (num_points_normal / total_num_points) * 100;

% Display results
fprintf('Average TAWSS (filtered): %.4f\n', average_tawss);
fprintf('Low threshold (%.0f%% below mean): %.4f\n', low_threshold_percentage, low_threshold_value);
fprintf('High threshold (%.0f%% above mean): %.4f\n', high_threshold_percentage, high_threshold_value);
fprintf('Points below threshold: %d (%.2f%%)\n', num_points_below_threshold, percentage_below_threshold);
fprintf('Points above threshold: %d (%.2f%%)\n', num_points_above_threshold, percentage_above_threshold);
fprintf('Points in normal range: %d (%.2f%%)\n', num_points_normal, percentage_normal);
