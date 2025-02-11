% Define the directory and case number
clear;
close all;
office_pc = 'C:\Users\67072\OneDrive - The Pennsylvania State University\Aneurysm CFD project\Post_export\';
case_number = '28';  %
save_dir = fullfile(office_pc, ['case' case_number]);
reshaped_save_dir = fullfile(save_dir, 'reshaped_velocity_data');

% Load the vorticity data (with coordinates included as 3D arrays)
load(fullfile(save_dir, 'time_averaged_vorticity.mat'));  % Adjust path as needed
load(fullfile(reshaped_save_dir, 'reshaped_velocity_time_505.mat'))
% Reshape the 3D arrays (Xn, Yn, Zn, and vorticity values) into 1D vectors
x = reshape(Xn, [], 1); % X-coordinates reshaped to 1D vector
y = reshape(Yn, [], 1); % Y-coordinates reshaped to 1D vector
z = reshape(Zn, [], 1); % Z-coordinates reshaped to 1D vector
vorticity_values = reshape(vorticity_magnitude_avg, [], 1); % Vorticity magnitudes reshaped to 1D vector

% Define the mask type (1 = Ball Mask, 2 = Range Mask)
mask_type = 2;  % Set to 1 for ball mask, 2 for range mask

% Define reference point for the ball mask
x_ref = 0;  
y_ref = 0e-3;  
z_ref = 0;  

% Calculate the distance from the arbitrary reference point (for ball mask)
distances = sqrt((x - x_ref).^2 + (y - y_ref).^2 + (z - z_ref).^2);
ball_mask = distances <= 5e-3;

% Define the range mask conditions
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

% Apply the mask to the reshaped data
masked_x = x(mask);
masked_y = y(mask);
masked_z = z(mask);
masked_vorticity = vorticity_values(mask);

% Filter out points where the vorticity is zero
non_zero_mask = masked_vorticity ~= 0;
masked_x = masked_x(non_zero_mask);
masked_y = masked_y(non_zero_mask);
masked_z = masked_z(non_zero_mask);
masked_vorticity = masked_vorticity(non_zero_mask);

% Compute the volume element (assuming uniform grid of size 0.1mm = 1e-4 m)
grid_size = 1e-4;  % in meters
volume_element = grid_size^3;  % volume of each grid cell in cubic meters

% Compute the total volume of the masked region
total_volume = sum(non_zero_mask) * volume_element;

% Compute the total vorticity (sum of vorticity magnitudes * volume element)
total_vorticity = sum(masked_vorticity) * volume_element;

% Normalize the vorticity by the total volume
normalized_vorticity = total_vorticity / total_volume;

% Save the masked vorticity data and the normalized vorticity for the current case
save(fullfile(save_dir, 'normalized_vorticity_data.mat'), 'masked_x', 'masked_y', 'masked_z', 'masked_vorticity', 'normalized_vorticity');

% Display the normalized vorticity for this case
fprintf('Normalized Vorticity for case %s: %.4f\n', case_number, normalized_vorticity);

%% Debugging figure: Plot the masked geometry, excluding zero vorticity values
figure;
scatter3(masked_x, masked_z, masked_y, 10, masked_vorticity, 'filled');
xlabel('X');
ylabel('Z');
zlabel('Y');
title(sprintf('Masked Geometry for Case %s', case_number));
axis equal;
grid on;
