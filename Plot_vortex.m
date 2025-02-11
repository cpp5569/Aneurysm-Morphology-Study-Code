clear;
close all;
% Define the directory and case number
office_pc = 'C:\Users\67072\OneDrive - The Pennsylvania State University\Aneurysm CFD project\Post_export\';
case_number = '28';
save_dir = fullfile(office_pc, ['case' case_number]);
stl_dir = 'C:\Users\67072\OneDrive - The Pennsylvania State University\Aneurysm CFD project\Geometry_STL';  % Update to the correct STL directory

vortex_method = 'qcrit';  % Choose 'qcrit' for Q-criterion or 'lambda' for lambda-2 criterion
filename = fullfile(save_dir, ['vortex_1' vortex_method '.mat']);

% Load the vortex field (3D array)
load(filename, 'vortex_field');

% Load the reshaped velocity field to get the x, y, z coordinates (3D arrays)
reshaped_mat_file = 'reshaped_velocity_field1.mat';  % Coordinates file
velocity_data = load(fullfile(save_dir, reshaped_mat_file));
x = velocity_data.Xn;
y = velocity_data.Yn;
z = velocity_data.Zn;

%% Masking logic (Ball Mask and Range Mask)

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

%% Apply the mask directly to the 3D arrays
masked_vortex_field = vortex_field .* mask;  % Apply the mask to vortex_field
masked_x = x .* mask;  % Apply the mask to X coordinates
masked_y = y .* mask;  % Apply the mask to Y coordinates
masked_z = z .* mask;  % Apply the mask to Z coordinates

%% Unit Conversion (Meters to Millimeters)
masked_x = masked_x * 1000;  % Convert X from meters to millimeters
masked_y = masked_y * 1000;  % Convert Y from meters to millimeters
masked_z = masked_z * 1000;  % Convert Z from meters to millimeters

%% Swap the Y and Z axes for plotting
swapped_y = masked_z;  % Now Z becomes Y for plotting
swapped_z = masked_y;  % Now Y becomes Z for plotting

%% Isosurface Rendering using the swapped Y and Z coordinates
figure;
% Define the threshold value for isosurface 
if strcmp(vortex_method, 'qcrit')
    %threshold = 0.07 * max(masked_vortex_field(:));
    threshold = 1.4e4;
else
    threshold = -10000;
end
% Create the isosurface using the swapped coordinates
fv = isosurface(masked_x, swapped_y, swapped_z, masked_vortex_field, threshold);

% Plot the isosurface
patch(fv, 'FaceColor', 'blue', 'EdgeColor', 'none');
camlight; lighting phong;

% Add labels and title
xlabel('X-axis (mm)');
ylabel('Z-axis (mm)');  % This was Y before, now swapped to Z
zlabel('Y-axis (mm)');  % This was Z before, now swapped to Y
title(sprintf('Isosurface Rendering for %s Vortex Field, Case %s', vortex_method, case_number));
axis equal;
view([40, 20]);
zoom(1.2);
grid off;

%% Compute and Print Isosurface Volume

% Define the grid size in meters
grid_size = 1e-4;  % Grid size is 0.1 mm or 1e-4 m
voxel_volume = grid_size^3;  % Volume of a single grid cell

% Create a binary mask where the vortex field exceeds the threshold
binary_isosurface_mask = masked_vortex_field >= threshold;

% Count the number of grid cells that meet the isosurface condition
num_voxels_in_isosurface = sum(binary_isosurface_mask(:));

% Calculate the total volume of the isosurface
isosurface_volume = num_voxels_in_isosurface * voxel_volume * 1e9;

% Print the volume in cubic millimeters
fprintf('The isosurface volume for the %s method in case %s is %.6f cubic mm.\n', ...
    vortex_method, case_number, isosurface_volume);

%% Import STL file and plot as an outer shell
stl_file = fullfile(stl_dir, ['7mm_case' case_number '.stl']);  % Path to your STL file
stl_data = stlread(stl_file);  % Import the STL file

% Swap Y and Z coordinates for the STL geometry and convert to millimeters
stl_points = stl_data.Points;  % Convert STL vertices to mm
swapped_stl_points = stl_points(:, [1, 3, 2]);  % Swap Y and Z axes for STL

% Plot the outer shell with transparency
hold on;
patch('Faces', stl_data.ConnectivityList, 'Vertices', swapped_stl_points, ...
      'FaceColor', [0.5, 0.5, 0.5], 'EdgeColor', 'none', 'FaceAlpha', 0.2);  % Set transparency to 20%
hold off;

% Adjust the plot
camlight; lighting phong;
