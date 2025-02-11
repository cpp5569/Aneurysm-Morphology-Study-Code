%% COMPUTE AND PLOT RELATIVE RESIDENCE TIMES (RRT)
% RRT is computed as: RRT = 1/((1 - 2*OSI) * TAWSS)
clear;
% Define the base directory and file names
office_pc = 'C:\Users\67072\OneDrive - The Pennsylvania State University\Aneurysm CFD project\Post_export\';
case_number = '28';
save_dir = fullfile(office_pc, ['case' case_number]);
osi_mat_file = 'OSI.mat';     % OSI file
tawss_mat_file = 'TAWSS.mat'; % TAWSS file

% Load TAWSS data
load(fullfile(save_dir, tawss_mat_file));
tawss_x = tawss_data.tawss_x;
tawss_y = tawss_data.tawss_y;
tawss_z = tawss_data.tawss_z;
tawss_mag = tawss_data.tawss_mag;

% Load OSI data
load(fullfile(save_dir, osi_mat_file));
%
osi_avg = osi_data.osi_avg;  
osi_x = osi_data.osi_x;
osi_y = osi_data.osi_y;
osi_z = osi_data.osi_z;
%}
%{
osi_x = osi_data(:, 1); % for case1
osi_y = osi_data(:, 2);
osi_z = osi_data(:, 3);
osi_avg = osi_data(:, 4);
%}
% Compute RRT for each component
rrt_x = 1 ./ ((1 - 2 * osi_x) .* tawss_x);
rrt_y = 1 ./ ((1 - 2 * osi_y) .* tawss_y);
rrt_z = 1 ./ ((1 - 2 * osi_z) .* tawss_z);

% Compute the averaged RRT
rrt_avg = (rrt_x + rrt_y + rrt_z) / 3;

% Save the RRT results as a .mat file
rrt_file = 'RRT.mat';
save(fullfile(save_dir, rrt_file), 'rrt_avg', 'rrt_x', 'rrt_y', 'rrt_z');

% Display completion message
fprintf('RRT calculation completed and saved to %s\n', fullfile(save_dir, rrt_file));

%% PLOT RRT
% Assuming the data contains coordinates (X, Y, Z) in the same files

% Use coordinates from the loaded OSI data (X, Y, Z columns from the file)
x = osi_data.x;  % X-coordinates
y = osi_data.y;  % Y-coordinates
z = osi_data.z;  % Z-coordinates

% Plot the RRT magnitude (use rrt_avg for plotting)
figure;
scatter3(x, z, y, 10, rrt_avg, 'filled');  % 3D scatter plot with color-coded RRT
colorbar;  % Show color scale
colormap('jet');  % Use the 'jet' colormap
view([-90, 0]);

title('3D Plot of Relative Residence Time (RRT)');
xlabel('X-coordinate');
ylabel('Y-coordinate');
zlabel('Z-coordinate');
grid on;

% Set limits for color map if needed (optional)
clim([0, 30]);

% Enable 3D rotation for better visualization
rotate3d on;
