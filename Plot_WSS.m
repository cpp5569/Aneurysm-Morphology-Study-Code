% Define the base directory and file name structure
clear;
office_pc = 'C:\Users\67072\OneDrive - The Pennsylvania State University\Aneurysm CFD project\Post_export\';
CFD_pc = 'C:\Users\67072\OneDrive - The Pennsylvania State University\Shared Documents - Brindise Lab\People\Cheng\Aneurysm CFD\Post_files\';
file_base_name = 'WSS-';
case_number = '28';
save_dir = fullfile(office_pc, ['case' case_number]);
base_dir = fullfile(CFD_pc, ['case' case_number]);
% Generate the file path for the specified time point
file_name = sprintf('%s%04d', file_base_name, time_point);
file_path = fullfile(base_dir, file_name);

% Open and read the file
fid = fopen(file_path, 'r');
if fid == -1
    error(['File ' file_name ' could not be opened.']);
end

% Read the data (assuming the format is consistent across all files)
header = fgetl(fid); % Skip header
data = fscanf(fid, '%f %f %f %f %f %f %f %f', [8 Inf]);
data = data'; % Transpose to match the format
fclose(fid);

% Extract the relevant WSS components and coordinates
WSS_x = data(:, 6); % x-component of WSS
WSS_y = data(:, 7); % y-component of WSS
WSS_z = data(:, 8); % z-component of WSS
WSS_mag = sqrt(WSS_x.^2 + WSS_y.^2 + WSS_z.^2); % Magnitude of WSS

% Extract coordinates
x = data(:, 2);  % X-coordinate
y = data(:, 3);  % Y-coordinate
z = data(:, 4);  % Z-coordinate

% Plot the WSS magnitude for the specified time point
figure;
scatter3(x, y, z, 20, WSS_mag, 'filled');
colorbar;
colormap('jet');
clim([0,2]);
title(['WSS Magnitude at Time Point: ', num2str(time_point)]);
xlabel('X-coordinate');
ylabel('Y-coordinate');
zlabel('Z-coordinate');
grid on;
view([-90, 0]);

% Enable 3D rotation
rotate3d on;
