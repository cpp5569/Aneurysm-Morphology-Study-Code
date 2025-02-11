% Define the base directory and file name structure
clear;
close all;
office_pc = 'C:\Users\67072\OneDrive - The Pennsylvania State University\Aneurysm CFD project\Post_export\';
CFD_pc = 'C:\Users\67072\OneDrive - The Pennsylvania State University\Shared Documents - Brindise Lab\People\Cheng\Aneurysm CFD\Post_files\';
file_base_name = 'velocity-';
case_number = '27';

start_index = 505;  % Start time index
end_index = 600;   % End time index
index_step = 5;    % Step size between time points
save_dir = fullfile(office_pc, ['case' case_number]);
base_dir = fullfile(CFD_pc, ['case' case_number]);
% Initialize variables to accumulate the sum of velocity components
velocity_x = [];
velocity_y = [];
velocity_z = [];
x = [];
y = [];
z = [];

% Loop through the specified file range (time points)
for curr_time = start_index:index_step:end_index
    % Generate the file path for the current time step
    file_name = sprintf('%s%04d', file_base_name, curr_time);
    file_path = fullfile(base_dir, file_name);
    
    % Open and read the file
    fid = fopen(file_path, 'r');
    if fid == -1
        error(['File ' file_name ' could not be opened.']);
    end
    
    % Read the data (assuming the format is consistent across all files)
    %
    header = fgetl(fid); % Skip header
    data = fscanf(fid, '%f %f %f %f %f %f %f %f', [8 Inf]);
    data = data'; % Transpose to match the format
    fclose(fid);
    %}
    %{
    header = fgetl(fid); % Skip header (for case10)
    data = fscanf(fid, '%f %f %f %f %f %f %f %f %f %f %f %f', [12 Inf]);
    data = data'; % Transpose to match the format
    fclose(fid);
    %}
    % Extract the relevant columns for velocity and coordinates
    x_coords = data(:, 2);  % X-coordinate (column 2)
    y_coords = data(:, 3);  % Y-coordinate (column 3)
    z_coords = data(:, 4);  % Z-coordinate (column 4)
    
    Vel_x = data(:, 6);     % X-component of velocity (column 6)
    Vel_y = data(:, 7);     % Y-component of velocity (column 7)
    Vel_z = data(:, 8);     % Z-component of velocity (column 8)
    
    % Initialize accumulation variables on the first time step
    if isempty(x)
        x = x_coords;
        y = y_coords;
        z = z_coords;
        velocity_x = zeros(size(Vel_x));
        velocity_y = zeros(size(Vel_y));
        velocity_z = zeros(size(Vel_z));
    end
    
    % Accumulate the velocity components (sum over time steps)
    velocity_x = velocity_x + Vel_x;
    velocity_y = velocity_y + Vel_y;
    velocity_z = velocity_z + Vel_z;
end

% Compute the average by dividing by the number of time steps
num_t = (end_index - start_index) / index_step + 1;
u = velocity_x / num_t;
v = velocity_y / num_t;
w = velocity_z / num_t;

% Compute the magnitude of the time-averaged velocity
velocity_mag = sqrt(u.^2 + v.^2 + w.^2);

% Save the time-averaged velocity results to a file
output_file_name = fullfile(save_dir, 'TAVEL.txt');
velocity_data = [x, y, z, velocity_mag, velocity_x, velocity_y, velocity_z]; % X, Y, Z, Velocity values
velocity_table = array2table(velocity_data, 'VariableNames', {'X', 'Y', 'Z', 'Velocity_Magnitude', 'Velocity_X', 'Velocity_Y', 'Velocity_Z'});
writetable(velocity_table, output_file_name, 'Delimiter', '\t');

% Save the results as a .mat file
save(fullfile(save_dir, 'TAVEL.mat'), 'x', 'y', 'z', 'velocity_mag', 'u', 'v', 'w');

% Display a message indicating completion
fprintf('Time-averaged velocity calculated and saved to %s\n', output_file_name);

%quiver3(x, y, z, u, v, w);
%% Plot Time-Averaged Velocity Magnitude
figure;
scatter3(x, y, z, 10, velocity_mag, 'filled');
colorbar; % Show color scale
colormap('jet');
clim([0, max(velocity_mag)]);  % Adjust the color limits dynamically based on the range of velocities
view([180,-90]);
title('3D Plot of Time-Averaged Velocity Magnitude');
xlabel('X-coordinate');
ylabel('Y-coordinate');
zlabel('Z-coordinate');
grid on;

% Enable 3D rotation
rotate3d on;
