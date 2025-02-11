% Define the base directory and file name structure
close all;
office_pc = 'C:\Users\67072\OneDrive - The Pennsylvania State University\Aneurysm CFD project\Post_export\';
CFD_pc = 'C:\Users\67072\OneDrive - The Pennsylvania State University\Shared Documents - Brindise Lab\People\Cheng\Aneurysm CFD\Post_files\';
file_base_name = 'WSS-';
case_number = '28';
save_dir = fullfile(office_pc, ['case' case_number]);
base_dir = fullfile(CFD_pc, ['case' case_number]);

start_index = 505;
end_index = 600;
index_step = 5;
% Initialize variables to accumulate the sum of WSS components
wss_x_avg = [];
wss_y_avg = [];
wss_z_avg = [];
wss_x_absavg = [];
wss_y_absavg = [];
wss_z_absavg = [];

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
    
    %{
    header = fgetl(fid); % Skip header (for case1 to case9)
    data = fscanf(fid, '%f %f %f %f %f %f %f %f %f %f %f %f', [12 Inf]);
    data = data'; % Transpose to match the format
    fclose(fid);
    
    % Extract the relevant WSS components
    WSS_x = data(:, 10); % x-component of WSS
    WSS_y = data(:, 11); % y-component of WSS
    WSS_z = data(:, 12); % z-component of WSS
    WSS_mag = sqrt(WSS_x.^2 + WSS_y.^2 + WSS_z.^2); % Magnitude of WSS
    %}
    %
    header = fgetl(fid); % Skip header (for case10 to )
    data = fscanf(fid, '%f %f %f %f %f %f %f %f', [8 Inf]);
    data = data'; % Transpose to match the format
    fclose(fid);
    
    % Extract the relevant WSS components
    WSS_x = data(:, 6); % x-component of WSS
    WSS_y = data(:, 7); % y-component of WSS
    WSS_z = data(:, 8); % z-component of WSS
    WSS_mag = sqrt(WSS_x.^2 + WSS_y.^2 + WSS_z.^2); % Magnitude of WSS
    %}
    
    % Initialize accumulation variables on the first time step
    if curr_time == start_index
        wss_x_avg = zeros(size(WSS_x));
        wss_y_avg = zeros(size(WSS_y));
        wss_z_avg = zeros(size(WSS_z));
        wss_x_absavg = zeros(size(WSS_x));
        wss_y_absavg = zeros(size(WSS_y));
        wss_z_absavg = zeros(size(WSS_z));
    end
    
    % Accumulate the WSS components
    wss_x_avg = wss_x_avg + WSS_x;
    wss_y_avg = wss_y_avg + WSS_y;
    wss_z_avg = wss_z_avg + WSS_z;
    
    % Accumulate the absolute WSS components
    wss_x_absavg = wss_x_absavg + abs(WSS_x);
    wss_y_absavg = wss_y_absavg + abs(WSS_y);
    wss_z_absavg = wss_z_absavg + abs(WSS_z);
end

x = data(:, 2);
y = data(:, 3);
z = data(:, 4);

% Compute the OSI for each component
osi_x = 0.5 * (1 - abs(wss_x_avg) ./ wss_x_absavg);
osi_y = 0.5 * (1 - abs(wss_y_avg) ./ wss_y_absavg);
osi_z = 0.5 * (1 - abs(wss_z_avg) ./ wss_z_absavg);

% Compute the magnitude of the OSI vector
osi_mag = sqrt(osi_x.^2 + osi_y.^2 + osi_z.^2);

% Compute the average OSI across all directions
osi_avg = (osi_x + osi_y + osi_z) / 3;

% Save the OSI results to a file
output_file_name = fullfile(save_dir, 'OSI_WSS_computed.txt');
osi_data = [x, y, z, osi_mag, osi_x, osi_y, osi_z, osi_avg]; % X, Y, Z, OSI values
osi_table = array2table(osi_data, 'VariableNames', {'X', 'Y', 'Z', 'OSI_Magnitude', 'OSI_X', 'OSI_Y', 'OSI_Z', 'OSI_Avg'});
writetable(osi_table, output_file_name, 'Delimiter', '\t');
% Save the results as a .mat file
output_mat_file = fullfile(save_dir, 'OSI.mat');

% Create a structure to store all the data
osi_data = struct('x', x, 'y', y, 'z', z, 'osi_mag', osi_mag, 'osi_avg', osi_avg,...
                      'osi_x', osi_x, 'osi_y', osi_y, 'osi_z', osi_z);

% Save the structure as a .mat file
save(output_mat_file, 'osi_data');% Display a message indicating completion
fprintf('OSI calculated and saved to %s\n', output_file_name);

%% Plot OSI
figure;
scatter3(x, z, y, 10, osi_mag, 'filled');
colorbar; % Show color scale
colormap('jet')
clim([0,0.6]);
view([-90, 0]);
title('3D Plot of Time-Averaged Wall Shear Stress (TAWSS)');
xlabel('X-coordinate');
ylabel('Z-coordinate');
zlabel('Y-coordinate');
grid on;

% Enable 3D rotation
rotate3d on;
