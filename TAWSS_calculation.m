% Note: case1, case4 
clear;
close all;
office_pc = 'C:\Users\67072\OneDrive - The Pennsylvania State University\Aneurysm CFD project\Post_export\';
CFD_pc = 'C:\Users\67072\OneDrive - The Pennsylvania State University\Shared Documents - Brindise Lab\People\Cheng\Aneurysm CFD\Post_files\';
file_base_name = 'WSS-';
case_number = '28';
save_dir = fullfile(office_pc, ['case' case_number]);
base_dir = fullfile(CFD_pc, ['case' case_number]);

% Check if save directory exists, if not, create it
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
    fprintf('Directory %s did not exist, so it was created.\n', save_dir);
end

start_index = 505;  % Start time index
end_index = 600;   % End time index
index_step = 5;    % Step size between time points

Plot_WSS = 1; % Decide whether to plot WSS
% Initialize variables to accumulate the sum of WSS components
tawss_x = [];
tawss_y = [];
tawss_z = [];
tawss_mag = [];
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
    % Extract coordinates if they haven't been initialized
    if isempty(x)
        x = data(:, 2);  % X-coordinate
        y = data(:, 3);  % Y-coordinate
        z = data(:, 4);  % Z-coordinate
        % Initialize accumulation variables
        tawss_x = zeros(size(WSS_x));
        tawss_y = zeros(size(WSS_y));
        tawss_z = zeros(size(WSS_z));
        tawss_mag = zeros(size(WSS_mag));
    end
    
    % Accumulate the WSS components (time-averaged)
    tawss_x = tawss_x + abs(WSS_x);
    tawss_y = tawss_y + abs(WSS_y);
    tawss_z = tawss_z + abs(WSS_z);
    tawss_mag = tawss_mag + WSS_mag;
end

% Compute the average by dividing by the number of time steps
num_t = (end_index - start_index) / index_step + 1;
tawss_x = tawss_x / num_t;
tawss_y = tawss_y / num_t;
tawss_z = tawss_z / num_t;
tawss_mag = tawss_mag / num_t;

% Save the TAWSS results to a file
output_file_name = fullfile(save_dir, 'TAWSS_computed.txt');
tawss_data = [x, y, z, tawss_mag, tawss_x, tawss_y, tawss_z]; % X, Y, Z, TAWSS values
tawss_table = array2table(tawss_data, 'VariableNames', {'X', 'Y', 'Z', 'TAWSS_Magnitude', 'TAWSS_X', 'TAWSS_Y', 'TAWSS_Z'});
writetable(tawss_table, output_file_name, 'Delimiter', '\t');

% Save the results as a .mat file
output_mat_file = fullfile(save_dir, 'TAWSS.mat');

% Create a structure to store all the data
tawss_data = struct('x', x, 'y', y, 'z', z, 'tawss_mag', tawss_mag, ...
                      'tawss_x', tawss_x, 'tawss_y', tawss_y, 'tawss_z', tawss_z);

% Save the structure as a .mat file
save(output_mat_file, 'tawss_data');

% Display a message indicating completion
fprintf('TAWSS calculated and saved to %s\n', output_file_name);

%% Plot TAWSS
if Plot_WSS
    load(fullfile(save_dir, 'TAWSS.mat'));
    x = tawss_data.x;
    y = tawss_data.y;
    z = tawss_data.z;
    tawss_mag = tawss_data.tawss_mag;
    
    figure;
    scatter3(x, z, y, 10, tawss_mag, 'filled');
    colorbar; % Show color scale
    colormap('jet');
    clim([0, 2]);  % Adjust the color limits dynamically based on the range of TAWSS
    view([-90, 0]);
    title('3D Plot of Time-Averaged Wall Shear Stress (TAWSS)');
    xlabel('X-coordinate');
    ylabel('Y-coordinate');
    zlabel('Z-coordinate');
    grid on;
    
    % Enable 3D rotation
    rotate3d on;
end
