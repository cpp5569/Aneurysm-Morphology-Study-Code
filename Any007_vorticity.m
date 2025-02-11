% MATLAB code to calculate time-averaged velocity gradients, volume normalized vorticity, and isosurface volume

% Base file path for velocity gradient files
base_file_path = ['C:\Users\67072\OneDrive - The Pennsylvania State University\Shared Documents - Brindise Lab\Projects\Cerebral-Aneurysm\4-output_data\' ...
    'postprocessed_data\ANY007\postprocessing_files\velocity_gradients\any007_STB_full_resolution_velgrads_'];

% Initialize variables for cumulative velocity gradients and vorticity
cumulative_dudx = 0; cumulative_dudy = 0; cumulative_dudz = 0;
cumulative_dvdx = 0; cumulative_dvdy = 0; cumulative_dvdz = 0;
cumulative_dwdx = 0; cumulative_dwdy = 0; cumulative_dwdz = 0;
num_time_points = 365; % Total number of time points

% Iterate through all time points
for t = 1:num_time_points
    % Generate file name for current time point
    file_name = sprintf('%s%05d.mat', base_file_path, t);
    
    % Load velocity gradient components from the current file
    load(file_name, 'dudx', 'dudy', 'dudz', 'dvdx', 'dvdy', 'dvdz', 'dwdx', 'dwdy', 'dwdz');
    
    % Accumulate velocity gradient components
    cumulative_dudx = cumulative_dudx + dudx;
    cumulative_dudy = cumulative_dudy + dudy;
    cumulative_dudz = cumulative_dudz + dudz;
    cumulative_dvdx = cumulative_dvdx + dvdx;
    cumulative_dvdy = cumulative_dvdy + dvdy;
    cumulative_dvdz = cumulative_dvdz + dvdz;
    cumulative_dwdx = cumulative_dwdx + dwdx;
    cumulative_dwdy = cumulative_dwdy + dwdy;
    cumulative_dwdz = cumulative_dwdz + dwdz;
end

% Calculate time-averaged velocity gradients
time_averaged_dudx = cumulative_dudx / num_time_points;
time_averaged_dudy = cumulative_dudy / num_time_points;
time_averaged_dudz = cumulative_dudz / num_time_points;
time_averaged_dvdx = cumulative_dvdx / num_time_points;
time_averaged_dvdy = cumulative_dvdy / num_time_points;
time_averaged_dvdz = cumulative_dvdz / num_time_points;
time_averaged_dwdx = cumulative_dwdx / num_time_points;
time_averaged_dwdy = cumulative_dwdy / num_time_points;
time_averaged_dwdz = cumulative_dwdz / num_time_points;

% Use the vortex computation function to calculate vortex
vortex_method = 'qcrit';
[cal_mat] = vortex_computation_3D(vortex_method, time_averaged_dudx, time_averaged_dudy, time_averaged_dudz, ...
    time_averaged_dvdx, time_averaged_dvdy, time_averaged_dvdz, time_averaged_dwdx, time_averaged_dwdy, time_averaged_dwdz);

% Calculate the total grid number
% Assuming the grid size is fixed as 0.3mm^3, the total grid number is the size of one gradient component array
num_grids = numel(cal_mat);

% Calculate the total volume (grid size * total grid number)
grid_volume = 0.3^3; % Grid size in mm^3
total_volume = grid_volume * num_grids;

% Calculate volume normalized vortex
sum_vortex = sum(cal_mat(:));
volume_normalized_vortex = sum_vortex / total_volume;

% Define the threshold value for isosurface
if strcmp(vortex_method, 'qcrit')
    threshold = 1.4e4;
else
    threshold = -10000;
end

% Create a binary mask where the vortex field exceeds the threshold
binary_isosurface_mask = cal_mat >= threshold;

% Count the number of grid cells that meet the isosurface condition
num_voxels_in_isosurface = sum(binary_isosurface_mask(:));

% Calculate the total volume of the isosurface
isosurface_volume = num_voxels_in_isosurface * grid_volume * 1e2; % Volume in cubic mm

% Save the results to a .mat file
save('any007_time_averaged_vortex_results.mat', 'time_averaged_dudx', 'time_averaged_dudy', 'time_averaged_dudz', ...
    'time_averaged_dvdx', 'time_averaged_dvdy', 'time_averaged_dvdz', 'time_averaged_dwdx', 'time_averaged_dwdy', 'time_averaged_dwdz', 'cal_mat', 'volume_normalized_vortex', 'isosurface_volume');

% Display the results
fprintf('Volume normalized time-averaged vortex: %.6f\n', volume_normalized_vortex);
fprintf('Isosurface volume: %.6f cubic mm\n', isosurface_volume);
fprintf('Time-averaged vortex calculation complete. Results saved to time_averaged_vortex_results.mat\n');
