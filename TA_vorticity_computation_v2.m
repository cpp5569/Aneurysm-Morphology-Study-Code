% Define directories
clear;
close all;
office_pc = 'C:\Users\67072\OneDrive - The Pennsylvania State University\Aneurysm CFD project\Post_export\';
CFD_pc = 'C:\Users\67072\OneDrive - The Pennsylvania State University\Shared Documents - Brindise Lab\People\Cheng\Aneurysm CFD\Post_files\';
file_base_name = 'velocity-';

% Define the range of cases
case_start = 29;  % Start case number case18 19 rerun
case_end = 29;   % End case number

start_index = 505;  % Start time index
end_index = 600;    % End time index
index_step = 5;     % Step size between time points

for case_number = case_start:case_end
    % Define case-specific directories
    case_str = num2str(case_number);  % Convert case number to string
    save_dir = fullfile(office_pc, ['case' case_str]);
    base_dir = fullfile(CFD_pc, ['case' case_str]);

    % Create folders to save velocity gradient and reshaped data
    gradient_save_dir = fullfile(save_dir, 'velocity_gradient');
    if ~exist(gradient_save_dir, 'dir')
        mkdir(gradient_save_dir);
    end
    
    reshaped_save_dir = fullfile(save_dir, 'reshaped_velocity_data');
    if ~exist(reshaped_save_dir, 'dir')
        mkdir(reshaped_save_dir);
    end
    
    % Initialize variables for accumulating vorticity
    vorticity_x_total = 0;
    vorticity_y_total = 0;
    vorticity_z_total = 0;
    vorticity_magnitude_total = 0;
    time_steps_count = 0;

    x = [];
    y = [];
    z = [];
    
    for curr_time = start_index:index_step:end_index
        % Generate the file path for the current time step
        file_name = sprintf('%s%04d', file_base_name, curr_time);
        file_path = fullfile(base_dir, file_name);
        
        % Open and read the ASCII file
        fid = fopen(file_path, 'r');
        if fid == -1
            error(['File ' file_name ' could not be opened.']);
        end
        
        % Read the data with different formats for case 10
        if case_number == 10
            header = fgetl(fid); % Skip header (for case 10)
            data = fscanf(fid, '%f %f %f %f %f %f %f %f %f %f %f %f', [12 Inf]);
            data = data'; % Transpose to match the format
        else
            header = fgetl(fid); % Skip header
            data = fscanf(fid, '%f %f %f %f %f %f %f %f', [8 Inf]);
            data = data'; % Transpose to match the format
        end
        fclose(fid);
        
        % Extract the relevant columns for velocity and coordinates
        x_coords = data(:, 2);  % X-coordinate (column 2)
        y_coords = data(:, 3);  % Y-coordinate (column 3)
        z_coords = data(:, 4);  % Z-coordinate (column 4)
        
        Vel_x = data(:, 6);     % X-component of velocity (column 6)
        Vel_y = data(:, 7);     % Y-component of velocity (column 7)
        Vel_z = data(:, 8);     % Z-component of velocity (column 8)
        
        % Reshape the data into 3D grid
        if isempty(x)
            x = x_coords;
            y = y_coords;
            z = z_coords;
            
            % Calculate grid sizes
            grid_size_x = 1e-4; % Adjust based on your data
            grid_size_y = grid_size_x;
            grid_size_z = grid_size_x;

            % Ensure grid sizes are uniform across all axes
            if abs(grid_size_x - grid_size_y) > 1e-6 || abs(grid_size_y - grid_size_z) > 1e-6
                error('Grid sizes are not uniform');
            end

            % Determine grid dimensions
            x_min = min(x); x_max = max(x);
            y_min = min(y); y_max = max(y);
            z_min = min(z); z_max = max(z);
            
            Nx = round((x_max - x_min) / grid_size_x) + 1;
            Ny = round((y_max - y_min) / grid_size_y) + 1;
            Nz = round((z_max - z_min) / grid_size_z) + 1;
            
            % Create 3D grids for velocity components
            U_3D = zeros(Ny, Nx, Nz); % 3D grid for u (y, x, z)
            V_3D = zeros(Ny, Nx, Nz); % 3D grid for v (y, x, z)
            W_3D = zeros(Ny, Nx, Nz); % 3D grid for w (y, x, z)
        end
        
        % Fill the 3D grids with velocity values
        for i = 1:length(x_coords)
            ix = round((x_coords(i) - x_min) / grid_size_x) + 1;
            iy = round((y_coords(i) - y_min) / grid_size_y) + 1;
            iz = round((z_coords(i) - z_min) / grid_size_z) + 1;
            if ix >= 1 && ix <= Nx && iy >= 1 && iy <= Ny && iz >= 1 && iz <= Nz
                U_3D(iy, ix, iz) = Vel_x(i); % Swap x and y for MATLAB ordering
                V_3D(iy, ix, iz) = Vel_y(i);
                W_3D(iy, ix, iz) = Vel_z(i);
            end
        end
        
        % Create 3D coordinate grids using meshgrid
        [Xn, Yn, Zn] = meshgrid(linspace(x_min, x_max, Nx), ...
                                linspace(y_min, y_max, Ny), ...
                                linspace(z_min, z_max, Nz));

        % Save reshaped velocity data for the current time point
        reshaped_file = fullfile(reshaped_save_dir, ['reshaped_velocity_time_' num2str(curr_time) '.mat']);
        save(reshaped_file, 'U_3D', 'V_3D', 'W_3D', 'Xn', 'Yn', 'Zn');
        disp(['Reshaped velocity data saved for time point ' num2str(curr_time)]);

        % Compute velocity gradients using the reshaped data
        velmask = ones(size(U_3D)); % Mask for the entire domain
        winsize = 3; % Window size for gradient computation

        % Compute velocity gradients using gradient_rbf_3D function
        [dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz] = gradient_rbf_3D(Xn, Yn, Zn, U_3D, V_3D, W_3D, velmask, winsize);

        % Compute vorticity components
        vorticity_x = dwdy - dvdz; % omega_x = dw/dy - dv/dz
        vorticity_y = dudz - dwdx; % omega_y = du/dz - dw/dx
        vorticity_z = dvdx - dudy; % omega_z = dv/dx - du/dy

        % Compute vorticity magnitude
        vorticity_magnitude = sqrt(vorticity_x.^2 + vorticity_y.^2 + vorticity_z.^2);

        % Accumulate vorticity for time averaging
        vorticity_x_total = vorticity_x_total + vorticity_x;
        vorticity_y_total = vorticity_y_total + vorticity_y;
        vorticity_z_total = vorticity_z_total + vorticity_z;
        vorticity_magnitude_total = vorticity_magnitude_total + vorticity_magnitude;
        time_steps_count = time_steps_count + 1;

        % Save gradients and vorticity for the current time point
        gradient_file = fullfile(gradient_save_dir, ['velocity_gradient_time_' num2str(curr_time) '.mat']);
        save(gradient_file, 'dudx', 'dudy', 'dudz', 'dvdx', 'dvdy', 'dvdz', 'dwdx', 'dwdy', 'dwdz', 'vorticity_x', 'vorticity_y', 'vorticity_z', 'vorticity_magnitude');
        disp(['Gradient and vorticity saved for time point ' num2str(curr_time)]);
    end

    % Compute time-averaged vorticity
    vorticity_x_avg = vorticity_x_total / time_steps_count;
    vorticity_y_avg = vorticity_y_total / time_steps_count;
    vorticity_z_avg = vorticity_z_total / time_steps_count;
    vorticity_magnitude_avg = vorticity_magnitude_total / time_steps_count;

    % Save time-averaged vorticity
    average_vorticity_file = fullfile(save_dir, 'time_averaged_vorticity.mat');
    save(average_vorticity_file, 'vorticity_x_avg', 'vorticity_y_avg', 'vorticity_z_avg', 'vorticity_magnitude_avg');
    disp(['Time-averaged vorticity calculation done for case ' case_str]);

end
