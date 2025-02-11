close all;
% Load the STL file
stlFileName = '7mm_case23.stl'; % Replace with the path to your STL file
fv = stlread(stlFileName);

% Swap y and z components
swappedPoints = fv.Points;
swappedPoints(:, [2, 3]) = swappedPoints(:, [3, 2]); % Swapping y and z columns

% Plot the STL using trisurf with swapped axes
figure;
trisurf(fv.ConnectivityList, swappedPoints(:, 1), swappedPoints(:, 2), swappedPoints(:, 3), ...
    'FaceColor', [0.5, 0.5, 0.5], 'EdgeColor', 'none', 'FaceAlpha', 0.8);

% Set up lighting and viewing options
axis on;
xlabel('X');
ylabel('Z'); % Adjusted label since y and z are swapped
zlabel('Y'); % Adjusted label since y and z are swapped
%title('STL File Visualization with Swapped Y and Z Axes');
camlight;
lighting gouraud;
grid off;
rotate3d on;
view([-90, 0]);
zoom(1);
% Adding some options to improve the visualization
material shiny;