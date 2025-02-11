close all;
filePath = 'C:\Users\67072\OneDrive - The Pennsylvania State University\Aneurysm CFD project\R square for heatmap.xlsx';

% Read the data from the Excel sheet
[num, txt, raw] = xlsread(filePath, 'Sheet1');

% Extract headers and data separately
hemodynamic_parameters = txt(2:end, 1); % Hemodynamic parameters (first column, starting from row 2)
morphology_parameters = txt(1, 2:end);  % Morphology parameters (first row, starting from column 2)
data = num;                             % Data values

% Create the heatmap without annotations
figure;
h = heatmap(morphology_parameters, hemodynamic_parameters, data);
% Customize colormap and labels
colormap('turbo');
h.Title = 'R-Squared Values Between Morphology and Hemodynamic';
h.XLabel = 'Morphology Parameter';
h.YLabel = 'Hemodynamic Parameter';
h.ColorLimits = [0 1]; 