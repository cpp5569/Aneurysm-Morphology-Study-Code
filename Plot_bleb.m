clear;
close all;

load('aneurysm_bleb_data.mat');
% Define shape mapping for bleb location and color mapping for bleb size
shapes = {'o', 's', '*'};  % Circle, Square, Star for locations (left, middle, right)
colors = {'r', 'g', 'b'};  % Red, Green, Blue for sizes (small, middle, big)

% Function to plot a figure
function plot_figure(x_data, y_data, x_label_text, y_label_text, title_text, shapes, colors, surface_area, low_WSS)
    figure;
    hold on;

    % Plot Small Bleb (color: red) for different locations
    plot(x_data(1), y_data(1), ['r' shapes{1}], 'MarkerSize', 10);  
    plot(x_data(2), y_data(2), ['r' shapes{2}], 'MarkerSize', 10);  
    plot(x_data(3), y_data(3), ['r' shapes{3}], 'MarkerSize', 10);  

    % Plot Middle Bleb (color: green) for different locations
    plot(x_data(4), y_data(4), ['g' shapes{1}], 'MarkerSize', 10);  
    plot(x_data(5), y_data(5), ['g' shapes{2}], 'MarkerSize', 10);  
    plot(x_data(6), y_data(6), ['g' shapes{3}], 'MarkerSize', 10);  

    % Plot Big Bleb (color: blue) for different locations
    plot(x_data(7), y_data(7), ['b' shapes{1}], 'MarkerSize', 10);  
    plot(x_data(8), y_data(8), ['b' shapes{2}], 'MarkerSize', 10);  
    plot(x_data(9), y_data(9), ['b' shapes{3}], 'MarkerSize', 10);  

    % Add the black solid dot for the reference point
    plot(x_data(10), y_data(10), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k');  

    % Set labels and title
    xlabel(x_label_text);
    ylabel(y_label_text);
    title(title_text);

    hold off;
end

% Define X and Y axis data and labels for each plot
x_data_list = {surface_area, volume, area_ratio};
x_labels = {'Surface Area', 'Volume', 'Area Ratio'};
y_data_list = {low_WSS, high_WSS, high_OSI, high_RRT, vorticity, qcrit_vortex};
y_labels = {'Low WSS', 'High WSS', 'High OSI', 'High RRT', 'Vorticity', 'Qcrit Vortex'};
titles = {
    'Low WSS vs Surface Area', 'Low WSS vs Volume', 'Low WSS vs Area Ratio', ...
    'High WSS vs Surface Area', 'High WSS vs Volume', 'High WSS vs Area Ratio', ...
    'High OSI vs Surface Area', 'High OSI vs Volume', 'High OSI vs Area Ratio', ...
    'High RRT vs Surface Area', 'High RRT vs Volume', 'High RRT vs Area Ratio', ...
    'Vorticity vs Surface Area', 'Vorticity vs Volume', 'Vorticity vs Area Ratio', ...
    'Qcrit Vortex vs Surface Area', 'Qcrit Vortex vs Volume', 'Qcrit Vortex vs Area Ratio'
};

% Loop through all combinations of X and Y axis data to create plots
for y_idx = 1:length(y_data_list)
    for x_idx = 1:length(x_data_list)
        title_text = titles{(y_idx-1)*3 + x_idx};
        plot_figure(x_data_list{x_idx}, y_data_list{y_idx}, x_labels{x_idx}, y_labels{y_idx}, title_text, shapes, colors, surface_area, low_WSS);
    end
end
