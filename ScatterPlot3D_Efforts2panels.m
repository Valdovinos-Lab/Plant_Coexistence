% Main script
% Create a figure and set its properties
figure;
set(gcf, 'Color', 'w', 'Position', [100, 100, 1600, 950]); % Set figure background and dimensions

% Define common font sizes
axisLFont = 36; % Label font size
axisTFont = 32; % Tick font size

% Define common tick positions
xTicks = linspace(0, 60, 4); % Example: 4 evenly spaced ticks on x-axis
yTicks = linspace(0, 150, 4); % Example: 4 evenly spaced ticks on y-axis
zTicks = linspace(0, 1, 6);  % Example: 4 evenly spaced ticks on z-axis

%% Panel 1: First scatter plot
subplot(1, 2, 1); % Divide figure into 1 row, 2 columns, and use the 1st panel
scatter3(T_CoexistRes_Efforts.degreeA, T_CoexistRes_Efforts.degreeP, T_CoexistRes_Efforts.Effort0, 10, 'k', 'filled');

% Set axis labels and fonts
xlabel('Animal degree', 'FontSize', axisLFont);
ylabel('Plant degree', 'FontSize', axisLFont);
zlabel('Foraging effort (\alpha_{ij})', 'FontSize', axisLFont +2, 'Interpreter', 'tex');

set(gca, 'Color', 'w', 'FontSize', axisTFont); % Set axes background and tick font size
grid on;

% Adjust axis labels
adjustAxisLabels(gca, xTicks, yTicks, zTicks);

% Adjust subplot position to make it wider
set(gca, 'Position', [0.12, 0.1, 0.38, 0.8]); % [left, bottom, width, height]

%% Panel 2: Second scatter plot
subplot(1, 2, 2); % Use the 2nd panel
scatter3(T_CoexistRes_Efforts.degreeA, T_CoexistRes_Efforts.degreeP, T_CoexistRes_Efforts.Effort1, 10, 'k', 'filled');

% Set axis labels and fonts
%xlabel('Animal degree', 'FontSize', axisLFont);
%ylabel('Plant degree', 'FontSize', axisLFont);
%zlabel('Effort', 'FontSize', axisLFont);
set(gca, 'Color', 'w', 'FontSize', axisTFont); % Set axes background and tick font size
grid on;

% Adjust axis labels
adjustAxisLabels(gca, xTicks, yTicks, zTicks);

% Adjust subplot position to make it wider
set(gca, 'Position', [0.59, 0.1, 0.38, 0.8]); % [left, bottom, width, height]


%% Helper Function: Adjust Axis Labels
function adjustAxisLabels(ax, xTicks, yTicks, zTicks)
    % Rotate axis labels
    ax.XLabel.Rotation = 25; % Rotate X-axis label
    ax.YLabel.Rotation = -35; % Rotate Y-axis label
    ax.ZLabel.Rotation = 90; % Rotate Z-axis label
    
    % Adjust X-axis label position
    ax.XLabel.Position = [mean(xlim) - 15, ax.XLabel.Position(2), ax.XLabel.Position(3) + 0.055]; % Move right and closer

    % Adjust Y-axis label position
    ax.YLabel.Position = [ax.YLabel.Position(1), mean(ylim) - 55, ax.YLabel.Position(3) + 0.12]; % Move left and closer

    % Adjust Z-axis label position (optional, if needed)
    ax.ZLabel.Position = [ax.ZLabel.Position(1) - 3, ax.ZLabel.Position(2), mean(zlim)]; % Centered along Z-axis

    % Set axis ticks
    xticks(ax, xTicks);
    yticks(ax, yTicks);
    zticks(ax, zTicks);
end

%exportgraphics(gcf, 'Efforts_without&withAF.pdf', 'ContentType', 'vector'); % Save as PDF (vector format)
%exportgraphics(gcf, 'Efforts_without&withAF2.emf', 'ContentType', 'vector'); % Save as EMF (vector format)
exportgraphics(gcf, 'Efforts_without&withAF2.png', 'Resolution', 300); % Save as PNG with high resolution