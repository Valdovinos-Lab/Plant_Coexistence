% Read input data from file
data = readtable('Coex_AllPsp_muAP2_1200m.txt');

% Split data into nested and non-nested subsets based on N_lowHigh indicator
data_lowN = data(data.N_lowHigh == 0, :);  % non-nested networks
data_highN = data(data.N_lowHigh == 1, :); % nested networks

% Initialize main figure with specific dimensions
figWidth = 1200;   % Width in pixels
figHeight = 800;   % Height in pixels
fig = figure('Color', 'w', ...                    % White background
             'Position', [100 100 figWidth figHeight], ... % [left bottom width height]
             'Units', 'normalized');              % Use normalized units for consistent sizing

% Define subplot creation function for consistent panel formatting
function ax = createSubplot(position, yticks, ylimits)
    % Creates and formats a subplot with standardized appearance
    ax = subplot(2, 3, position, 'Position', getSubplotPosition(position));
    hold on;
    xlim([0 100]);
    % Apply custom y-axis ticks if provided
    if nargin > 1 && ~isempty(yticks)
        set(ax, 'YTick', yticks);
    end
    % Apply custom y-axis limits if provided
    if nargin > 2 && ~isempty(ylimits)
        ylim(ylimits);
    end
    % Set standard appearance properties
    set(ax, 'FontSize', 16, 'Box', 'on', 'LineWidth', 1);
    grid on;
end

% Define helper function for calculating subplot positions
function pos = getSubplotPosition(plotNumber)
    % Calculate precise positioning for each subplot
    % Parameters for panel layout
    leftMargin = 0.06;          % Distance from left edge
    rightMargin = 0.97;         % Distance from right edge
    topMargin = 1;             % Distance from top edge
    bottomMargin = 0.25;        % Distance from bottom edge
    horizontalSpacing = 0.089;  % Space between columns
    verticalSpacing = 0.1;     % Space between rows
    
    % Calculate dimensions for each panel
    panelWidth = (rightMargin - leftMargin - 2 * horizontalSpacing) / 3;
    panelHeight = (topMargin - bottomMargin - verticalSpacing) / 2;
    
    % Determine panel position based on its number
    row = ceil(plotNumber / 3);
    col = mod(plotNumber - 1, 3) + 1;
    
    % Calculate exact position coordinates
    x = leftMargin + (col - 1) * (panelWidth + horizontalSpacing);
    y = topMargin - row * (panelHeight + verticalSpacing);
    
    pos = [x y panelWidth panelHeight];
end

% Define plotting functions for data visualization
function plotScatterAndFit(x, y, markerColor, markerStyle, fitColor)
    % Creates scatter plot with linear fit
    scatter(x, y, 20, markerColor, markerStyle, 'filled');
    if ~isempty(x) && ~isempty(y)
        fitResult = fit(x, y, 'poly1');
        xFit = linspace(min(x), max(x), 100);
        yFit = feval(fitResult, xFit);
        h = plot(xFit, yFit, 'Color', fitColor, 'LineWidth', 2);
        uistack(h, 'top');
    end
end

function h = plotScatterOnly(x, y, markerColor, markerStyle)
    % Creates scatter plot without fit line and returns the scatter handle
    h = scatter(x, y, 20, markerColor, markerStyle, 'filled');
end

%% Create all subplots
% Panel 1: Demographic potential for non-nested networks
ax1 = createSubplot(1, [], [0 32]);
plotScatterAndFit(data_lowN.sVisits_perP0, data_lowN.fitness0, 'b', 'o', 'b');   % Without AF
plotScatterAndFit(data_lowN.sVisits_perP1, data_lowN.fitness1, [1 0.5 0], '^', [0.8 0.3 0]); % With AF

% Panel 2: Plant density for non-nested networks
ax2 = createSubplot(2, [], []);
plotScatterOnly(data_lowN.sVisits_perP0, data_lowN.abund0, 'b', 'o');           % Without AF
plotScatterOnly(data_lowN.sVisits_perP1, data_lowN.abund1, [1 0.5 0], '^');     % With AF

% Panel 3: Critical value for non-nested communities
ax3 = createSubplot(3, [], [0 0.05]);
% Capture the handles when plotting
h1 = plotScatterOnly(data_lowN.sVisits_perP0, data_lowN.sigma_c0, 'b', 'o');         % Without AF
h2 = plotScatterOnly(data_lowN.sVisits_perP1, data_lowN.sigma_c1, [1 0.5 0], '^');   % With AF

% Add legend to panel 3
legend([h1, h2], {'Without AF', 'With AF'}, ...
    'FontSize', 16, ...
    'Location', 'southeast', ...
    'Box', 'off');

% Panel 4: Demographic potential for nested communities
ax4 = createSubplot(4, [], [0 32]);
plotScatterAndFit(data_highN.sVisits_perP0, data_highN.fitness0, 'b', 'o', 'b');          % Without AF
plotScatterAndFit(data_highN.sVisits_perP1, data_highN.fitness1, [1 0.5 0], '^', [0.8 0.3 0]); % With AF

% Panel 5: Plant density for nested communities
ax5 = createSubplot(5, [], []);
plotScatterOnly(data_highN.sVisits_perP0, data_highN.abund0, 'b', 'o');         % Without AF
plotScatterOnly(data_highN.sVisits_perP1, data_highN.abund1, [1 0.5 0], '^');   % With AF

% Panel 6: Critical value for nested communities
ax6 = createSubplot(6, [], [0 0.05]);
plotScatterOnly(data_highN.sVisits_perP0, data_highN.sigma_c0, 'b', 'o');       % Without AF
plotScatterOnly(data_highN.sVisits_perP1, data_highN.sigma_c1, [1 0.5 0], '^'); % With AF

%% Format axes and add labels
% Set consistent x-axis tick marks for all panels
for ax = [ax1, ax2, ax3, ax4, ax5, ax6]
    set(ax, 'XTick', 0:20:100);
end

%% Add text annotations for axis labels and panel identifiers
% Add x-axis label (shared across bottom panels)
annotation('textbox', [0.35, 0.042, 0.3, 0.05], ...
    'String', 'Total visits per-plant', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', ...
    'FontSize', 22, ...
    'FitBoxToText', 'on');

% Add y-axis label for left panels (Demographic potential)
annotation('textbox', [0.12, 0.09, 0.1, 0.3], ...
    'String', 'Demographic potential', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', ...
    'FontSize', 22, ...
    'Rotation', 90, ...
    'FitBoxToText', 'on');

% Add y-axis label for middle panels (Plant density)
annotation('textbox', [0.425, 0.132, 0.05, 0.3], ...
    'String', 'Final density', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', ...
    'FontSize', 22, ...
    'Rotation', 90, ...
    'FitBoxToText', 'on');

% Add y-axis label for right panels (Critical value)
annotation('textbox', [0.797, 0.05, 0.05, 0.3], ...
    'String', 'Critical quality of visits', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', ...
    'FontSize', 22, ...
    'Rotation', 90, ...
    'FitBoxToText', 'on');

% Add labels to identify nested vs non-nested panels
% Label for top row (Non-nested communities)
annotation('textbox', [0.968, 0.8, 0.1, 0.1], ...
    'String', 'Non-Nested', ...
    'EdgeColor', 'none', ...
    'FontSize', 22, ...
    'Rotation', -90, ...
    'FitBoxToText', 'on');

% Label for bottom row (Nested communities)
annotation('textbox', [0.968, 0.34, 0.1, 0.1], ...
    'String', 'Nested', ...
    'EdgeColor', 'none', ...
    'FontSize', 22, ...
    'Rotation', -90, ...
    'FitBoxToText', 'on');

% exportgraphics(gcf, 'Fig2.pdf', 'ContentType', 'vector'); % Save as PDF (vector format)
exportgraphics(gcf, 'Fig2.emf', 'ContentType', 'vector'); % Save as EMF (vector format)
%exportgraphics(gcf, 'Fig2.png', 'Resolution', 300); % Save as PNG with high resolution