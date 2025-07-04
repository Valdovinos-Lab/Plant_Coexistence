% readtable
T=readtable('Coex_AllPsp_muAP2_1200m_v2.txt');

% Generate some sample data
x = T.abund1;
y = T.pred_abund1_noNaN;

% Create a scatter plot with black dots
figure;
scatter(x, y, 'filled', 'k');

% Add regression line in dark orange
hold on;
p = polyfit(x, y, 1);
f = polyval(p, x);
plot(x, f, 'Color', [0.8, 0.4, 0], 'LineWidth', 2);
hold off;

% Set background color to white
set(gca, 'Color', 'w');

% Add labels and title with specified font settings
xlabel('Plant abundance', 'FontName', 'Segoe UI', 'FontSize', 20);
ylabel('Predicted plant abundance', 'FontName', 'Segoe UI', 'FontSize', 20);
%title('Plant abundance at equilibrium vs Predicted abundance without AF', 'FontName', 'Segoe UI', 'FontSize', 16);

% Set font size for tick labels
set(gca, 'FontSize', 20);

% Round values before calculating R-squared
x_rounded = round(x, 10); % Round to 10 decimal places
y_rounded = round(y, 10); % Round to 10 decimal places

% Calculate R-squared value
R_squared = corr(x_rounded, y_rounded)^2;

% Display regression line statistics
%text(0.1, 0.9, sprintf('y = %.2fx + %.2f\nR^2 = %.2f', p(1), p(2), R_squared), 'Units', 'normalized', 'FontName', 'Segoe UI', 'FontSize', 14);

% Adjust figure properties
set(gcf, 'Color', 'w');

% Set the figure size and aspect ratio to make it square
fig = gcf; % Get the current figure handle
%fig.Position = [0, 0, 600, 600]; % Set the figure size (width x height) in pixels
fig.Position = [0, 0, 300, 300]; % Set the figure size (width x height) in pixels
axis square; % Make the axis aspect ratio 1:1

% Save the figure in high resolution as a PDF file
saveas(gcf, 'pred_abund_vs_abundance_withoutAF_small.emf', 'emf'); % Save the figure as a PDF file

