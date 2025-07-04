
% Create an empty table with variable names and data types
T = table('Size', [0, 6], 'VariableTypes', {'string', 'double', 'double', 'double', 'double', 'double'}, ...
    'VariableNames', {'File_Name', 'P', 'A', 'S', 'L', 'C'});

% Load data and populate the table
d = dir('M_PL*');
metrics_empirical = zeros(length(d), 5);

% Initialize a cell array to store all empirical matrices
M_empirical = cell(1, length(d));

for i = 1:length(d)
    M = load(d(i).name);
    M(M>0)=1;% making M binary (there are a few networks with number of visits)
    [P, A] = size(M);
    metrics_empirical(i, 1) = P;
    metrics_empirical(i, 2) = A;
    metrics_empirical(i, 3) = P + A;
    L = sum(sum(M)); % total number of links
    metrics_empirical(i, 4) = L;
    metrics_empirical(i, 5) = L / (P * A);
    
    % Add the file name to the first column
    T.File_Name(end+1) = d(i).name;
    
    % Add matrices to the cell
    M_empirical{i} = M;
end

% Populate the table with the calculated metrics
T.P = metrics_empirical(:, 1);
T.A = metrics_empirical(:, 2);
T.S = metrics_empirical(:, 3);
T.L = metrics_empirical(:, 4);
T.C = metrics_empirical(:, 5);

% Write the populated table
writetable(T, 'network_properties_empirical.csv');

% Save all empirical matrices in one .mat file
save('all_empirical.mat', 'M_empirical');

