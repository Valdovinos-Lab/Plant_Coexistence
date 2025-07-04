% This code accompanies the manuscript "Equalizing Effect of Pollinator Adaptive 
% Foraging on Plant Network Coexistence" by Fernanda S. Valdovinos, Taranjot Kaur, 
% and Robert Marsland III. Developed by Fernanda S. Valdovinos.

dataset = 1200;
muAP = 2;

% Load network data and properties
load(sprintf('%dm.mat',dataset)) % Networks sorted by degree
T = readtable('network_properties_1200m.csv');

% Classify networks by nestedness and size
N_lowHigh = T.NODFst > 0.2; % 1=nested (NODFst>0.2), 0=non-nested
S = T.S;
S_category = [S>10 S>60 S>150];
S_category = sum(S_category,2);

% Store network properties
Network_properties = [T.matrix T.P T.A T.S S_category T.C T.NODFst N_lowHigh T.nodfc];

% Initialize results storage
CoexistSummary_AllSp_dataset = [];
CoexistRes_mean_dataset = [];
NetProp_Allsp_i_dataset = [];

% Main loop over all networks
for i = 1%:dataset
    In = cell2mat(m1200(i));
    [plant_qty, animal_qty] = size(In);
    
    % Initialize arrays for storing results
    rho_i01 = zeros(plant_qty,2);  % One row per plant, two columns for frG=0,1
    k_i01 = zeros(plant_qty,2);
    P01 = zeros(plant_qty,2);
    R01 = zeros(plant_qty,2);
    pred_p_i01 = zeros(plant_qty,2);
    persistedP_i01 = zeros(plant_qty,2);
    criterion_check_i01 = zeros(plant_qty,2);
    degree_Pi = sum(In,2);
    
    % Arrays to store results for both foraging conditions
    fPersP_01 = zeros(1,2);
    rho_mean01 = zeros(1,2);
    criterion_check_01 = zeros(1,2);
    
    % Run simulations for fixed (frG=0) and adaptive (frG=1) foraging
    for frG = [0,1]
        vectG = frG*ones(1,animal_qty);
        metadata = create_metadata(vectG,In,muAP,i);
        
        % Run dynamics simulation
        [plantsf, nectarf, animalsf, alphasf] = IntegrateValdovinos2013_dataset(metadata);
        
        % Calculate coexistence metrics
        %Alpha=alphasf; p=plantsf; a=animalsf;
        [k_i, rho_ik, pred_p_i, persistedP_i, criterion_check_i] = ...
            calCoexist(alphasf,plantsf,animalsf,metadata);

        % Taking the mean across the rho vealues for each plant
        rho_i01(:,frG+1) = mean(rho_ik, 'omitnan')';

        % Mean rho for the network (upper triangle given that rho_ik is symmetric)
        % That is, the niche overlap between plant species i and k is the same as
        % that of plant species k and i.
        upper_triangle = triu(rho_ik, 1); % Get upper triangle (1 means exclude diagonal)
        rho_mean01(1,frG+1) = mean(upper_triangle(upper_triangle ~= 0)); % Mean of non-zero elements
        
        % Store results
        k_i01(:,frG+1) = k_i;
        pred_p_i01(:,frG+1) = pred_p_i;
        persistedP_i01(:,frG+1) = persistedP_i;
        criterion_check_i01(:,frG+1) = criterion_check_i;
        P01(:,frG+1) = plantsf;
        R01(:,frG+1) = nectarf;
        
        % Calculate and store network-level metrics
        fPersP_01(frG+1) = sum(persistedP_i)/plant_qty;
        criterion_check_01(frG+1) = sum(criterion_check_i)/plant_qty;
    end
    
    % Repeat network properties for plant species in the same network
    Network_properties_i = Network_properties(i,:);
    NetProp_Allsp_i = ones(plant_qty,length(Network_properties_i)) * diag(Network_properties_i);
    NetProp_Allsp_i_dataset = [NetProp_Allsp_i_dataset; NetProp_Allsp_i];
    
    % Store results for all plants species
    Plant_ID = (1:plant_qty)';
    CoexistSummary_AllSp = [Plant_ID degree_Pi k_i01 rho_i01 pred_p_i01 P01 R01 ...
        persistedP_i01 criterion_check_i01];
    CoexistSummary_AllSp_dataset = [CoexistSummary_AllSp_dataset; CoexistSummary_AllSp];
    
    % Store network-level results
    network_results = [rho_mean01 fPersP_01 criterion_check_01];  % [fPersP0 fPersP1 criterion0 criterion1]
    CoexistRes_mean_dataset = [CoexistRes_mean_dataset; network_results];
end

% Create final results tables
CoexistRes_mean = [Network_properties CoexistRes_mean_dataset];
T_CoexistRes_mean = array2table(CoexistRes_mean,'VariableNames',...
    {'matID','P','A','S','S_category','C','NODFst','N_lowHigh','NODFc',...
    'rho_mean0','rho_mean1','fPersP0','fPersP1','criterion_check0', 'criterion_check1'});

CoexistSummary_AllSp_dataset2 = [NetProp_Allsp_i_dataset CoexistSummary_AllSp_dataset];
T_CoexistSummary_AllSp_dataset = array2table(CoexistSummary_AllSp_dataset2,'VariableNames',...
    {'matID','P','A','S','S_category','C','NODFst','N_lowHigh','NODFc',...
    'Plant_ID','degree','k0','k1','mean_rho0','mean_rho1','pred_p0','pred_p1','p0','p1','R0','R1',...
    'PersYoN0','PersYoN1','criterion_check0','criterion_check1'});

% Save results
writetable(T_CoexistRes_mean,sprintf('Coex_mean_muAP%d_%dm.csv',muAP,dataset),'Delimiter',' ')
writetable(T_CoexistSummary_AllSp_dataset,sprintf('Coex_AllPsp_muAP%d_%dm.csv',muAP,dataset),'Delimiter',' ')