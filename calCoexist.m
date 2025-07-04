function [k_i, rho_ik, pred_p_i, persistedP_i, criterion_check_i] = calCoexist(Alpha, p, a, metadata)
% This function calculates coexistence metrics for plant-pollinator networks, including:
% - Fitness (k_i)
% - Niche overlap (rho_ik) between each plant species pair in the network
%   -> symmetric matrix: the niche overlap between species i and k 
%      is the same asthat of plant species k and i.
% - Predicted abundances (p_i)
% - Recruitment potential discounting for inter-specific competition (U_i)
% - Visit quality metrics (wa_sigma)
% - Pollination events (Q_i)
% - Recruitment success (Gamma)
% - Seed production (S_i)
% - Visits per plant (sVisits_perP_i)
% - Persistence indicators (persistedP_i)
% - Coexistence criterion check (criterion_check_i)

% Extract parameters from metadata
plant_qty = metadata.plant_qty;
e = full(metadata.e);
mu_p = metadata.mu_p;
u = metadata.u;
w = metadata.w;
g = metadata.g;
tau = metadata.tau;
epsilon = metadata.epsilon;

%% Calculate basic metrics
% Fitness calculation
k_i = g .* sum(e .* Alpha * diag(a .* tau.^2), 2) - mu_p;

% Competition and visits
U_i = (1 - u'*p + u.*p);  % Recruitment potential discounting for inter-specific competition
Visits = diag(p) * Alpha * diag(a.*tau);
Visits_perP = Alpha * diag(a.*tau);
sVisits_perP_i = sum(Visits_perP, 2);

% Quality metrics
sigma = diag(p.*epsilon) * Alpha;
sigma = sigma * diag(1./(sum(sigma)+realmin));
wa_sigma = sum(sigma.*Alpha.*Visits,2)./sum(Alpha.*Visits,2);

% Predicted abundances
pred_p_i = (U_i)./w - mu_p./(g.*w.*e.*wa_sigma.*sVisits_perP_i);
Q_i = wa_sigma.*sVisits_perP_i;  % Pollination events

% Seed production and recruitment
S_i = sum(e .* sigma .* Visits_perP, 2);
Gamma = g .* (1 - u'*p - w.*p + u.*p);

%% Calculate niche overlap (rho)
% Initialize matrices for direct effects
Aii_1stTerm = -w .* g .* S_i;
Aii_2ndTerm = (Gamma./ p) .* sum(e .* Visits_perP .* sigma .* (1-sigma), 2);

% Initialize matrices for inter-specific effects
Aik_1stTerm = zeros(plant_qty, plant_qty);
Aik_2ndTerm = zeros(plant_qty, plant_qty);
rho_ik = zeros(plant_qty, plant_qty);

% Calculate inter-specific effects
for k = 1:plant_qty
    for i = 1:plant_qty
        Aik_1stTerm(i,k) = -g(i) * u(k) * S_i(i);
        Aik_2ndTerm(i,k) = -Gamma(i)/p(k) * sum(e(i,:) .* Visits_perP(i,:) .* sigma(i,:) .* sigma(k,:));
        
        % Calculate rho components
        rho_num = Aik_1stTerm(i,k) * Aik_1stTerm(k,i) + ...
                  Aik_1stTerm(k,i) * Aik_2ndTerm(i,k) + ...
                  Aik_1stTerm(i,k) * Aik_2ndTerm(k,i) + ...
                  Aik_2ndTerm(i,k) * Aik_2ndTerm(k,i);
        
        rho_den = Aii_1stTerm(i) * Aii_1stTerm(k) - ...
                  Aii_1stTerm(i) * Aii_2ndTerm(k) - ...
                  Aii_1stTerm(k) * Aii_2ndTerm(i) + ...
                  Aii_2ndTerm(k) * Aii_2ndTerm(i);
        
        rho_ik(i,k) = sqrt(rho_num/rho_den);
    end
end

% Set diagonal elements to NaN
rho_ik(eye(size(rho_ik))==1) = nan;

% Calculate fitness differences matrix efficiently
k_i2 = k_i + realmin; % Add small constant to avoid division by zero
fitness_diff_ik = k_i2 .\ k_i2'; % Outer division creates the full matrix

% Calculate invasibility criterion: rho < k_k/k_i < 1/rho
one_over_rho_ik = 1./rho_ik;
Mutual_Inv = (rho_ik < fitness_diff_ik) & (fitness_diff_ik < one_over_rho_ik);
Mutual_Inv = double(Mutual_Inv); % converting logical to double
Mutual_Inv(eye(size(Mutual_Inv))==1) = nan;

% Calculate persistence and criterion check
persistedP_i = double(p > 0.01); % Plants below this value are considered extinct
nMutInvaded_i = sum(Mutual_Inv, 'omitnan')';
nPersisted = sum(persistedP_i);
criterion_check_i = double((nPersisted*persistedP_i)==(nMutInvaded_i+persistedP_i));


end
