function [k_i, p_i, U_i, wa_sigma, sigma_ci, Q_i, Q_ci, Gamma, S_i, sVisits_perP_i, persistedP_i, criterion_check_i, D_j, wd_i]= calCoexist(Alpha,p,a,metadata)

% Called by RunValdovinos2013_1200M.m 
% Calculates the partial derivatives (direct effect) of plant species k on
% plant species i, as well as the effect of plant species i on itself.
% Finally, this function calculates niche
% differentation (rho, sensu Chesson 1989), niche overlap (1-rho), and
% fitness diferrences.
% Mutual Invasibility criterion: rho < k_i/k_k < 1/rho ; ensures
% coexistence.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Last Modification 07/014/2024, Davis
% Removed effective degree.

%% Last Modification 07/06/2024, Davis
% Added effective degree, weigthed average sigma, and an improve calculation
% of sigma_c that includes inter-specific competition for recruitment among
% plant species.

%% Last modification 3/27/2024 in Davis, CA
% Added calculations that determine whether a plant population has a
% negative population growth which starts the non-reversible process to
% extinction (at least in networks with AF), which are:
% Per-capita mortality rate of a plant species over the fraction of its
% seeds produced that become adults (mort_recruit, below).
% Per-capita seed production must be equal or higher than that quantity for
% the plant population to keep or increase its density.
% => this was not very useful...

%% Last modification 3/26/2024 in Davis, CA
% Added calculations from Valdovinos & Marsland 2020. Specifically those
% that determine when a plant species go extinct:
% sij: # of seedlings produced per plant lifetime under optiomal conditions
% -> no other plant species nearby to contaminate the pollen and the field
% is clear of all competing plants.
% di: Fraction of total rewards lost to dilution when plants are at
% carrying capacity 1/wi.
% sigma_c: critical value of quality of visits under which a plant goes
% extinct under ideal conditions, which value is approximately 4/sij.

%% Last modification 5/29/2023 in Davis, CA
% Modified from calValMechs.m to calculate all the components of the
% partial derivatives indicating the direct effect of plant species k on
% plant species i, as well as the effect of plant species i on itself. Then
% those derivatives are used to identify the intra- and inter-specific
% direct effects caused by the plant species on each other and itsels via
% recruitment or via pollination. Finally, this function calculates niche
% differentation (rho, sensu Chesson 1989), niche overlap (1-rho), and
% fitness diferrences.

%% Calling all parameters stored in structure network_metadata 
% To run bits of this code after running RunValdovinos2013_singleM_Coexist.m 
% use: Alpha=alphasf; p=plantsf; a=animalsf;

plant_qty  = metadata.plant_qty ;
e       = full(metadata.e) ;
mu_p    = metadata.mu_p ;

u       = metadata.u ;
w       = metadata.w ;

g       = metadata.g ;
tau     = metadata.tau ;

epsilon = metadata.epsilon ; % This parameter has value 1, but when there is added variance to it (varp>0),
                                     % some plant species get below 1 others above 1, and it then affects
                                     % a bit the dynamics and the calculated metrics. I'm keeping because 
                                     % the equations have it so I can trakc the effects when its varp > 0

%% Demographic potential, weighted-average and critical sigma, and predicted plant abundances at equilibrium

k_i = g .* sum(e .* Alpha * diag(a .* tau.^2),2) - mu_p; % Fitness of each plant species assuming 
                                                         % no effect of any type of competition .
                                                         
% Effective degree of pollinators and weigthed degree of plants
[D_j, wd_i]=effective_degree(Alpha);

U_i=(1 - u'*p + u.*p) ;% Total interspecific competition.
                       % The added u.*p is substracting the corresponding
                       % u to self (i,i) since it's inter-specific competition
                                             
Visits = diag(p)* Alpha * diag(a.*tau);
Visits_perP = Alpha * diag(a .*tau) ; % Per-plant (per-capita) visits of pollinator sp j
                                      % assigned to each plant sp i. It's a
                                      % matrix (i,j).
sVisits_perP_i = sum(Visits_perP,2) ;
                                      
sigma = diag(p.*epsilon) * Alpha ; % generally, epsilon is 1 and its variance is 0, so it's like it's not there.
                                   % But there might be cases in which want to turn epsilon on.
                                   
sigma = sigma * diag(1./(sum(sigma)+realmin)) ; % Quality of visits. Matrix (i,j)

wa_sigma = sum(sigma.*Alpha.*Visits,2)./sum(Alpha.*Visits,2);

p_i=(U_i)./w - mu_p./(g.*w.*e.*wa_sigma.*sVisits_perP_i);

sigma_ci=mu_p./(g.*e.*sVisits_perP_i.*U_i);

Q_i=wa_sigma.*sVisits_perP_i;
Q_ci = mu_p./(g.*e.*U_i);

%% Starting calculations for niche overlap sensu Chesson

S_i = sum(e .* sigma .* Visits_perP, 2); % Per-capita total production of seeds i

Gamma = g .* (1 - u'*p - w.*p + u.*p) ;% Fraction of seeds that recruit to adult plants.
                                       % The added u.*p is substracting the corresponding
                                       % u to self (i,i) since it is
                                       % inter-specific competition

s_ij = e .* Visits_perP ; % per-plant production of seeds resulting from visits of j
                          % assuming that the pollinator is a specialist on i
                          % without sigma because sigma is accounted later
                          % in the calculations of Aii_2ndTerm as F_ij.

Visits_j = sum(Visits) ; % Total visits performed by pollinator sp j

F_ij = Visits*diag(1./Visits_j) ; % Fraction of j's visits assiged to i (same as sigma!)

%% Calculating intra-specific direct effects (diagonal of the Jacobian)

Aii_1stTerm = - w.* g.* S_i ; % Direct effect of plant species i on itself via 
                        % intra-specific competition for seed recruitment.

Aii_2ndTerm = (Gamma./ p) .* sum(s_ij.* F_ij.*(1-F_ij), 2) ; % Direct effect of 
                                % plant species i on itself via pollination

%Aii = Aii_1stTerm + Aii_2ndTerm ; % Total direct effect of plant species i on itself via 
                      % both seed recruitment and pollination.
                      
%% Allocating zero matrices for calculating inter-specific direct effects

fitness_diff_ik = zeros(plant_qty,plant_qty);

Aik_1stTerm = zeros(plant_qty, plant_qty);
Aik_2ndTerm = zeros(plant_qty, plant_qty);
S_ik = zeros(plant_qty, plant_qty);
Omegaj_ik = cell(plant_qty,plant_qty);
Omega_ik = zeros(plant_qty,plant_qty);

for k=1:plant_qty  % the effecter in columns
    
    for i=1:plant_qty % the receiver of the effect in the rows
        
        fitness_diff_ik(i,k) = k_i(k)/k_i(i); % fitness differences between the effecter and the receiver. Effecter is the focal species
                                              
        Aik_1stTerm(i,k) = - g(i)* u(k)* S_i(i); % all three terms are vectors with as many elements
                                           % as plants from which only one plant is taken
        Aik_2ndTerm(i,k) = - Gamma(i)/ p(k) * sum(s_ij(i,:).* F_ij(i,:).*F_ij(k,:)); % Each time a species goes extinct, i.e. p=0,
                                                                                     % Gamma/ p(k) becomes Inf, and when the overlap 
                                                                                     % Fij* Fkj is zero, then the entire term becomes NaN
                                                                                     % because Inf * 0 = NaN.
        
        S_ik(i,k) = sum(s_ij(i,:).* F_ij(i,:).*F_ij(k,:));
        Omegaj_ik{i,k} = F_ij(i,:).*F_ij(k,:); % sigma_ij * sigma_kj - Overlap for pollinator j's visits
        Omega_ik(i,k) = sum(F_ij(i,:).*F_ij(k,:),2); % summed overlap across all pollinators
                          
    end
end
        
%Aik = Aik_1stTerm + Aik_2ndTerm;

%rho_num_ik = zeros(plant_qty,plant_qty);
%rho_den_ik = zeros(plant_qty,plant_qty);
rho_ik = zeros(plant_qty,plant_qty);

for k=1:plant_qty
    
    for i=1:plant_qty
        
        rho_num1 =   Aik_1stTerm(i,k) * Aik_1stTerm(k,i);
        rho_num2 =   Aik_1stTerm(k,i) * Aik_2ndTerm(i,k);
        rho_num3 =   Aik_1stTerm(i,k) * Aik_2ndTerm(k,i);
        rho_num4 =   Aik_2ndTerm(i,k) * Aik_2ndTerm(k,i);
        
        rho_den1 =   Aii_1stTerm(i) * Aii_1stTerm(k);
        rho_den2 = - Aii_1stTerm(i) * Aii_2ndTerm(k);
        rho_den3 = - Aii_1stTerm(k) * Aii_2ndTerm(i);
        rho_den4 =   Aii_2ndTerm(k) * Aii_2ndTerm(i);
        
        rho_num = rho_num1 + rho_num2 + rho_num3 + rho_num4;
        rho_den = rho_den1 + rho_den2 + rho_den3 + rho_den4;
        
        %rho_num_ik(i,k) = rho_num;
        %rho_den_ik(i,k) = rho_den;
        rho_ik(i,k) = sqrt( rho_num/rho_den );
    
    end
    
end

one_over_rho_ik = 1./rho_ik;

% Invasibility criterion: rho < k_i/k_k < 1/rho
Mutual_Inv=(rho_ik < fitness_diff_ik)&(fitness_diff_ik < one_over_rho_ik);
Mutual_Inv=double(Mutual_Inv); % converting logical to double

% Adding NaN to the diagonal of all matrices
Mutual_Inv(eye(size(Mutual_Inv))==1) = nan ; % Making the diagonal (invasibility of i on i) NAN
%fitness_diff_ik(eye(size(fitness_diff_ik))==1) = nan ; % Making the diagonal (difference between i and i) NAN
%S_ik(eye(size(S_ik))==1) = nan ; % Making the diagonal (effect of i on i) NAN
%Omega_ik(eye(size(Omega_ik))==1) = nan ; % Making the diagonal (overlap between i and i) NAN
%rho_ik(eye(size(rho_ik))==1) = nan ; % Making the diagonal (overlap between i and i) NAN

% Checking if Chesson's criterion worked
persistedP_i=double(p > 0);
nMutInvaded_i = sum(Mutual_Inv, 'omitnan')';
nPersisted=sum(persistedP_i);
criterion_check_i=double((nPersisted*persistedP_i)==(nMutInvaded_i+persistedP_i));


end