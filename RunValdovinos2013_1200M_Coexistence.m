% Last Modification 07/06/2024, Davis
% Added effective degree, weigthed average sigma, and an improve calculation
% of sigma_c that includes inter-specific competition for recruitment among
% plant species. Also, deleted hand pollination as it is not providing
% added understanding. Finally, I deleted all about inter-specific
% effects.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last Modification 02/02/2024, Davis
% Added species persistence, changed frCoexist for number of mutually
% invaded plants, and added check on whether Chesson criterion worked.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last Modification 07/08/2023, Davis
% Also calls IntegrateValdovinos2013_dataset_handPollination.m, which calls
% Valdovinos2013_rhs_handPollination.m which has modifications needed
% for simulating hand pollination, that is, it assumes:
% 1. Sigma and tau are equal to 1.
% 2. A slightly lower abundance of pollinator vectors (humans)
%    than the average abundance of pollinators in bening environmental
%    conditions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last Modification 06/25/2023, Davis
% Saves network_metadata for each run as it is hard to retrieve correctly
% later when running other analysis such as species coexistence or
% efficiencies.
% Also runs species coexistence metrics to make sure they are the correct
% ones, matching the correct randomly generated parameters saved in
% network_metadata.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last Modification 08/03/2019, Ann Arbor
% Cleaning up my codes
% Only run the dynamics, without species invasions or removals
% Runs the dynamics for 1200 dataset
% RunValdovinos2013_1200M.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
global J_pattern network_metadata

dataset=1200;
muAP=2;
sem=0;

load(sprintf('%dm.mat',dataset)) % Already sorted by degree
                    %to get matrix i: cell2mat(m1200(i))
registro=load(sprintf('registro%dm4Matlab.csv',dataset));             
registro=[registro registro(:,end)>0.2]; % Notting the networks whose NODFst > 0.2 as nested (1)
                                         % and the ones whose NODFst < 0.2 as non-nested (1)
S=registro(:,4);
S_category=[S>10 S>60 S>150];
S_category=sum(S_category,2);

registro=[registro(:,1:4) S_category registro(:,5:end)];

A_dataset=[];
Alphas=cell(dataset,2);
N_metadata=cell(dataset,2);
CoexistSummary_AllSp_dataset=[];
Reg_i_Res_Inter_dataset=[];
CoexistRes_mean_dataset=[];
CoexistRes_ci_dataset=[];
registro_dataset=[];
Reg_Allsp_i_dataset=[];
Reg_Allsp_j_dataset=[];

for i=1%:dataset
    rand('seed',sem+i);
    In=cell2mat(m1200(i));% change for every dataset!!!!
    %In=cell2mat(m800(i));% change for every dataset!!!!
    [plant_qty, animal_qty]= size(In);
    J_pattern = J_zero_pattern(In) ;
    
    P01=zeros(plant_qty,2);
    R01=zeros(plant_qty,2);
    A01=zeros(animal_qty,2);
    
    %[k_i, D_i, w_degree, degree, sigma_ci, wa_sigma, p_i,
    ki_01=zeros(plant_qty,2);
    Di_01=zeros(plant_qty,2);
    wdegree_01=zeros(plant_qty,2);
    sigmac_01=zeros(plant_qty,2);
    wasigma_01=zeros(plant_qty,2); % weighted average quality of visits received by i
    pi_01=zeros(plant_qty,2); % predicted plant abundances at equilibrium
    
    Gamma_01=zeros(plant_qty,2);
    mort_recruit_01=zeros(plant_qty,2);
    
    PollenDepi_01=zeros(plant_qty,2);% Per-capita total pollen deposition on stigmas of i
    Si_01=zeros(plant_qty,2);        % Per-capita total production of seeds i
    
    PersPi_01=zeros(plant_qty,2);
   
    fitness_diff_i_01=zeros(plant_qty,2);
    rho_k_i_01=zeros(plant_qty,2);
    Mutual_Invi_01=zeros(plant_qty,2);
    
    criterion_check_i_01=zeros(plant_qty,2);
    
    sVisits_perP_i_01=zeros(plant_qty,2); % Per-capita total visits to i
    
    Omega_ik_i_01=zeros(plant_qty,2);
    S_ik_01=zeros(plant_qty^2,2);
    Omega_ik_01=zeros(plant_qty^2,2);
    
    rho_ik_01=zeros(plant_qty^2,2);
    fitness_diff_ik_01=zeros(plant_qty^2,2);
    matMutual_Inv_01=zeros(plant_qty^2,2);
    
                   
    %% Without and with AF
    for frG=[0,1]
            
        vectG=frG*ones(1,animal_qty);
    
        [plantsf, nectarf, animalsf, alphasf, network_metadata]=IntegrateValdovinos2013_dataset(vectG,In,muAP);
        
        [k_i, D_i, w_degree, degree, sigma_ci, wa_sigma, p_i, Gamma, S_i, mort_recruit, PollenDep_i, sVisits_perP_i, persistedP_i,...
        S_ik, Omega_ik, rho_ik, fitness_diff_ik, Mutual_Inv, criterion_check_i]= calCoexist(alphasf,plantsf,animalsf,network_metadata);
        
        ki_01(:,frG+1)=k_i;
        Di_01(:,frG+1)=D_i;
        wdegree_01(:,frG+1)=w_degree;
        sigmac_01(:,frG+1)=sigma_ci;
        wasigma_01(:,frG+1)=wa_sigma; % weighted average quality of visits received by i
        pi_01(:,frG+1)=p_i;
        
        Gamma_01(:,frG+1)=Gamma;
        mort_recruit_01(:,frG+1)= mort_recruit;
        
        PollenDepi_01(:,frG+1)=PollenDep_i;
        Si_01(:,frG+1)=S_i;
        sVisits_perP_i_01(:,frG+1)=sVisits_perP_i;
                
        Omega_ik_i_01(:,frG+1)=nanmean(Omega_ik)'; % Omega is symmetrical, that is: Omega_ik = Omega_ki 
        S_ik_01(:,frG+1)= reshape(S_ik,[],1);
        Omega_ik_01(:,frG+1)= reshape(Omega_ik,[],1);
                   
        PersPi_01(:,frG+1)=persistedP_i; % Species persisted: yes = 1, no = 0.
        fitness_diff_i_01(:,frG+1)=nanmean(fitness_diff_ik)';% Average fitness diff of each plant sp with each of the other plant sp. 
                                                          % The mean should be taken across rows for each column because
                                                          % the differences are calculated k_k/k_i, where k is in the columns
                                                          % and i is in the rows.
        rho_k_i_01(:,frG+1)=nanmean(rho_ik)'; % Also symmetrical, that is: rho_ik = rho_ki 
        Mutual_Invi_01(:,frG+1)=nansum(Mutual_Inv)'; % Also symmetrical in terms of who mutually invades whom, denoted by a 1.
                                                     % Thus, the sum indicates how many species coexisted with the focal sp.
        rho_ik_01(:,frG+1)= reshape(rho_ik,[],1);
        fitness_diff_ik_01(:,frG+1)=reshape(fitness_diff_ik,[],1);
        matMutual_Inv_01(:,frG+1)=reshape(Mutual_Inv,[],1);
        criterion_check_i_01(:,frG+1)=criterion_check_i;        
        
        P01(:,frG+1)=plantsf;
        Alphas{i,frG+1}=alphasf;
        R01(:,frG+1)=nectarf;
        A01(:,frG+1)=animalsf;
        
    end
 
    degree_Pi=sum(In,2);
    Plant_ID=(1:plant_qty)';
    registro_i= registro(i,:);
    
    registro_i_useful=registro(i,[1:6,19,20]);
           
    Reg_Allsp_i= ones( plant_qty,length(registro_i_useful)) * diag(registro_i_useful);
    Reg_Allsp_j= ones( animal_qty,length(registro_i_useful)) * diag(registro_i_useful);
    
    Reg_Allsp_i_dataset=[Reg_Allsp_i_dataset; Reg_Allsp_i]; % For all plants
    Reg_Allsp_j_dataset=[Reg_Allsp_j_dataset; Reg_Allsp_j]; % For all pollinators
    
    P_dataset_i=registro_i_useful(:,2);
    
    fPersP_01 = sum(PersPi_01)/plant_qty;
    fPersPi_01 = fPersP_01 .* ones(plant_qty,2);% Species persistence repeated i times
                                                % only for storing in the spreadsheet for AllP
    
 %% Data type 1: Saving All sp, per each network, concatenated across networks in the loop 
    CoexistSummary_AllSp= [Plant_ID degree_Pi ki_01 Di_01 wdegree_01 pi_01 P01 R01 wasigma_01 sigmac_01...
        Gamma_01 sVisits_perP_i_01 PollenDepi_01 Si_01 mort_recruit_01 PersPi_01...
        Omega_ik_i_01 rho_k_i_01 fitness_diff_i_01 Mutual_Invi_01 fPersPi_01 criterion_check_i_01];
    
    CoexistSummary_AllSp_dataset=[CoexistSummary_AllSp_dataset; CoexistSummary_AllSp];
    A_dataset=[A_dataset; A01];
       
        
 %% Data type 3: Saving one summary descriptor (sum or mean across species) per metric per network
 
    CoexistSummary=[degree_Pi ki_01 Di_01 wdegree_01 pi_01 P01 R01 wasigma_01 sigmac_01...
        Gamma_01 sVisits_perP_i_01 PollenDepi_01 Si_01 mort_recruit_01];
    CoexistSummary(CoexistSummary==0)=NaN;
    CoexistSummary_mean=nanmean(CoexistSummary);
    CoexistSummary_ci=interv_conf(CoexistSummary,1);
    
    Inter = [Omega_ik_01 rho_ik_01 fitness_diff_ik_01];
    Inter(Inter==0)=NaN;
    CoexistInter_mean=nanmean(Inter);
    CoexistInter_ci=interv_conf(Inter,1);
    
    criterion_check_01 = sum(criterion_check_i_01)./plant_qty;
    
    CoexistRes_mean_dataset=[CoexistRes_mean_dataset; CoexistSummary_mean mean(A01) CoexistInter_mean...
        fPersP_01 criterion_check_01];
    CoexistRes_ci_dataset=[CoexistRes_ci_dataset; CoexistSummary_ci CoexistInter_ci];
    
    registro_dataset=[registro_dataset;registro_i_useful]; 

    
end

CoexistRes_mean=[registro_dataset CoexistRes_mean_dataset];

CoexistSummary_AllSp_dataset2=[Reg_Allsp_i_dataset CoexistSummary_AllSp_dataset];

%% Data type 1 -> Writing table for all plant species across all networks
T_CoexistSummary_AllSp_dataset= array2table(CoexistSummary_AllSp_dataset2,'VariableNames',{'matID','P','A','S','S_category','C','NODFst','N_lowHigh',...
    'Plant_ID','degree','fitness0','fitness1','Effective_deg0','Effective_deg1','wdegree0','wdegree1','predicted_abund0','predicted_abund1',...
    'abund0','abund1','R0','R1','wasigma0','wasigma1','sigma_c0','sigma_c1','Gamma0','Gamma1','sVisits_perP0','sVisits_perP1',...
    'PollenDep_perP0','PollenDep_perP1','SeedProd_perP0','SeedProd_perP1','mort_recruit0','mort_recruit1','PersYoN0','PersYoN1'...
    'Omega0','Omega1','rho0','rho1','fitness_diff0','fitness_diff1','nMutInv0','nMutInv1','fPersP0', 'fPersP1', 'criterion_check0', 'criterion_check1'});

%writetable(T_CoexistSummary_AllSp_dataset,sprintf('Coex_AllPsp_muAP%d_%dm.txt',muAP,dataset),'Delimiter',' ')

%% Data type 2 -> Writing table of all summary stats of each network

T_CoexistRes_mean= array2table(CoexistRes_mean,'VariableNames',{'matID','P','A','S','S_category','C','NODFst','N_lowHigh',...
    'mdegree','mfitness0','mfitness1','mEffective_deg0','mEffective_deg1','wdegree0','wdegree1','mpredicted_abund0','mpredicted_abund1',...
    'mabund0','mabund1','mR0','mR1','mwasigma0','mwasigma1','msigma_c0','msigma_c1','mGamma0','mGamma1','sVisits_perP0','sVisits_perP1',...
    'mPollenDep_perP0','mPollenDep_perP1','mSeedProd_perP0','mSeedProd_perP1','mmort_recruit0','mmort_recruit1','PersYoN0','PersYoN1'...
    'mOmega0','mOmega1','mrho0','mrho1','mfitness_diff0','mfitness_diff1','fPersP0', 'fPersP1', 'criterion_check0', 'criterion_check1'});

%writetable(T_CoexistRes_mean,sprintf('Coex_mean_muAP%d_%dm.txt',muAP,dataset),'Delimiter',' ')

toc
