% Main code for simulations by Valdovinos et al 2025 on plant species coexistence
% It runs the Valdovinos et al's 2013 model for 1200 algorithm-generated
% networks, without (frG=0) and with (frG=1) adaptive foraging; by calling  
tic

dataset=1200;
muAP=2;

load(sprintf('%dm.mat',dataset)) % Already sorted by degree
                    %to get matrix i: cell2mat(m1200(i))
                    
T=readtable('network_properties_1200m.csv');
N_lowHigh=T.NODFst>0.2;% Notting the networks whose NODFst > 0.2 as nested (1)
                       % and the ones whose NODFst < 0.2 as non-nested (0)
S=T.S;
S_category=[S>10 S>60 S>150];
S_category=sum(S_category,2);

Network_properties=[T.matrix T.P T.A T.S S_category T.C T.NODFst N_lowHigh T.nodfc] ;

A_dataset=[];
Alphas=cell(dataset,2);
N_metadata=cell(dataset,2);
CoexistSummary_AllSp_dataset=[];
Reg_i_Res_Inter_dataset=[];
CoexistRes_mean_dataset=[];
CoexistRes_ci_dataset=[];
registro_dataset=[];
NetProp_Allsp_i_dataset=[];
Reg_Allsp_j_dataset=[];
didjEff_dataset=[];
parsP_dataset=[];

for i=1:dataset

    In=cell2mat(m1200(i));
    [plant_qty, animal_qty]= size(In);
    
    P01=zeros(plant_qty,2);
    R01=zeros(plant_qty,2);
    A01=zeros(animal_qty,2);
    Dj01=zeros(animal_qty,2);
    
    ki_01=zeros(plant_qty,2);
    wdi_01=zeros(plant_qty,2);
    sigma_ci01=zeros(plant_qty,2);
    Qi_01=zeros(plant_qty,2);
    Q_ci01=zeros(plant_qty,2);
    pi_01=zeros(plant_qty,2); % predicted plant abundances at equilibrium
    sVisits_perP_i_01=zeros(plant_qty,2); % Per-capita total visits to i
    PersPi_01=zeros(plant_qty,2);
    
    Gamma_01=zeros(plant_qty,2);
    Si_01=zeros(plant_qty,2);
    Ui_01=zeros(plant_qty,2);
    wa_sigma_01=zeros(plant_qty,2);
    criterion_check_i_01=zeros(plant_qty,2);
    Uidwi_01=zeros(plant_qty,2);
        
    degree_Pi=sum(In,2);
    degree_Aj=sum(In);
    [rowIdx, colIdx] = ind2sub([plant_qty, animal_qty], (1:numel(In))');
    Eff01=zeros(numel(In), 2 );
    
    parsP_01=zeros(plant_qty,8);
                   
    %% Without and with AF
    for frG=[0,1]
            
        vectG=frG*ones(1,animal_qty);
        metadata = create_metadata(vectG,In,muAP,i);

        [plantsf, nectarf, animalsf, alphasf]=IntegrateValdovinos2013_dataset(metadata);
        
        [k_i, p_i, U_i, wa_sigma, sigma_ci, Q_i, Q_ci, Gamma, S_i, sVisits_perP_i, persistedP_i, criterion_check_i, D_j, wd_i]= calCoexist(alphasf,plantsf,animalsf,metadata);
        
        parsP_01( :, (1:4) + frG*4 )=[metadata.mu_p metadata.g metadata.w metadata.e];

        ki_01(:,frG+1)=k_i;
        pi_01(:,frG+1)=p_i;
        wdi_01(:,frG+1)=wd_i;
        sigma_ci01(:,frG+1)=sigma_ci;
        Qi_01(:,frG+1)=Q_i;
        Q_ci01(:,frG+1)=Q_ci;
        sVisits_perP_i_01(:,frG+1)=sVisits_perP_i;
        PersPi_01(:,frG+1)=persistedP_i; % Species persisted: yes = 1, no = 0.
        
        Gamma_01(:,frG+1)=Gamma;
        Si_01(:,frG+1)=S_i;
        Ui_01(:,frG+1)=U_i;
        Uidwi_01(:,frG+1)=U_i./metadata.w;
        wa_sigma_01(:,frG+1)=wa_sigma;
        criterion_check_i_01(:,frG+1)=criterion_check_i;        
        
        P01(:,frG+1)=plantsf;
        Alphas{i,frG+1}=alphasf;
        R01(:,frG+1)=nectarf;
        A01(:,frG+1)=animalsf;
        Dj01(:,frG+1)=D_j;

        Eff01(:,frG+1)=alphasf(:);
        
    end
   
%% Indicating which species are specialist or generalist using terciles.
    % These terciles indicate the 1/3 most specialist (1) and 1/3 most
    % generalist (3) plant species. IMPORTANT: these terciles rely in networks 
    % already ordered in ascending degree order for plants.
    tercile1=round(plant_qty/3);
    tercile3=round(plant_qty/3);
    tercile2=plant_qty-tercile1-tercile3;
    vtercile=[ones(tercile1,1); repmat(2,tercile2,1); repmat(3,tercile3,1)]; 
    
    specialist=degree_Pi==1;
    
    Network_properties_i= Network_properties(i,:);
           
    NetProp_Allsp_i= ones( plant_qty,length(Network_properties_i)) * diag(Network_properties_i);
    
    NetProp_Allsp_i_dataset=[NetProp_Allsp_i_dataset; NetProp_Allsp_i]; % For all plants
          
    fPersP_01_All = sum(PersPi_01)/plant_qty;
    
 %% Data type 1: Saving All Plant sp, per each network, concatenated across networks in the loop 

    Plant_ID=(1:plant_qty)';

    CoexistSummary_AllSp= [Plant_ID degree_Pi specialist vtercile ki_01 wdi_01 pi_01 P01 R01 Qi_01 Q_ci01...
        Gamma_01 Si_01 Ui_01 Uidwi_01 wa_sigma_01  sigma_ci01 sVisits_perP_i_01 PersPi_01 criterion_check_i_01];
    
    CoexistSummary_AllSp_dataset=[CoexistSummary_AllSp_dataset; CoexistSummary_AllSp];
       
        
 %% Data type 2: Saving one summary descriptor (sum or mean across species) per metric per network
 
    criterion_check_01 = sum(criterion_check_i_01)./plant_qty;
    
    CoexistRes_mean_dataset=[CoexistRes_mean_dataset; fPersP_01_All criterion_check_01];

%% Data type 3: Saving All Animal sp, per each network, concatenated across networks in the loop 
    
    A_dataset=[A_dataset; degree_Aj' A01 Dj01];
    
%% Data type 4: Efforts with the corresponding plant and pollinator degrees
    didjEff01=[degree_Pi(rowIdx) degree_Aj(colIdx)' Eff01];
    posEff=In>0;
    didjEff01=didjEff01(posEff,:);
    didjEff_dataset=[didjEff_dataset; didjEff01];

%% Data type 5: Saving key parameters, concatenated across networks in the loop 
    parsP_dataset=[parsP_dataset; parsP_01];

end

%CoexistRes_mean=[Network_properties(1,:) CoexistRes_mean_dataset];
CoexistRes_mean=[Network_properties CoexistRes_mean_dataset];

CoexistSummary_AllSp_dataset2=[NetProp_Allsp_i_dataset CoexistSummary_AllSp_dataset];

%% Data type 1 -> Writing table for all plant species across all networks
T_CoexistSummary_AllSp_dataset= array2table(CoexistSummary_AllSp_dataset2,'VariableNames',{'matID','P','A','S','S_category','C','NODFst','N_lowHigh','NODFc',...
    'Plant_ID','degree','specialist','vtercile','fitness0','fitness1','wdegree0','wdegree1','predicted_abund0','predicted_abund1','abund0','abund1','R0','R1', 'Q0','Q1','Q_c0','Q_c1',...
    'Gamma0','Gamma1','sSeedProd_perP0','sSeedProd_perP1','Ui0','Ui1','Ui/wi0','Ui/wi1','wa_sigma0','wa_sigma1', 'sigma_c0',' sigma_c1', 'sVisits_perP0','sVisits_perP1',...
    'PersYoN0','PersYoN1','criterion_check0','criterion_check1'});

fQc_0=T_CoexistSummary_AllSp_dataset.Q0<T_CoexistSummary_AllSp_dataset.Q_c0; % we predict that they go extinct
fwasigmaNaN_0=isnan(T_CoexistSummary_AllSp_dataset.wa_sigma0); % is NaN when visits to plant i are zero, which happens when they go extinct
                                                                                                        % wa_sigma has in its calculation visits/visits, which in the case of extinct plants is 0/0=NaN.
fQc_1=T_CoexistSummary_AllSp_dataset.Q1<T_CoexistSummary_AllSp_dataset.Q_c1;
fwasigmaNaN_1=isnan(T_CoexistSummary_AllSp_dataset.wa_sigma1);

% For the derivation based on the steady-state of the rewards dynamics: p* goes to -Inf when wa_sigma = 0, so the other equilibrium p*=0 applies.
T_CoexistSummary_AllSp_dataset.predicted_abund0(fQc_0)=0;
T_CoexistSummary_AllSp_dataset.predicted_abund0(fwasigmaNaN_0)=0;

T_CoexistSummary_AllSp_dataset.predicted_abund1(fQc_1)=0;
T_CoexistSummary_AllSp_dataset.predicted_abund1(fwasigmaNaN_1)=0;

writetable(T_CoexistSummary_AllSp_dataset,sprintf('Coex_AllPsp_muAP%d_%dm.txt',muAP,dataset),'Delimiter',' ')

%% Data type 2 -> Writing table of all summary stats of each network

T_CoexistRes_mean= array2table(CoexistRes_mean,'VariableNames',{'matID','P','A','S','S_category','C','NODFst','N_lowHigh','NODFc',...
       'fPersP0','fPersP1','criterion_check0', 'criterion_check1'});

writetable(T_CoexistRes_mean,sprintf('Coex_mean_muAP%d_%dm.txt',muAP,dataset),'Delimiter',' ')

%% Data type 3 -> Writing table for all animal species across all networks

T_CoexistRes_A_dataset= array2table(A_dataset,'VariableNames',{'degree','abundA0','abundA1','Dj0','Dj1'});

writetable(T_CoexistRes_A_dataset,sprintf('Coex_A_muAP%d_%dm.txt',muAP,dataset),'Delimiter',' ')

%% Data type 4 -> Writing table for efforts with the corresponding plant and pollinator degrees

T_CoexistRes_Efforts= array2table(didjEff_dataset,'VariableNames',{'degreeP','degreeA','Effort0','Effort1'});

writetable(T_CoexistRes_Efforts,sprintf('Coex_Efforts_muAP%d_%dm.txt',muAP,dataset),'Delimiter',' ')

%% Data type 5 -> Writing table for important parameter values for plants (the ones that appear in calculations of p*)

T_CoexistRes_parsP= array2table(parsP_01,'VariableNames',{'mu_p0','g0','w0','e0','mu_p1','g1','w1','e1'});

writetable(T_CoexistRes_parsP,sprintf('Coex_parsP_muAP%d_%dm.txt',muAP,dataset),'Delimiter',' ')

%% Save array with final alpha matrices
save(sprintf('Alphas%dm.mat',dataset), 'Alphas');

toc
