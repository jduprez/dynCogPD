%% 1. DEFINE GLOBAL VARIABLES (CONTROLS=HC + PATIENTS=PD) (By User***)
% This version of the code works on the concatenation of the PD and HC data

clear all

path_base     = 'E:\DynCogPD'; % path of the main folder 

% HC
n_controls    = 10; % number of controls HC
NCs_HC        = 5; % number of ICA components for HC
path_conn_HC  = 'E:\DynCogPD\Results\conn-PLV\HC'; % path of saved cmat connectivities results for each subject of control grp
path_state_HC = 'E:\DynCogPD\Results\state-ICA\HC'; % path of saved ICA states results for each subject of control grp

% PD
n_parks       = 21; % number of patients PD
NCs_PD        = 5; % number of ICA components for PD
path_conn_PD  = 'E:\DynCogPD\Results\conn-PLV\PD'; % path of saved cmat connectivities results for each subject of park grp
path_state_PD = 'E:\DynCogPD\Results\state-ICA\PD'; % path of saved ICA states results for each subject of park grp

band_interval=[30 40]; % band of interest


%% 2. DEFINE CMAT LIST AND LOAD GROUP ICA + PERMS RESULTS FOR HC + PD

folder1=[path_base '\Code']; addpath(genpath(folder1)); % add Code folder and subfolders

% HC
for i = 1:n_controls
    cmat_list_HC{i} = [path_conn_HC '\cmat_PLV_' int2str(i) '.mat']; % cmat_list
end
load([path_state_HC '\results_ICA.mat']); % ICA results on HC 
load([path_state_HC '\perms.mat']); % perms results on HC 
results_ICA_HC = results_ICA; 
results_ICA_HC.NCs = NCs_HC;
perms_HC = perms; 
clear results
clear perms

% PD
for i = 1:n_parks
    cmat_list_PD{i} = [path_conn_PD '\cmat_PLV_' int2str(i) '.mat']; % cmat_list
end
load([path_state_PD '\results_ICA.mat']); % ICA results on PD 
load([path_state_PD '\perms.mat']); % perms results on PD 
results_ICA_PD = results_ICA;
results_ICA_PD.NCs = NCs_PD;
perms_PD = perms;
clear results
clear perms

% Determine the index of onset time, should be the same between grps
ind_0s = find(results_ICA_HC.time==0); 
if(isempty(ind_0s))
    [mmin,ind_0s] = min(abs(results_ICA_HC.time));
end


%% 3. EXTRACT AUTOMATICALLY SIGNIFICANT STATES FOR HC + PD

% 3.1. Define minimum duration for significance. 

ncycles = 3;
d_cy = ncycles*(round(1000/band_interval(1)));


% 3.2. Extract significant states with corresponding significance time (based on null distribution + 3 cycles surviving)

[isSignif_NCs_HC,timeSignif_HC] = isSignif(results_ICA_HC,perms_HC,NCs_HC,ind_0s, d_cy);
[isSignif_NCs_PD,timeSignif_PD] = isSignif(results_ICA_PD,perms_PD,NCs_PD,ind_0s, d_cy);
kept_net                        = [isSignif_NCs_HC,isSignif_NCs_PD];


% 3.3. Store kept significant maps for each group

cHC=0; cPD=0;
% HC
for i=1:NCs_HC
    if(isSignif_NCs_HC{i})
        cHC = cHC+1;
        states_maps_HC(:,:,cHC) = results_ICA_HC.maps(:,:,i);
    end
end
% PD
for i=1:NCs_PD
    if(isSignif_NCs_PD{i})
        cPD = cPD+1;
        states_maps_PD(:,:,cPD) = results_ICA_PD.maps(:,:,i);
    end
end


% 3.4. Combine all significant states maps (for HC followed by PD) in one variable: states_maps_all
nROI = size(results_ICA_HC.maps,1);
states_maps_all = zeros(nROI,nROI,cHC+cPD);
states_maps_all(:,:,1:cHC) = states_maps_HC(:,:,1:cHC);
states_maps_all(:,:,cHC+1:cHC+cPD) = states_maps_PD(:,:,1:cPD);


%% 4. APPLY BACKFITTING ALGORITHM FOR HC + PD

% 4.1. Configuration for backfitting algo

cfg_algo                 = [];
cfg_algo.threshnet_meth  = 'no'; % only one choice is implemented in the code
cfg_algo.corr_meth       = 'corr2'; % correlation as spatial similarity measure
cfg_algo.cHC             = cHC;
cfg_algo.cPD             = cPD;
cfg_algo.states_maps_all = states_maps_all;


% 4.2. Run backfitting algo for HC + PD

[corr_tw_HC,max_tw_HC,ind_tw_HC,cmat_allHC] = do_backfitting(cfg_algo,cmat_list_HC,n_controls);
[corr_tw_PD,max_tw_PD,ind_tw_PD,cmat_allPD] = do_backfitting(cfg_algo,cmat_list_PD,n_parks);


%% 5. CALCULATE MICROSTATS PARAMS FOR HC + PD

% 5.1. Configuration for Microstats extraction

cfg_ms           = [];
cfg_ms.cHC       = cHC; % number of significant kept states in HC grp
cfg_ms.cPD       = cPD; % number of significant kept states in PD grp
cfg_ms.ind_0s    = ind_0s; % onset time
cfg_ms.totaltime = results_ICA_HC.time(end); % total time duration after onset (in sec)
cfg_ms.deltat    = results_ICA_HC.time(end)-results_ICA_HC.time(end-1); % difference time between two time windows (in sec)


% 5.2. Extract Microstats for HC + PD

microparams_HC = extract_microstates(cfg_ms,corr_tw_HC,ind_tw_HC,cmat_allHC,n_controls);
microparams_PD = extract_microstates(cfg_ms,corr_tw_PD,ind_tw_PD,cmat_allPD,n_parks);

