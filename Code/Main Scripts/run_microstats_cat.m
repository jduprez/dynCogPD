%% 1. DEFINE GLOBAL VARIABLES
% This version of the code works on the concatenation of the PD and HC data

clear all

% add functions to path

addpath(genpath('C:\GitHub\dynCogPD\Code'));

frequency = 'beta'; % or gamma
CATfile = dir(['F:\WJD\Simon Dynamic FC\Results\FCmat\CAT\', frequency, '\', '*_incong_wplv1225.mat']);

path_base     = 'F:\DynCogPD'; % path of the main folder

% Define the controls and patients indices

HC = [1, 3, 6, 10, 11, 13, 16, 17, 18, 19];
PD = [2, 4, 5, 7, 8, 9, 12, 14, 15, 20];

nsub = size(HC, 2) + size(PD, 2);

NCs = 5;

path_states = 'F:\WJD\Simon Dynamic FC\Results\ICA\HC_PD_CAT'; % path of the ICA results
path_conn = ['F:\WJD\Simon Dynamic FC\Results\FCmat\CAT\', frequency];
%
%
% % HC
% n_controls    = 10; % number of controls HC
% NCs_HC        = 5; % number of ICA components for HC
% path_conn_HC  = 'E:\DynCogPD\Results\conn-PLV\HC'; % path of saved cmat connectivities results for each subject of control grp
% path_state_HC = 'E:\DynCogPD\Results\state-ICA\HC'; % path of saved ICA states results for each subject of control grp
%
% % PD
% n_parks       = 21; % number of patients PD
% NCs_PD        = 5; % number of ICA components for PD
% path_conn_PD  = 'E:\DynCogPD\Results\conn-PLV\PD'; % path of saved cmat connectivities results for each subject of park grp
% path_state_PD = 'E:\DynCogPD\Results\state-ICA\PD'; % path of saved ICA states results for each subject of park grp
%

if strcmp(frequency, 'beta')
    band_interval = [12 25];
else
    band_interval=[30 45]; % band of interest
end


%% 2. DEFINE CMAT LIST AND LOAD GROUP ICA + PERMS RESULTS FOR HC + PD

% folder1=[path_base '\Code']; addpath(genpath(folder1)); % add Code folder and subfolders

% HC
for i = 1:size(HC, 2)
    cmat_list_HC{i} = [path_conn '\' CATfile(HC(i)).name]; % cmat_list
end

% PD
for i = 1:size(PD, 2)
    cmat_list_PD{i} = [path_conn '\' CATfile(PD(i)).name]; % cmat_list
end

if strcomp(frequency, 'beta')
    load([path_state '\CAT_IC_plvdyn1225_5IC.mat']); % ICA results on HC
    load([path_state_HC '\perms_cat_beta.mat']); % perms results on HC
else
    load([path_state '\CAT_IC_plvdyn3045_5IC.mat']); % ICA results on HC
    load([path_state_HC '\perms_cat_gamma.mat']); % perms results on HC
end


% Determine the index of onset time, should be the same between grps
ind_0s = find(results.time==0);
if(isempty(ind_0s))
    [mmin,ind_0s] = min(abs(results.time));
end


%% 3. EXTRACT AUTOMATICALLY SIGNIFICANT STATES FOR HC + PD

% 3.1. Define minimum duration for significance.

ncycles = 3;
d_cy = ncycles*(round(1000/band_interval(1)));


% 3.2. Extract significant states with corresponding significance time (based on null distribution + 3 cycles surviving)

[isSignif_NCs,timeSignif] = isSignif(results,perms,NCs,ind_0s, d_cy);

kept_net = [isSignif_NCs];


% 3.3. Store kept significant maps

cCAT=0; 
% HC
for i=1:NCs
    if(isSignif_NCs{i})
        cCAT = cCAT+1;
        states_maps(:,:,cCAT) = results.maps(:,:,i);
    end
end

% 3.4. Combine all significant states maps (for HC followed by PD) in one variable: states_maps_all
nROI = size(results.maps,1);
states_maps_all = zeros(nROI,nROI,cCAT);
states_maps_all(:,:,1:cCAT) = states_maps(:,:,1:cCAT);
% states_maps_all(:,:,cHC+1:cHC+cPD) = states_maps_PD(:,:,1:cPD);


%% 4. APPLY BACKFITTING ALGORITHM FOR HC + PD

% 4.1. Configuration for backfitting algo

cfg_algo                 = [];
cfg_algo.threshnet_meth  = 'no'; % only one choice is implemented in the code
cfg_algo.corr_meth       = 'corr2'; % correlation as spatial similarity measure
cfg_algo.cCAT            = cCAT;
% cfg_algo.cPD             = cPD;
cfg_algo.cHC = cCAT;
cfg_algo.cPD = 0;
cfg_algo.states_maps_all = states_maps_all;


% 4.2. Run backfitting algo for HC + PD

[corr_tw_HC,max_tw_HC,ind_tw_HC,cmat_allHC] = do_backfitting(cfg_algo,cmat_list_HC,size(HC, 2));
[corr_tw_PD,max_tw_PD,ind_tw_PD,cmat_allPD] = do_backfitting(cfg_algo,cmat_list_PD,size(PD, 2));


%% 5. CALCULATE MICROSTATS PARAMS FOR HC + PD

% 5.1. Configuration for Microstats extraction

cfg_ms           = [];
cfg_ms.cHC       = cCAT; % number of significant kept states in HC grp
cfg_ms.cPD       = 0; % number of significant kept states in PD grp
cfg_ms.cCAT       = cCAT; % number of significant kept states in PD grp

cfg_ms.ind_0s    = ind_0s; % onset time
cfg_ms.totaltime = results.time(end); % total time duration after onset (in sec)
cfg_ms.deltat    = results.time(end)-results.time(end-1); % difference time between two time windows (in sec)


% 5.2. Extract Microstats for HC + PD

microparams_HC = extract_microstates(cfg_ms,corr_tw_HC,ind_tw_HC,cmat_allHC,size(HC, 2));
microparams_PD = extract_microstates(cfg_ms,corr_tw_PD,ind_tw_PD,cmat_allPD,size(PD, 2));



% 5.3. Extract parameters

n_net = size(microparams_HC.fraction_covtime, 2);

frac_covtime_HC = zeros(10, n_net);

for neti = 1:n_net
    frac_covtime_HC(:, neti) = microparams_HC.fraction_covtime{1, neti};
end

frac_covtime_PD = zeros(21, n_net);

for neti = 1:n_net
    frac_covtime_PD(:, neti) = microparams_PD.fraction_covtime{1, neti};
end


%%

freq_occurence_HC = zeros(10, n_net);

for neti = 1:n_net
    freq_occurence_HC(:, neti) = microparams_HC.freq_occurence{1, neti};
end

freq_occurence_PD = zeros(21, n_net);

for neti = 1:n_net
    freq_occurence_PD(:, neti) = microparams_PD.freq_occurence{1, neti};
end

%%

avg_lifespan_HC = zeros(10, n_net);

for neti = 1:n_net
    avg_lifespan_HC(:, neti) = microparams_HC.avg_lifespan{1, neti};
end

avg_lifespan_PD = zeros(21, n_net);

for neti = 1:n_net
    avg_lifespan_PD(:, neti) = microparams_PD.avg_lifespan{1, neti};
end%%

%%
GEV_HC = zeros(10, n_net);

for neti = 1:n_net
    GEV_HC(:, neti) = microparams_HC.GEV{1, neti};
end

GEV_PD = zeros(21, n_net);

for neti = 1:n_net
    GEV_PD(:, neti) = microparams_PD.GEV{1, neti};
end

%%
TR_HC = zeros(n_net, n_net, 10);

for subi = 1:10
    TR_HC(:,:,subi) = microparams_HC.TR{1, subi};
end

TR_PD = zeros(n_net, n_net, 21);

for subi = 1:21
    TR_PD(:,:,subi) = microparams_PD.TR{1, subi};
end
%%
TRsym0_HC = zeros(n_net, n_net, 10);

for subi = 1:10
    TRsym0_HC(:,:,subi) = microparams_HC.TRsym0{1, subi};
end

TRsym0_PD = zeros(n_net, n_net, 21);

for subi = 1:21
    TRsym0_PD(:,:,subi) = microparams_PD.TRsym0{1, subi};
end

