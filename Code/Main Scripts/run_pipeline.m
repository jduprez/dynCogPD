%% 1. DEFINE VARIABLES + AUTOMATICALLY ADD NECESSARY PATHS AND VARIABLES

clear all

% 1.1. DEFINE GENERAL PARAMS (By User***)

path_base     = 'E:\DynCogPD'; % path of the main folder 
path_conn     = 'E:\DynCogPD\Results\conn-PLV\HC'; % path of connectivity (PLV) results folder to save in 
path_state    = 'E:\DynCogPD\Results\state-ICA\HC'; % path of states (ICA) results folder to save in
mri_template  = 'icbm'; % 'icbm' for ICBM152 (only this option is implemented in code)
atlas         = 'destrieux'; % 'destrieux' for Destrieux148 atlas, format of BrainStorm (only this option is implemented in code)
nROIs         = 148; % number of ROIs in destrieux atlas


% 1.2. AUTOMATICALLY ADD NECESSARY FOLDERS

folder1=[path_base '\Inputs']; addpath(genpath(folder1)); % add Inputs folder and subfolders
folder2=[path_base '\Code']; addpath(genpath(folder2)); % add Code folder and subfolders
setenv('PATH', [path_base '\TOOLBOXES\OpenMEEG\bin']);
addpath([path_base '\Toolboxes\fieldtrip-20190224\external\openmeeg']); % add openmeeg toolbox
tmp = which('ft_defaults');
if isempty(tmp)
    addpath([path_base '\Toolboxes\fieldtrip-20190224']); % add defaults fieldtrip
    ft_defaults
end


% 1.3. LOAD NECESSARY INPUTS AND BUILD 'INPARAM' STRUCTURE TO BE USED 
% (already computed and saved variables, as we are dealing with template MRI;
% these variables were computed for once (refer to code_for_inputs.m script in Code Folder)

load([path_base '\Inputs\' mri_template '\elec_BS_' mri_template '.mat']); % fieldtrip format electrode extracted from BrainStorm, see Section 4 in code_for_inputs.m
load([path_base '\Inputs\' mri_template '\' atlas '\' atlas '_scs_' mri_template '.mat']); %scout structure to extract sources centroids positions orientations (in SCS or CTF coord): calculated in Section 2 of 'code_for_inputs'
load([path_base '\Inputs\' mri_template '\' atlas '\' atlas '_subgrid.mat']); % fieldtrip format of the grid specific to the mri and scout used (refer to code_for_inputs.m for details about computation): calculated in Section 3 of 'code_for_inputs'
load([path_base '\Inputs\' mri_template '\' atlas '\' atlas '_labels.mat']); % cell labels for the scout regions: extracted from BrainStorm
cfg_inparam.path_base    = path_base;
cfg_inparam.path_conn    = path_conn;
cfg_inparam.path_state   = path_state;
cfg_inparam.mri_template = mri_template;
cfg_inparam.atlas        = atlas; 
cfg_inparam.nROIs        = nROIs;
cfg_inparam.elec_BS_mm   = elec_BS_mm;
cfg_inparam.subgrid      = subgrid;
cfg_inparam.scout_labels = scout_labels;
cfg_inparam.scout_scs    = scout_scs;


%% 2. CREATE DATA STRUCTURE  

% 2.1. Configuration for database necessary to create our data Structure (By User***)

load([path_base '\Inputs\indata.mat']); % refer to code_for_inputs.m, Section 5
load([path_base '\Inputs\nb_trials_persub.mat']); % refer to code_for_inputs.m, Section 5

cfg_database.indata       = indata; % cell of length (nb_sub), each cell is a cell set of trials represented by a 2D matrix of dim [nb_chan,nb_samples]
cfg_database.nb_sub       = 2; % number of subjects
cfg_database.nb_trials    = nb_trials_persub; % number of trials for each subject
cfg_database.nb_chan      = 32; % number of channels (electrodes)
cfg_database.fs           = 1024; % sampling frequency
cfg_database.onset        = 25; % sample index for onset 
cfg_database.pre_samples  = 25; % nb samples before trial onset 
cfg_database.post_samples = 100; % nb samples after trial onset 


% 2.2. Create your dataStruct (fieldtrip format)

dataStruct=create_dataStruct(cfg_database,cfg_inparam);
cfg_database=rmfield(cfg_database,'indata');
clear indata


%% 3. COMPUTE HEADMODEL (OpenMEEG)

cfg_inparam.subvol = compute_headmodel(cfg_inparam); %subvol is the headmodel struct


%% 4. COMPUTE SOURCES WMNE + DYNAMIC CONNECTIVITY PLV DYNAMIQUE

% 4.1. Configuration for source and dFC (dynamic Functional Connectivity) structures (By User***)

% Configuration for source structure 
cfg_source.weightExp   = 0.5; % param for depth weighting
cfg_source.weightLimit = 10; % param for depth weighting limit
cfg_source.SNR         = 3; % Signal to Noise Ratio
% Configuration for dFC structure
cfg_dFC.conn_method = 'plv_dyn'; % 'plv_dyn' for windowedPLV  (only this option is implemented in code)
cfg_dFC.bpfreq      = [30 40]; % frequency band of interest (gamma/beta..)
cfg_dFC.window.size = 0.17; % sliding window length in seconds (for example calculated for 6cycles,CentralFreq=35 ==> 6/35=0.17s) 
cfg_dFC.window.step = 0.017; % step between windows in seconds (for example 90% overlapping=10/100*window_size)


% 4.2. Compute source connectivity 

[filters_grp,cmat_grp] = compute_source_connectivity(cfg_source,cfg_dFC,cfg_database,cfg_inparam,dataStruct);
cfg_inparam=rmfield(cfg_inparam,'subvol');


%% 5. COMPUTE ICA 

% 5.0. Determine optimal number of ICs to look for usign DIFFIT



% 5.1. Configuration for states (ICA algo) (By User***)

NCs                 = 5; % to specify by the user
cfg_state           = [];
cfg_state.NCs       = NCs; 
cfg_state.conn      = cmat_grp;
cfg_state.n_parcels = cfg_inparam.nROIs;


% 5.2. Run ICA algo

tic
results_ICA = go_decomposeConnectome_SS_ICA(cfg_state); %all
timeElapsed = toc


% 5.3. Save ICA states results in the specified path_state folder

save([path_state '\results_ICA.mat'],'results_ICA');


% 5.4. Uncomment below if you want to visualize ICA states results

% load([cfg_inparam.path_base '\Inputs\' cfg_inparam.mri_template '\Surfmatrix_' cfg_inparam.mri_template '.mat']); %structure extracted from brainnetviewer to draw the brain cortex
% load([cfg_inparam.path_base '\Inputs\' cfg_inparam.mri_template '\' cfg_inparam.atlas '\' cfg_inparam.atlas '_mni_' cfg_inparam.mri_template '.mat']); %scout structure to extract sources centroids positions orientations (in MNI coord)
% cfg_inparam.scout_mni=scout_mni;
% cfg_inparam.Surfmatrix=Surfmatrix;
% cfg_plot              = [];
% cfg_plot.meth         = 'thresh_conn'; % 4 choices for thresholding: 'thresh_conn' or 'thresh_node' or 'thresh_abs' or 'thresh_val'
% cfg_plot.threshold    = 0.995; % threshold value
% cfg_plot.label        = 0; %1:if yes display labels, 0:if no display labels
% cfg_plot.scout_labels = cfg_inparam.scout_labels;
% cfg_plot.scout_mni    = cfg_inparam.scout_mni;
% cfg_plot.Surfmatrix   = cfg_inparam.Surfmatrix;
% go_viewNetworkComponents_eeg_concat_interface(cfg_plot,results_ICA); 


%% 6. Apply NULL DISTRIBUTION 

% 6.1. Configuration for null distribution (Sign-Flip approach) (By User***)

cfg_null                   = [];
cfg_null.p                 = 0.05;
cfg_null.bonferroni_factor = 2*NCs; % 2 tails x NCs ICs


% 6.2. Run permutation

perms = go_testNetworks_general(cfg_null,results_ICA);


% 6.3. Save permutations results in the specified path_state

save([path_state '\perms.mat'],'perms');


% 6.4. Uncomment below if you want to visualize ICA states results

% load([cfg_inparam.path_base '\Inputs\' cfg_inparam.mri_template '\Surfmatrix_' cfg_inparam.mri_template '.mat']); %structure extracted from brainnetviewer to draw the brain cortex
% load([cfg_inparam.path_base '\Inputs\' cfg_inparam.mri_template '\' cfg_inparam.atlas '\' cfg_inparam.atlas '_mni_' cfg_inparam.mri_template '.mat']); %scout structure to extract sources centroids positions orientations (in MNI coord)
% cfg_inparam.scout_mni=scout_mni;
% cfg_inparam.Surfmatrix=Surfmatrix;
% cfg_plot              = [];
% cfg_plot.meth         = 'thresh_conn'; %***  4 choices for thresholding: 'thresh_conn' or 'thresh_node' or 'thresh_val' or 'thresh_abs'
% cfg_plot.threshold    = 0.995; %*** threshold value
% cfg_plot.label        = 0; %1:if yes display labels, 0:if no display labels
% cfg_plot.scout_labels = cfg_inparam.scout_labels;
% cfg_plot.scout_mni    = cfg_inparam.scout_mni;
% cfg_plot.Surfmatrix   = cfg_inparam.Surfmatrix;
% go_viewNetworkComponents_eeg_concat_interface(cfg_plot,results_ICA,perms);

