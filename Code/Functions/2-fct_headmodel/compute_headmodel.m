function subvol = compute_headmodel(inparam)

% This function calculates headmodel from MRI, head, inner/outerskull structures using OpenMEEG FieldTrip toolbox
%
% Input:
%   inparam : structure of useful input parameters (namely for path_base and mri_realign variables)
%
% Output:
%   subvol  : headmodel structure, fieldTrip format



% 1. Load realigned MRI saved in Inputs (Refer to code_for_inputs)
load([inparam.path_base '\Inputs\' inparam.mri_template '\mri_' inparam.mri_template '_realign.mat']); 

% 2. Segment MRI
cfg           = [];
cfg.output    = {'brain','skull','scalp'};
segmentedmri  = ft_volumesegment(cfg, mri_realign); %ctf/mm

% 3. Add path for OpenMEEG toolbox
setenv('PATH', [inparam.path_base '\TOOLBOXES\OpenMEEG\bin']); %add path for openmeeg

% 4. Prepare mesh using 'headshape' method
cfg=[];
cfg.method='headshape';
brain = ft_read_headshape([path_base '\input4code\' mri_template '\tess_innerskull_' mri_template '.mat']); % load the innerskull template used in Brainstorm 
brain_mm=ft_convert_units(brain,'mm');
cfg.headshape=brain_mm;
cfg.numvertices = [3000];
bnd(1)=ft_prepare_mesh(cfg,segmentedmri);

cfg=[];
cfg.method='headshape';
skull=ft_read_headshape([inparam.path_base '\input4code\' inparam.mri_template '\tess_outerskull_' inparam.mri_template '.mat']); %load the outerskull template used in Brainstorm 
skull_mm=ft_convert_units(skull,'mm');
cfg.headshape=skull_mm;
cfg.numvertices = [3000];
bnd(2)=ft_prepare_mesh(cfg,segmentedmri);

cfg=[];
cfg.method='headshape';
head=ft_read_headshape([inparam.path_base '\input4code\' inparam.mri_template '\tess_head_' inparam.mri_template '.mat']); %load the head template used in Brainstorm 
head_mm=ft_convert_units(head,'mm');
cfg.headshape=head_mm;
cfg.numvertices = [3000];
bnd(3)=ft_prepare_mesh(cfg,segmentedmri);

% %plot mesh
% figure();
% ft_plot_mesh(bnd(1), 'edgecolor', 'none', 'facecolor', 'r')
% ft_plot_mesh(bnd(2), 'edgecolor', 'none', 'facecolor', 'g')
% ft_plot_mesh(bnd(3), 'edgecolor', 'none', 'facecolor', 'b')
% alpha 0.3

% 5. Use OpenMEEG to build subvol (headmodel)
cfg        = [];
cfg.method ='openmeeg'; 
subvol     = ft_prepare_headmodel(cfg, bnd); %ctf/mm

end

