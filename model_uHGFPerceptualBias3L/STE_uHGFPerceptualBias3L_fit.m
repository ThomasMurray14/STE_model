% Function to run parameter recovery for PerceptualBias model

clear; close all; clc
addpath('..\helper_functions\')

toolboxroot = 'C:\Users\Tom\Documents\MATLAB\Toolboxes\hgf-toolbox-main\hgf-toolbox-main';
run(fullfile(toolboxroot, 'setup.m'));

%% Get inputs
% example data (to get contingencies etc)
sub_data = readtable('..\STE_data\10369536_A_Threat.csv');
[u,y] = data_prep(sub_data);


%% Get configuration structures
[prc_config, obs_config] = STE_uHGFPerceptualBias3L_config;
optim_config     = quasinewton_optim_config(); % optimisation algorithm
optim_config.nRandInit = 5;


%% Fit data
model_fits = fit_master_parallel(u, prc_config, obs_config, optim_config, 8);
save('model_uHGFPerceptualBias3L_fit.mat', 'model_fits');
