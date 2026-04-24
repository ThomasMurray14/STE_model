% Function to run parameter recovery for PerceptualBias model

clear; close all; clc
addpath('../helper_functions/')


%% Get inputs
% example data (to get contingencies etc)
sub_data = readtable('..\STE_data\10369536_A_Threat.csv');
[u,y] = data_prep(sub_data);


%% Get configuration structures
[prc_config, obs_config] = STE_PerceptualBias3L_config;
optim_config     = tapas_quasinewton_optim_config(); % optimisation algorithm
optim_config.nRandInit = 5;


%% Fit data
model_fits = fit_master(u, prc_config, obs_config, optim_config);
save('model_PerceptualBias3L_fit2.mat', 'model_fits');
