% Function to run parameter recovery for PerceptualBias model

clear; close all; clc
addpath('../helper_functions/');

%% Get inputs
% example data (to get contingencies etc)
sub_data = readtable('..\STE_data\10369536_A_Threat.csv');
[u,~] = data_prep(sub_data);


%% Get configuration structures
[prc_config, obs_config] = STE_ResponseBias_config;
optim_config     = tapas_quasinewton_optim_config(); % optimisation algorithm
optim_config.nRandInit = 5;


%% Run parameter recovery

% run recovery
N=200;

% Parameters to recover
prc_param_names = {'om2', 'al'};
prc_param_idx   = [13, 15];
prc_param_space = {'native', 'log'};
obs_param_names = {'zeta0', 'zeta1', 'beta0', 'beta1', 'beta2', 'beta3', 'beta4', 'sa'};
obs_param_idx   = [1, 2, 3, 4, 5, 6, 7, 8];
obs_param_space = {'log', 'native', 'native', 'native', 'native', 'native', 'native', 'log'};

recov = parameter_recovery_master( ...
    u,...
    prc_config,...
    obs_config,...
    optim_config,...
    N,...
    prc_param_names,...
    prc_param_idx,...
    prc_param_space,...
    obs_param_names,...
    obs_param_idx,...
    obs_param_space);

save('ResponseBias_recovery.mat', 'recov');
recovery_figures(recov);



