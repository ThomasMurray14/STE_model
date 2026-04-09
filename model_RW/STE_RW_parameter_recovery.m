% Parameter recovery for RW combined model

%%
clear; close all; clc
addpath('../helper_functions/');

%% Get inputs
% example data (to get contingencies etc)
sub_data = readtable('..\STE_data\10369536_A_Threat.csv');
[u,~] = data_prep(sub_data);


%% Get configuration structures
[prc_config, obs_config] = STE_RW_config;
optim_config     = tapas_quasinewton_optim_config();
optim_config.nRandInit = 10; % Annoying but better chance of fitting


%% Run parameter recovery

% number of iterations
N=200;

% Parameters to recover
prc_param_names = {'al'};
prc_param_idx   = [2];
prc_param_space = {'logit'};
obs_param_names = {'ze', 'beta0', 'beta1', 'beta2', 'sa'};
obs_param_idx   = [1, 2, 3, 4, 5, 6];
obs_param_space = {'log', 'native', 'native', 'native',  'log'};

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

save('model_RW_recovery.mat', 'recov');
recovery_figures(recov);

