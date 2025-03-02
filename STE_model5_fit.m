% Script to fit model 5 to the data

clear;
close all;

%% Get model input
addpath('custom_hgf_2');
addpath('vkf')

% example data (to get contingencies etc)
sub_data = readtable('STE_data\10369536_A_Threat.csv');

% input in contingency space
state = double(sub_data.Cue_idx == sub_data.Outcome_idx);

% stimulus intensity
sub_data.p_sad = sub_data.Outcome_p_sad/100;

% model input
u = [state, sub_data.p_sad];


%% Get model configs and set priors
% load model configs
prc_model_config = prc_vkf_binary_config(); % perceptual model
obs_model_config = obs3_psychometric_config(); % response model
optim_config     = tapas_quasinewton_optim_config(); % optimisation algorithm
optim_config.nRandInit = 10;

% set priors
prc_model_config.logv0mu = log(.2);
prc_model_config.logv0sa = 0;
prc_model_config.logw0mu = log(.2);
prc_model_config.logw0sa = 0;
prc_model_config = tapas_align_priors(prc_model_config);

obs_model_config.b0sa=2;
obs_model_config.b1sa=2;
obs_model_config.alphasa=1;
obs_model_config = tapas_align_priors(obs_model_config);



%% Loop through data and fit model to each

STE_dir = dir('STE_data\*.csv');
N_files = numel(STE_dir);
model_fits(N_files) = struct('ID', [], 'group', '', 'condition', '', 'model_fit', struct());

for i = 1:N_files
    % get ID, group, condition
    f_name = fullfile(STE_dir(i).folder, STE_dir(i).name);
    [~,name,~] = fileparts(f_name);
    tokens = split(name, '_');
    model_fits(i).ID           = str2double(tokens{1});
    model_fits(i).group        = tokens{2};
    model_fits(i).condition    = tokens{3};
    
    % load data
    sub_data = readtable(f_name);

    % model input
    u_sub = u;

    % subject responses
    y = sub_data.Response_idx; % 1 = sad, 0 = happy

    % remove missing
    missed = isnan(sub_data.Response_idx);
    u_sub = u_sub(~missed, :);
    y = y(~missed, :);

    % fit model
    try
        model_fits(i).est = tapas_fitModel(...
            y,...
            u_sub,...
            prc_model_config,...
            obs_model_config,...
            optim_config);
    catch
    end
end

save('STE_model5_fits.mat', 'model_fits');








