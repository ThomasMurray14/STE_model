% script to fit models 1-4 using the split obs scripts



%%
clear;
close all;
addpath('prc1');
addpath('prc2');
addpath('split_obs');


%% Load configs
model1_prc_model_config = prc1_ehgf_binary_pu_tbt_config();
model1_resp_obs_model_config = obs1_unitsq_sgm_tbt_config();
model1_RT_obs_model_config = obs1_logrt_linear_binary_config();

model2_prc_model_config = prc2_ehgf_binary_pu_tbt_config(); 
model2_resp_obs_model_config = obs1_unitsq_sgm_tbt_config(); 
model2_RT_obs_model_config = obs1_logrt_linear_binary_config();

model3_prc_model_config = prc1_ehgf_binary_pu_tbt_config(); 
model3_resp_obs_model_config = obs2_unitsq_sgm_tbt_config(); 
model3_RT_obs_model_config = obs1_logrt_linear_binary_config(); 

model4_prc_model_config = tapas_ehgf_binary_config();
model4_resp_obs_model_config = obs3_psychometric_config();
model4_RT_obs_model_config = obs3_gaussianRT_config();

optim_config     = tapas_quasinewton_optim_config(); % optimisation algorithm
optim_config.nRandInit = 5;


%% Set priors (model 1 - perceptual)
model1_prc_model_config.ommu(2)    = -2;
model1_prc_model_config.omsa(2)    = 4;
model1_prc_model_config.rhomu(2)   = 0; % bias towards sad - does work in terms of # responses, but psychometric functions look wrong
model1_prc_model_config.rhosa(2)   = 4;
model1_prc_model_config.logalmu    = log(.1); % perceptual uncertainty
model1_prc_model_config.logalsa    = 2;
model1_prc_model_config = tapas_align_priors(model1_prc_model_config);

%% Set priors (model 1 - responses)
model1_resp_obs_model_config.logzemu = log(1);
model1_resp_obs_model_config.logzesa = 2;
model1_resp_obs_model_config = tapas_align_priors(model1_resp_obs_model_config);

%% Set priors (model 1 - RT)
model1_RT_obs_model_config.beta0mu = 6.5000;
model1_RT_obs_model_config.beta0sa = 4;
model1_RT_obs_model_config.beta1mu = 0;
model1_RT_obs_model_config.beta1sa = 4;
model1_RT_obs_model_config.beta2mu = 0;
model1_RT_obs_model_config.beta2sa = 4;
model1_RT_obs_model_config.beta3mu = 0;
model1_RT_obs_model_config.beta3sa = 4;
model1_RT_obs_model_config.beta4mu = 2;
model1_RT_obs_model_config.beta4sa = 4;
model1_RT_obs_model_config.logsasa = log(.1);
model1_RT_obs_model_config.logsasa = 2;
model1_RT_obs_model_config = tapas_align_priors(model1_RT_obs_model_config);

%% Set priors (model 2 - perceptual)
model2_prc_model_config.ommu(2)    = -2;
model2_prc_model_config.omsa(2)    = 4;
model2_prc_model_config.rhomu(2)   = 0; % bias towards sad
model2_prc_model_config.rhosa(2)   = 4;
model2_prc_model_config.logalmu    = log(.1); % perceptual uncertainty
model2_prc_model_config.logalsa    = 2;
model2_prc_model_config = tapas_align_priors(model2_prc_model_config);

%% Set priors (model 2 - resp)
model2_resp_obs_model_config.logzemu = log(1);
model2_resp_obs_model_config.logzesa = 2;
model2_resp_obs_model_config = tapas_align_priors(model2_resp_obs_model_config);

%% Set priors (model 2 - RT)
model2_RT_obs_model_config.beta0mu = 6.5000;
model2_RT_obs_model_config.beta0sa = 4;
model2_RT_obs_model_config.beta1mu = 0;
model2_RT_obs_model_config.beta1sa = 4;
model2_RT_obs_model_config.beta2mu = 0;
model2_RT_obs_model_config.beta2sa = 4;
model2_RT_obs_model_config.beta3mu = 0;
model2_RT_obs_model_config.beta3sa = 4;
model2_RT_obs_model_config.beta4mu = 2;
model2_RT_obs_model_config.beta4sa = 4;
model2_RT_obs_model_config.logsasa = log(.1);
model2_RT_obs_model_config.logsasa = 2;
model2_RT_obs_model_config = tapas_align_priors(model2_RT_obs_model_config);


%% Set priors (model 3 - perceptual)
model3_prc_model_config.ommu(2)    = -2;
model3_prc_model_config.omsa(2)    = 4;
model3_prc_model_config.rhomu(2)   = 0; % bias towards sad
model3_prc_model_config.rhosa(2)   = 0; %% IMPORTANT TO FIX IN THIS SETUP
model3_prc_model_config.logalmu    = log(.1); % perceptual uncertainty
model3_prc_model_config.logalsa    = 2;
model3_prc_model_config = tapas_align_priors(model3_prc_model_config);


%% Set priors (model 3 - resp)
model3_resp_obs_model_config.logzeta0mu = log(1);
model3_resp_obs_model_config.logzeta0sa = 2;
model3_resp_obs_model_config.zeta1mu = 0;
model3_resp_obs_model_config.zeta1sa = 2;
model3_resp_obs_model_config = tapas_align_priors(model3_resp_obs_model_config);

%% Set priors (model 3 - RT)
model3_RT_obs_model_config.beta0mu = 6.5000;
model3_RT_obs_model_config.beta0sa = 4;
model3_RT_obs_model_config.beta1mu = 0;
model3_RT_obs_model_config.beta1sa = 4;
model3_RT_obs_model_config.beta2mu = 0;
model3_RT_obs_model_config.beta2sa = 4;
model3_RT_obs_model_config.beta3mu = 0;
model3_RT_obs_model_config.beta3sa = 4;
model3_RT_obs_model_config.beta4mu = 2;
model3_RT_obs_model_config.beta4sa = 4;
model3_RT_obs_model_config.logsasa = log(.1);
model3_RT_obs_model_config.logsasa = 2;
model3_RT_obs_model_config = tapas_align_priors(model3_RT_obs_model_config);


%% Set priors (model 4 - perceptual)
model4_prc_model_config.ommu(2) = -3;
model4_prc_model_config.ommu(3) = -3;
model4_prc_model_config.omsa(2) = 0;
model4_prc_model_config.omsa(2) = 0;
model4_prc_model_config = tapas_align_priors(model4_prc_model_config);

%% Set priors (model 4 - resp)
model4_resp_obs_model_config.logitalphamu = tapas_logit(.5, 1);
model4_resp_obs_model_config.logitalphasa = 4;
model4_resp_obs_model_config.logbetamu = log(15);
model4_resp_obs_model_config.logbetasa = 4;
model4_resp_obs_model_config.logitgamma0mu = tapas_logit(.0001, 1);
model4_resp_obs_model_config.logitgamma0sa = 0;
model4_resp_obs_model_config.logitlambda0mu = tapas_logit(.0001, 1);
model4_resp_obs_model_config.logitlambda0sa = 0;
model4_resp_obs_model_config = tapas_align_priors(model4_resp_obs_model_config);

%% Set priors (model 4 - RT)
model4_RT_obs_model_config.logitmumu = tapas_logit(.5, 1);
model4_RT_obs_model_config.logitmusa = 4;
model4_RT_obs_model_config.logsigma0mu = log(.7); % width
model4_RT_obs_model_config.logsigma0sa = 4;
model4_RT_obs_model_config.loggamma1mu = log(0); % baseline logRT - fix to 0
model4_RT_obs_model_config.loggamma1sa = 0;
model4_RT_obs_model_config.loglambda1mu = log(7); % peak of logRT (above baseline)
model4_RT_obs_model_config.loglambda1sa = 4;
model4_RT_obs_model_config.logsigma1mu = log(1); % logRT noise
model4_RT_obs_model_config.logsigma1sa = 4;
model4_RT_obs_model_config = tapas_align_priors(model4_RT_obs_model_config);




%% Fit the data!


% example data (to get contingencies etc)
sub_data = readtable('STE_data\10369536_A_Threat.csv');

% Contingency space
cue = sub_data.Cue_idx;
cue(cue==0) = -1; % balanced contrast coding for the "happy bias"
outcome = sub_data.Outcome_p_sad/100;
u_al = 1 - (0.5*(1+cue)-cue.*outcome);

% model inputs
u_state = [u_al, cue]; % models 1-3
u_stim = sub_data.Outcome_p_sad/100; % model 4

% data
STE_dir = dir('STE_data\*.csv');
N_files = numel(STE_dir);

model_fits(N_files) = struct( ...
    'ID', [], 'group', '', 'condition', '', ...
    'model1_resp_est', struct(), ...
    'model1_RT_est', struct(), ...
    'model2_resp_est', struct(), ...
    'model2_RT_est', struct(), ...
    'model3_resp_est', struct(), ...
    'model3_RT_est', struct(), ...
    'model4_resp_est', struct(), ...
    'model4_RT_est', struct());


% Just to be fancy
N_files_batch = 100;
% completion_times = zeros(N_files_batch, 1);

% loop through datasets
for i = 101:400
    % tic;
    fprintf('\nFitting dataset %i\n', i);
    % if i>1
    %     avg_iter_time = mean(completion_times(1:i));
    %     fprintf('\tAverage iteration time = %1.2fs', avg_iter_time)
    %     estimated_total_time = avg_iter_time * ((N_files_batch-i) + 1);
    %     fprintf('\n\tEstimated completion time = %im, %1.2fs\n\n', floor(estimated_total_time/60), rem(estimated_total_time,60));
    % end
    

    % get ID, group, condition
    f_name = fullfile(STE_dir(i).folder, STE_dir(i).name);
    [~,name,~] = fileparts(f_name);
    tokens = split(name, '_');
    model_fits(i).ID           = str2double(tokens{1});
    model_fits(i).group        = tokens{2};
    model_fits(i).condition    = tokens{3};
    
    % load data
    sub_data = readtable(f_name);

    % model inputs
    u_state_sub = u_state;
    u_stim_sub = u_stim;

    % subject responses
    y_RT = log(sub_data.Response_RT);
    y_resp_state = double(sub_data.Cue_idx == sub_data.Response_idx);
    y_resp_stim = sub_data.Response_idx;

    % remove missing
    missed = isnan(sub_data.Response_idx);
    u_state_sub = u_state_sub(~missed, :);
    u_stim_sub = u_stim_sub(~missed, :);
    y_RT = y_RT(~missed, :);
    y_resp_state = y_resp_state(~missed, :);
    y_resp_stim = y_resp_stim(~missed, :);

    % fit model 1 responses
    try
        model_fits(i).model1_resp_est = tapas_fitModel(...
            y_resp_state,...
            u_state_sub,...
            model1_prc_model_config,...
            model1_resp_obs_model_config,...
            optim_config);
    catch
    end

    % fit model 1 RT
    try
        model_fits(i).model1_RT_est = tapas_fitModel(...
            y_RT,...
            u_state_sub,...
            model1_prc_model_config,...
            model1_RT_obs_model_config,...
            optim_config);
    catch
    end

    % fit model 2 responses
    try
        model_fits(i).model2_resp_est = tapas_fitModel(...
            y_resp_state,...
            u_state_sub,...
            model2_prc_model_config,...
            model2_resp_obs_model_config,...
            optim_config);
    catch
    end

    % fit model 2 RT
    try
        model_fits(i).model2_RT_est = tapas_fitModel(...
            y_RT,...
            u_state_sub,...
            model2_prc_model_config,...
            model2_RT_obs_model_config,...
            optim_config);
    catch
    end

    % fit model 3 responses
    try
        model_fits(i).model3_resp_est = tapas_fitModel(...
            y_resp_state,...
            u_state_sub,...
            model3_prc_model_config,...
            model3_resp_obs_model_config,...
            optim_config);
    catch
    end

    % fit model 3 RT
    try
        model_fits(i).model3_RT_est = tapas_fitModel(...
            y_RT,...
            u_state_sub,...
            model3_prc_model_config,...
            model3_RT_obs_model_config,...
            optim_config);
    catch
    end


    % fit model 4 responses
    try
        model_fits(i).model4_resp_est = tapas_fitModel(...
            y_resp_stim,...
            u_stim_sub,...
            model4_prc_model_config,...
            model4_resp_obs_model_config,...
            optim_config);
    catch
    end

    % fit model 4 RT
    try
        model_fits(i).model4_RT_est = tapas_fitModel(...
            y_RT,...
            u_stim_sub,...
            model4_prc_model_config,...
            model4_RT_obs_model_config,...
            optim_config);
    catch
    end

    % completion_times(i) = toc;
end


save('split_obs_101-400.mat', 'model_fits');




