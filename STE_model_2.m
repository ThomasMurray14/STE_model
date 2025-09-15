%% MODEL 2
% Perceptual model 
%   = [prc2] 2 level eHGF Nace remix (bias on perception)
% Response model 
%   = [obs1] combined unitsq_sgm and logRT



%%
clear;
close all;
addpath('prc2');
addpath('obs1');

%% 

% example data (to get contingencies etc)
sub_data = readtable('STE_data\10369536_A_Threat.csv');

% Contingency space
cue = sub_data.Cue_idx;
cue(cue==0) = -1; % balanced contrast coding for the "happy bias"
outcome = sub_data.Outcome_p_sad/100;
u_al = 1 - (0.5*(1+cue)-cue.*outcome);
state = double(sub_data.Cue_idx == sub_data.Outcome_idx);

sub_data.u_al = u_al;
sub_data.state=state;
sub_data.p_sad=outcome;


% Responses
sub_data.logRT = log(sub_data.Response_RT);
sub_data.resp_state = double(sub_data.Cue_idx == sub_data.Response_idx);

% model input
u = [sub_data.u_al, cue];
y = [sub_data.resp_state, sub_data.logRT];%, sub_data.Confidence_idx];


%% Get configuration structures
prc_model_config = prc2_ehgf_binary_pu_tbt_config(); % perceptual model
obs_model_config = obs1_comb_obs_config(); % response model
optim_config     = tapas_quasinewton_optim_config(); % optimisation algorithm
optim_config.nRandInit = 5;

%% Set parameters and priors
prc_model_config.ommu(2)    = -2;
prc_model_config.omsa(2)    = 4;

prc_model_config.rhomu(2)   = 0; % bias towards sad
prc_model_config.rhosa(2)   = 4;

prc_model_config.logalmu    = log(.1); % perceptual uncertainty
prc_model_config.logalsa    = 2;

prc_model_config = tapas_align_priors(prc_model_config);

obs_model_config.logzemu = log(1);
obs_model_config.logzesa = 2;

obs_model_config.beta0mu = 6.5000;
obs_model_config.beta0sa = 4;

obs_model_config.beta1mu = 0;
obs_model_config.beta1sa = 4;

obs_model_config.beta2mu = 0;
obs_model_config.beta2sa = 4;

obs_model_config.beta3mu = 0;
obs_model_config.beta3sa = 4;

obs_model_config.beta4mu = 2;
obs_model_config.beta4sa = 4;

obs_model_config.logsasa = log(.1);
obs_model_config.logsasa = 2;

obs_model_config = tapas_align_priors(obs_model_config);

%% simulate responses

r_temp = [];
r_temp.c_prc.n_levels = 3;
prc_params = prc1_ehgf_binary_pu_tbt_transp(r_temp, prc_model_config.priormus);

obs_params = obs_model_config.priormus;
obs_params(1) = exp(obs_params(1));
obs_params(7) = exp(obs_params(7));
% obs_params(8) = 3;%exp(obs_params(8));
sim = tapas_simModel(u,...
    'prc2_ehgf_binary_pu_tbt',...
    prc_params,...
    'obs1_comb_obs',...
    obs_params,...
    123456789);


sim_sad = (sub_data.Cue_idx == 1 & sim.y(:,1) == 1) + (sub_data.Cue_idx == 0 & sim.y(:,1) == 0);
N_sad = sum(sim_sad);


visualise_psychometric(u, sub_data, 'prc2_ehgf_binary_pu_tbt', prc_params, 'obs1_comb_obs', obs_params, 20)



%% Plot trajectory
prc2_ehgf_binary_tbt_plotTraj(sim);


%% recover parameters

est = tapas_fitModel(...
    sim.y,...
    sim.u,...
    prc_model_config,...
    obs_model_config,...
    optim_config);



%% Full parameter recovery

% set N iterations
N = 200;

% Parameters to recover
prc_param_names = {'om2', 'rho2', 'al'};
prc_param_idx   = [13, 8, 15];
prc_param_space = {'native', 'native', 'log'};
obs_param_names = {'ze', 'beta0', 'beta1', 'beta2', 'beta3', 'beta4', 'sa'};
obs_param_idx   = [1, 2, 3, 4, 5, 6, 7];
obs_param_space = {'log', 'native', 'native', 'native', 'native', 'native', 'log'};

recov = parameter_recovery_master(u,...
    prc_model_config,...
    obs_model_config,...
    optim_config,...
    N,...
    prc_param_names,...
    prc_param_idx,...
    prc_param_space,...
    obs_param_names,...
    obs_param_idx,...
    obs_param_space);
save('model2_recovery.mat', 'recov');
recovery_figures(recov);


%% Fit actual data

model_fits = fit_master(u, prc_model_config, obs_model_config, optim_config);
save('model2_fit.mat', 'model_fits');





