%% MODEL 1
% Perceptual model = 2 level eHGF Nace remix
% Response model = combined unitsq_sgm and logRT
% 
% TODO
% 
% - Try including confidence
% 
% Notes:
% 
% At the moment, rho(2) is getting multiplied by stim_noise - expectation
% that biases will be greater for ambiguous expressions (and reduced for
% clear expressions)



%%
clear;
close all;
addpath('STE_model1');

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
prc_model_config = prc1_ehgf_binary_pu_tbt_config(); % perceptual model
obs_model_config = obs1_comb_obs_config(); % response model
optim_config     = tapas_quasinewton_optim_config(); % optimisation algorithm
optim_config.nRandInit = 5;

%% Set parameters and priors
prc_model_config.ommu(2)    = -2;
prc_model_config.omsa(2)    = 4;

prc_model_config.rhomu(2)   = 0; % bias towards sad - does work in terms of # responses, but psychometric functions look wrong
prc_model_config.rhosa(2)   = 4;

prc_model_config.logalmu    = log(0.1); % perceptual uncertainty
prc_model_config.logalsa    = 4;

prc_model_config = tapas_align_priors(prc_model_config);

obs_model_config.logzemu = log(1);
obs_model_config.logzesa = 2;

obs_model_config.beta0mu = 6.5000;
obs_model_config.beta0sa = 2;

obs_model_config.beta1mu = 0;
obs_model_config.beta1sa = 2;

obs_model_config.beta2mu = 10;
obs_model_config.beta2sa = 2;

obs_model_config.beta3mu = 0;
obs_model_config.beta3sa = 2;

obs_model_config.beta4mu = 0;
obs_model_config.beta4sa = 2;

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
    'prc1_ehgf_binary_pu_tbt',...
    prc_params,...
    'obs1_comb_obs',...
    obs_params,...
    123456789);


sim_sad = (sub_data.Cue_idx == 1 & sim.y(:,1) == 1) + (sub_data.Cue_idx == 0 & sim.y(:,1) == 0);
N_sad = sum(sim_sad);


visualise_psychometric(u, sub_data, 'prc1_ehgf_binary_pu_tbt', prc_params, 'obs1_comb_obs', obs_params, 20)



%% Plot trajectory
prc1_ehgf_binary_tbt_plotTraj(sim);


%% recover parameters

est = tapas_fitModel(...
    sim.y,...
    sim.u,...
    prc_model_config,...
    obs_model_config,...
    optim_config);


%% Simulate psychometric functions

% figure('name', 'simulated psychometric'); hold on;
% prc_params_sim = prc_params;
% obs_params_sim = obs_params;
% for rho = [-2, 0, 2]
% 
%     prc_params_sim(8) = rho;
% 
%     sim = tapas_simModel(u,...
%         'prc1_ehgf_binary_pu_tbt',...
%         prc_params_sim,...
%         'obs1_comb_obs',...
%         obs_params_sim,...
%         123456789);
% 
% 
%     % Visualise Psychometric
%     sim_sad = (sub_data.Cue_idx == 1 & sim.y(:,1) == 1) + (sub_data.Cue_idx == 0 & sim.y(:,1) == 0);
%     sim_psychometric = arrayfun(@(x) mean(sim_sad(sub_data.Outcome_p_sad==x, 1)), 0:20:100);
%     plot(0:20:100, sim_psychometric, 'linewidth', 3, 'DisplayName', sprintf('\\rho = %1.1f', rho));
% end
% 
% set(gca, 'Ylim', [0,1], 'Xtick', 0:20:100)
% ylabel('p(Sad)')
% xlabel('%Sad')
% legend()

%% Full parameter recovery

% set N iterations
N = 200;

% Parameters to recover
prc_param_names = {'om2', 'rho2', 'al'};
obs_param_names = {};

recov = parameter_recovery_master(u, prc_model_config, obs_model_config, optim_config, N, prc_param_names, obs_param_names, true);
save('model1_recovery.mat', 'recov');

% recovery looks good for rho2 and omega2. For alpha, some huge outliers,
% but looks good (in log space) when they are removed

