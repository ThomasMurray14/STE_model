%% MODEL 2
% Perceptual model = 2 level eHGF Nace remix (Tom edit)
% Response model = combined unitsq_sgm and logRT
% 
% In the perceptual model, I have changed the way the happy bias (rho2) and
% perceptual uncertainty (al) are coded. I have kept the parameter names
% for consistency
% 
% TODO
% 
% - Figure out how to deal with al and rho - rho needs to be 0-1, al needs
% to be positive. Use logit with rho, exp/log with al
%   Specify rho in logit space (e.g. tapas_logit(rho, 1), where 0<rho<1)
%   and al in log space (e.g. log(al)) in config file and priors. Then in
%   transp, use tapas_sgm on rho (e.g. tapas_sgm(logitrho, 1)) and exp on 
%   alpha (e.g. exp(al))
% - And I am going to change rho to just be a single parameter (i.e. not
% rho1 and rho2), so make sure idx of parameters are adjusted in
% parameter_recovery
% 
% Notes:
% 




%%
clear;
close all;
addpath('STE_model2');

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
obs_model_config = obs2_comb_obs_config(); % response model
optim_config     = tapas_quasinewton_optim_config(); % optimisation algorithm
optim_config.nRandInit = 5;

%% Set parameters and priors
prc_model_config.ommu(2)    = -2;
prc_model_config.omsa(2)    = 4;

prc_model_config.logitrhomu   = tapas_logit(0.9, 1); % Discrimination threshold
prc_model_config.logitrhosa   = 2;

prc_model_config.logalmu    = log(1); % perceptual uncertainty
prc_model_config.logalsa    = 2;

prc_model_config = tapas_align_priors(prc_model_config);

obs_model_config.logzemu = log(.1);
obs_model_config.logzesa = 2;

obs_model_config.beta0mu = 6.5000;
obs_model_config.beta0sa = 4;

obs_model_config.beta1mu = 0;
obs_model_config.beta1sa = 4;

obs_model_config.beta2mu = 0;
obs_model_config.beta2sa = 4;

obs_model_config.beta3mu = 0;
obs_model_config.beta3sa = 4;

obs_model_config.beta4mu = 0;
obs_model_config.beta4sa = 4;

obs_model_config.logsasa = log(.1);
obs_model_config.logsasa = 2;

obs_model_config = tapas_align_priors(obs_model_config);

%% simulate responses

r_temp = [];
r_temp.c_prc.n_levels = 3;
prc_params = prc2_ehgf_binary_pu_tbt_transp(r_temp, prc_model_config.priormus);

obs_params = obs_model_config.priormus;
obs_params(1) = exp(obs_params(1));
obs_params(7) = exp(obs_params(7));
% obs_params(8) = 3;%exp(obs_params(8));
sim = tapas_simModel(u,...
    'prc2_ehgf_binary_pu_tbt',...
    prc_params,...
    'obs2_comb_obs',...
    obs_params,...
    123456789);


sim_sad = (sub_data.Cue_idx == 1 & sim.y(:,1) == 1) + (sub_data.Cue_idx == 0 & sim.y(:,1) == 0);
N_sad = sum(sim_sad);


visualise_psychometric(u, sub_data, 'prc2_ehgf_binary_pu_tbt', prc_params, 'obs2_comb_obs', obs_params, 20)



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
% 
% % set N iterations
% N = 200;
% 
% % Parameters to recover
% prc_param_names = {'om2', 'rho2', 'al'};
% prc_param_idx   = [13, 8, 15];
% obs_param_names = {'ze', 'beta0', 'beta1', 'beta2', 'beta3', 'beta4', 'sa'};
% obs_param_idx   = [1, 2, 3, 4, 5, 6, 7];
% 
% recov = parameter_recovery_master(u,...
%     prc_model_config,...
%     obs_model_config,...
%     optim_config,...
%     N,...
%     prc_param_names,...
%     prc_param_idx,...
%     obs_param_names,...
%     obs_param_idx,...
%     true);
% save('model2_recovery.mat', 'recov');


%% Fit actual data

% model_fits = fit_master(u, prc_model_config, obs_model_config, optim_config);
% save('model2_fit.mat', 'model_fits');





