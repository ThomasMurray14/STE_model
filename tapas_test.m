% script to test out models in the tapas toolbox


% example data (to get contingencies etc)
sub_data = readtable('STE_data\10369536_A_Threat.csv');

% model input
u = sub_data.Outcome_idx == sub_data.Cue_idx;
y = log(sub_data.Response_RT);

% mode configs
prc_model_config = tapas_hgf_binary_config(); % perceptual model
obs_model_config = tapas_logrt_linear_binary_config(); % response model
optim_config     = tapas_quasinewton_optim_config(); % optimisation algorithm

% simulate responses
bopars = tapas_fitModel([],...
                         u,...
                         'tapas_hgf_binary_config',...
                         'tapas_bayes_optimal_binary_config',...
                         'tapas_quasinewton_optim_config');

%%
prc_params = bopars.p_prc.p;
obs_params = tapas_logrt_linear_binary_transp(obs_model_config, obs_model_config.priormus);
sim = tapas_simModel(u, 'tapas_hgf_binary', prc_params, 'tapas_logrt_linear_binary', obs_params, 123456789);


tapas_hgf_binary_plotTraj(sim);