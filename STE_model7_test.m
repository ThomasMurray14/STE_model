% Psychometric-pse, but using surprise trajectory and 2 level HGF


close all; clear;
addpath('custom_hgf_4');

% example data (to get contingencies etc)
sub_data = readtable('STE_data\10369536_A_Threat.csv');

% input in contingency space
state = double(sub_data.Cue_idx == sub_data.Outcome_idx);

% stimulus intensity
sub_data.p_sad = sub_data.Outcome_p_sad/100;

% model input
u = [state, sub_data.p_sad];
u_sub = u;




%% simulate responses

prc_model_config = prc2_ehgf_binary_config(); % perceptual model
obs_model_config = obs4_psychometric_config(); % response model
optim_config     = tapas_quasinewton_optim_config(); % optimisation algorithm
optim_config.nRandInit = 10;


prc_model_config.ommu(3) = 0;
prc_model_config.omsa(3) = 0;
prc_model_config.logkamu(2) = -inf;

prc_model_config = tapas_align_priors(prc_model_config);
r_temp = [];
r_temp.c_prc.n_levels = 3;
prc_params = prc2_ehgf_binary_transp(r_temp, prc_model_config.priormus);


obs_model_config.logitb0mu = tapas_logit(.5, 1);
obs_model_config.logitzetamu = tapas_logit(.2, 1);
obs_model_config = tapas_align_priors(obs_model_config);
obs_params = obs4_psychometric_transp([], obs_model_config.priormus);


sim = tapas_simModel(u_sub,...
    'prc2_ehgf_binary',...
    prc_params,...
    'obs4_psychometric',...
    obs_params,...
    123456789);

% visualise psychometric
sim_psychometric = arrayfun(@(x) mean(sim.y(sub_data.p_sad==x), 'omitnan'), 0:.2:1);
% figure('name', 'simulated psychometric'); hold on;
plot(0:.2:1, sim_psychometric, 'linewidth', 3);
set(gca, 'Ylim', [0,1], 'Xtick', 0:.2:1)


%% Recover
est = tapas_fitModel(...
    sim.y,...
    sim.u,...
    prc_model_config,...
    obs_model_config,...
    optim_config);

% Check parameter identifiability
tapas_fit_plotCorr(est)
prc1_ehgf_binary_tbt_plotTraj(est)
