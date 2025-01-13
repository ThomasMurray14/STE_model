% Perceptual = binary HGF; response = binary (Tom's psychometric VIP)

% prc model is normal binary HGF, but edited to ignore column 2 of u


close all; clear;

% example data (to get contingencies etc)
sub_data = readtable('STE_data\10369536_A_Threat.csv');

% input in contingency space
state = double(sub_data.Cue_idx == sub_data.Outcome_idx);

% stimulus intensity
sub_data.p_sad = sub_data.Outcome_p_sad/100;

% model input
u = [state, sub_data.p_sad];

% u = u(:,1)

%% Models


%% simulate responses

prc_model_config = prc2_hgf_binary_config(); % perceptual model
obs_model_config = obs2_psychometric_config(); % response model
% obs_model_config = tapas_unitsq_sgm_config(); % response model
optim_config     = tapas_quasinewton_optim_config(); % optimisation algorithm
optim_config.nRandInit = 5;

r_temp = [];
r_temp.c_prc.n_levels = 3;
prc_params = prc2_hgf_binary_transp(r_temp, prc_model_config.priormus);
obs_params = obs_model_config.priormus;


sim = tapas_simModel(u,...
    'prc2_hgf_binary',...
    prc_params,...
    'obs2_psychometric',...
    obs_params,...
    123456789);

% visualise psychometric
sim_psychometric = arrayfun(@(x) mean(sim.y(sub_data.p_sad==x, 1)), 0:.2:1);
% figure('name', 'simulated psychometric'); hold on;
plot(0:.2:1, sim_psychometric, 'linewidth', 3);
set(gca, 'Ylim', [0,1], 'Xtick', 0:.2:1)


%% recover parameters

prc_model_config = prc2_hgf_binary_config(); % perceptual model
obs_model_config = obs2_psychometric_config(); % response model
% obs_model_config = tapas_unitsq_sgm_config(); % response model

optim_config.nRandInit = 50;

prc_model_config.omsa(3)=0;
prc_model_config.priorsas(14)=prc_model_config.omsa(3);

obs_model_config.b0sa=1;
obs_model_config.b1sa=1;
obs_model_config.priorsas(1)=obs_model_config.b0sa;
obs_model_config.priorsas(2)=obs_model_config.b1sa;

est = tapas_fitModel(...
    sim.y,...
    sim.u,...
    prc_model_config,...
    obs_model_config,...
    optim_config);







