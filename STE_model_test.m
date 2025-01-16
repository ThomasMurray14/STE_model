
% close all; clear;
addpath('custom_hgf_2');

% example data (to get contingencies etc)
sub_data = readtable('10369536_A_Safe.csv');

% input in contingency space
state = double(sub_data.Cue_idx == sub_data.Outcome_idx);

% stimulus intensity
sub_data.p_sad = sub_data.Outcome_p_sad/100;

% model input
u = [state, sub_data.p_sad];


u = state; % just contingency space




%% simulate responses

prc_model_config = tapas_ehgf_binary_config(); % perceptual model
obs_model_config = tapas_unitsq_sgm_config(); % response model
optim_config     = tapas_quasinewton_optim_config(); % optimisation algorithm



%%


prc_model_config.ommu(2) = -2.2;
prc_model_config.priormus(13)=prc_model_config.ommu(2);

prc_model_config.ommu(3) = 0.5;
prc_model_config.priormus(14)=prc_model_config.ommu(3);

prc_model_config.logkamu(2) = log(.1);
prc_model_config.priormus(11)=prc_model_config.logkamu(2);


r_temp = [];
r_temp.c_prc.n_levels = 3;
prc_params = tapas_ehgf_binary_transp(r_temp, prc_model_config.priormus);

obs_params = obs_model_config.priormus;

sim = tapas_simModel(u,...
    'tapas_ehgf_binary',...
    prc_params,...
    'tapas_unitsq_sgm',...
    obs_params,...
    123456789);


%% recover parameters
prc_model_config = tapas_ehgf_binary_config(); % perceptual model
obs_model_config = tapas_unitsq_sgm_config(); % response model

optim_config.nRandInit = 5;

prc_model_config.ommu(2) = -3;
prc_model_config.priormus(13)=prc_model_config.ommu(2);

prc_model_config.omsa(2) = 16;
prc_model_config.priorsas(13)=prc_model_config.omsa(2);

prc_model_config.ommu(3) = 0;
prc_model_config.priormus(14)=prc_model_config.ommu(3);

prc_model_config.omsa(3) = 0;
prc_model_config.priorsas(14)=prc_model_config.omsa(3);

prc_model_config.logkasa(1) = 0;
prc_model_config.priorsas(10) = prc_model_config.logkasa(1);

prc_model_config.logkasa(2) = 4;
prc_model_config.priorsas(11) = prc_model_config.logkasa(2);


est = tapas_fitModel(...
    sim.y,...
    sim.u,...
    prc_model_config,...
    obs_model_config,...
    optim_config);

% tapas_fit_plotCorr(est)
tapas_hgf_binary_plotTraj(est)


mu2 = est.traj.mu(:,2);
mu3 = est.traj.mu(:,3);
pv = tapas_sgm(mu2, 1).*(1-tapas_sgm(mu2, 1)).*exp(mu3);

figure; plot(1:192, exp(mu3) )









