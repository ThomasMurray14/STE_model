% Perceptual = binary HGF; 
% response = combined psychometric (prediction of PSE) / logRT

% prc model is normal binary HGF, but edited to ignore column 2 of u


% All the parameters are called beta... b0, b1, and zeta are for the
% psychometric; beta0, beta1, etc are for the logRT


close all; clear;
addpath('custom_hgf_2');

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
obs_model_config = obs2_comb_obs_config(); % response model
optim_config     = tapas_quasinewton_optim_config(); % optimisation algorithm
optim_config.nRandInit = 10;


r_temp = [];
r_temp.c_prc.n_levels = 3;
prc_params = prc2_ehgf_binary_transp(r_temp, prc_model_config.priormus);

% obs_model_config.b1mu = .2;
obs_params = obs2_comb_obs_transp([], obs_model_config.priormus);

sim = tapas_simModel(u_sub,...
    'prc2_ehgf_binary',...
    prc_params,...
    'obs2_comb_obs',...
    obs_params,...
    123456789);


%% recover parameters

optim_config.nRandInit = 5;

% estimate regression coefficients
obs_model_config.b0sa=2;
obs_model_config.b1sa=2;
obs_model_config = tapas_align_priors(obs_model_config);

est = tapas_fitModel(...
    sim.y,...
    sim.u,...
    prc_model_config,...
    obs_model_config,...
    optim_config);

tapas_fit_plotCorr(est)
tapas_hgf_binary_plotTraj(est)



