% Model 8:
% perceptual = eHGF (Nace remix)
% response = logRT


%%
clear;
close all;

addpath('STE_model8');

%% 

% example data (to get contingencies etc)
sub_data = readtable('STE_data\10369536_A_Threat.csv');

% Contingency space
sub_data.p_sad = sub_data.Outcome_p_sad/100;
state = double(sub_data.Cue_idx == sub_data.Outcome_idx);
u_al = sub_data.p_sad;
tr_temp = or((state ==1 & sub_data.Outcome_idx== 0), (state == 0 & sub_data.Outcome_idx == 1));
u_al(tr_temp) = 1 - sub_data.p_sad(tr_temp);

sub_data.u_al = u_al;
sub_data.state=state;

cue = sub_data.Cue_idx;
cue(cue==0) = -1; % balanced contrast coding for the "happy bias"

% Responses
sub_data.logRT = log(sub_data.Response_RT);
sub_data.correct = sub_data.Outcome_idx == sub_data.Response_idx;

% model input
u = [sub_data.u_al, cue];
y = sub_data.logRT;

u_sub = u(~isnan(sub_data.Response_idx),:);
y_sub = y(~isnan(sub_data.Response_idx),:);


%% Get configuration structures
prc_model_config = prc3_ehgf_binary_pu_tbt_config(); % perceptual model
obs_model_config = obs_logrt_linear_binary_config(); % response model


optim_config     = tapas_quasinewton_optim_config(); % optimisation algorithm
optim_config.nRandInit = 5;



%% simulate responses

prc_model_config.logalmu = log(0.1);
prc_model_config.rhomu(2) = 0;
prc_model_config.ommu(2)=-2;
prc_model_config = tapas_align_priors(prc_model_config);

r_temp = [];
r_temp.c_prc.n_levels = 3;
prc_params = prc3_ehgf_binary_pu_tbt_transp(r_temp, prc_model_config.priormus);

obs_model_config.be1mu = 0;
obs_model_config = tapas_align_priors(obs_model_config);
obs_logrt_linear_binary_transp([], obs_model_config.priormus)

sim = tapas_simModel(u_sub,...
    'prc3_ehgf_binary_pu_tbt',...
    prc_params,...
    'obs_logrt_linear_binary',...
    obs_params);
    

%% Recover
est = tapas_fitModel(...
    sim.y,...
    sim.u,...
    prc_model_config,...
    obs_model_config,...
    optim_config);

% Check parameter identifiability
tapas_fit_plotCorr(est)

