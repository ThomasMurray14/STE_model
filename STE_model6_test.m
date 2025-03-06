% Ok so Nace's model actually works if I use posteriors instead of
% predictions... Well sort of - alpha seems to have the desired effect on
% the psychometric function. However rho (sad bias) does not. So could I
% use the trial-wise mu2 (which is determined by rho) to determine
% psychometric PSE?

% I think the problem is that rho (any bias) should affect the ambiguous
% faces more.


%%
clear;
close all;

addpath('custom_hgf_3');

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
sub_data.correct = sub_data.Outcome_idx == sub_data.Response_idx; %%% Not sure why we model correct

% model input
u = [sub_data.u_al, cue];
y = [sub_data.correct,sub_data.logRT];

u_sub = u(~isnan(sub_data.Response_idx),:);
y_sub = y(~isnan(sub_data.Response_idx),:);


%% Get configuration structures
prc_model_config = prc3_ehgf_binary_pu_tbt_config(); % perceptual model
obs_model_config = tapas_unitsq_sgm_config(); % response model

obs_model_config.predorpost = 2; % use posterior

optim_config     = tapas_quasinewton_optim_config(); % optimisation algorithm
optim_config.nRandInit = 5;



%% simulate responses

prc_model_config.logalmu = log(0.1);
prc_model_config.rhomu(2) = 2;
prc_model_config.ommu(2)=-2;
prc_model_config = tapas_align_priors(prc_model_config);

r_temp = [];
r_temp.c_prc.n_levels = 3;
prc_params = prc3_ehgf_binary_pu_tbt_transp(r_temp, prc_model_config.priormus);

obs_params = 5; % zeta (inverse decision noise)


figure('name', 'simulated psychometric'); hold on;

for i=1:20
sim = tapas_simModel(u_sub,...
    'prc3_ehgf_binary_pu_tbt',...
    prc_params,...
    'tapas_unitsq_sgm',...
    obs_params);

    sim_sad = (sub_data.Cue_idx == 1 & sim.y(:,1) == 1) + (sub_data.Cue_idx == 0 & sim.y(:,1) == 0);
    
    sim_psychometric = arrayfun(@(x) mean(sim_sad(sub_data.Outcome_p_sad==x, 1)), 0:20:100);
    plot(0:20:100, sim_psychometric, 'linewidth', 1);
end

% Visualise Psychometric
% y is "correct"

set(gca, 'Ylim', [0,1], 'Xtick', 0:20:100)







