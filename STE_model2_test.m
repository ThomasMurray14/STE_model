% Perceptual = binary HGF; 
% response = psychometric (prediction of PSE)

% prc model is normal binary HGF, but edited to ignore column 2 of u


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
obs_model_config = obs2_psychometric_config(); % response model
optim_config     = tapas_quasinewton_optim_config(); % optimisation algorithm
optim_config.nRandInit = 10;


% prc_model_config.ommu(3) = 0;
% prc_model_config.omsa(3) = 0;
% prc_model_config.logkamu(2) = -inf;

prc_model_config = tapas_align_priors(prc_model_config);
r_temp = [];
r_temp.c_prc.n_levels = 3;
prc_params = prc2_ehgf_binary_transp(r_temp, prc_model_config.priormus);


% obs_model_config.b1mu = .2;
obs_model_config = tapas_align_priors(obs_model_config);
obs_params = obs2_psychometric_transp([], obs_model_config.priormus);


sim = tapas_simModel(u_sub,...
    'prc2_ehgf_binary',...
    prc_params,...
    'obs2_psychometric',...
    obs_params,...
    123456789);

% visualise psychometric
sim_psychometric = arrayfun(@(x) mean(sim.y(sub_data.p_sad==x), 'omitnan'), 0:.2:1);
% figure('name', 'simulated psychometric'); hold on;
plot(0:.2:1, sim_psychometric, 'linewidth', 3);
set(gca, 'Ylim', [0,1], 'Xtick', 0:.2:1)


%% recover parameters

optim_config.nRandInit = 5;

% prc_model_config.ommu(3) = 0;
% prc_model_config.omsa(3)=0;
prc_model_config = tapas_align_priors(prc_model_config);

% estimate regression coefficients
obs_model_config.b0sa=2;
obs_model_config.b1sa=0;
obs_model_config = tapas_align_priors(obs_model_config);

est = tapas_fitModel(...
    sim.y,...
    sim.u,...
    prc_model_config,...
    obs_model_config,...
    optim_config);

% tapas_fit_plotCorr(est)
% tapas_hgf_binary_plotTraj(est)


%% Psychometric functions

% prc_model_config = prc2_ehgf_binary_config(); % perceptual model
% obs_model_config = obs2_psychometric_config(); % response model
% optim_config     = tapas_quasinewton_optim_config(); % optimisation algorithm
% optim_config.nRandInit = 10;
% 
% 
% figure; hold on;
% 
% for zeta = [1, 10, 50]
% 
%     prc_model_config.ommu(2) = -3;
%     prc_model_config.ommu(3) = -2;
%     prc_model_config = tapas_align_priors(prc_model_config);
%     r_temp = [];
%     r_temp.c_prc.n_levels = 3;
%     prc_params = prc2_ehgf_binary_transp(r_temp, prc_model_config.priormus);
% 
%     b1 = 0;
%     b0 = .5;
%     obs_params = [b0, b1, zeta];
% 
%     sim = tapas_simModel(u_sub,...
%         'prc2_ehgf_binary',...
%         prc_params,...
%         'obs2_psychometric',...
%         obs_params,...
%         123456789);
% 
% 
%     % visualise psychometric
%     sim_psychometric = arrayfun(@(x) mean(sim.y(sub_data.p_sad==x), 'omitnan'), 0:.2:1);
%     % figure('name', 'simulated psychometric'); hold on;
%     plot(0:.2:1, sim_psychometric, 'linewidth', 3, 'DisplayName', sprintf('\\zeta = %1.1f', zeta));
% 
% end
% 
% set(gca, 'Ylim', [0,1], 'Xtick', 0:.2:1)
% xlabel('%Sad');
% ylabel('p(Sad)')
% legend()

