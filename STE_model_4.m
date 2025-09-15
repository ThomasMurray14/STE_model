% model 4 (no learning)


% I think it's working, but for some reason there is no noise in the
% psychometric part


%%
clear;
close all;
addpath('prc1');
addpath('obs3');

%% 

% example data (to get contingencies etc)
sub_data = readtable('STE_data\10369536_A_Threat.csv');

% input as %sad
u = sub_data.Outcome_p_sad/100;

% responses
y = [sub_data.Response_idx, log(sub_data.Response_RT)];


%% Get configuration structures
prc_model_config = tapas_ehgf_binary_config(); % perceptual model
obs_model_config = obs3_comb_obs_config(); % response model
optim_config     = tapas_quasinewton_optim_config(); % optimisation algorithm
optim_config.nRandInit = 5;


%% Set parameters and priors
% Obs model only
obs_model_config.logitalphamu = tapas_logit(.5, 1);
obs_model_config.logitalphasa = 4;

obs_model_config.logbetamu = log(10);
obs_model_config.logbetasa = 4;

obs_model_config.logitgamma0mu = tapas_logit(.0001, 1);
obs_model_config.logitgamma0sa = 4;

obs_model_config.logitlambda0mu = tapas_logit(.0001, 1);
obs_model_config.logitlambda0sa = 4;

obs_model_config.logitmumu = tapas_logit(.5, 1);
obs_model_config.logitmusa = 4;

obs_model_config.logsigma0mu = log(.2);
obs_model_config.logsigma0sa = 4;

obs_model_config.loggamma1mu = log(6); % baseline logRT
obs_model_config.loggamma1sa = 4;

obs_model_config.loglambda1mu = log(1); % peak of logRT (above baseline)
obs_model_config.loglambda1sa = 4;

obs_model_config.logsigma1mu = log(1); % logRT noise
obs_model_config.logsigma1sa = 4;

obs_model_config = tapas_align_priors(obs_model_config);


%% simulate responses

r_temp = [];
r_temp.c_prc.n_levels = 3;
prc_params = tapas_ehgf_binary_transp(r_temp, prc_model_config.priormus);

obs_params = obs3_comb_obs_transp([], obs_model_config.priormus);

sim = tapas_simModel(u,...
    'tapas_ehgf_binary',...
    prc_params,...
    'obs3_comb_obs',...
    obs_params,...
    123456789);

%% Visualise functions
% visualise_psychometric uses contingency space
N = 20;
all_y = nan(size(u,1), 2, N);
for i = 1:N
    sim = tapas_simModel(u,...
        'tapas_ehgf_binary',...
        prc_params,...
        'obs3_comb_obs',...
        obs_params);
    all_y(:,:,i) = sim.y;
end

y_resp = squeeze(all_y(:,1,:)); % uxN
y_RT = squeeze(all_y(:,2,:));

% visualise psychometric
figure('name', 'simulated response distribution'); hold on;
all_resp_dists = zeros(N, 6);
for i = 1:N
    sim_resp_dist = arrayfun(@(x) mean(y_resp(sub_data.Outcome_p_sad==x, 1)), 0:20:100);
    all_resp_dists(i, :) = sim_resp_dist;
    plot(0:20:100, sim_resp_dist, 'linewidth', 1, 'color', [.5,.5,.5]);
end
mean_resp = mean(all_resp_dists, 1);
plot(0:20:100, mean_resp, 'linewidth', 3, 'Color', [0    0.4470    0.7410]);
set(gca, 'Ylim', [0,1], 'Xtick', 0:20:100)

% visualise RT
figure('name', 'simulated RT distributino'); hold on;
all_RT_dists = zeros(N, 6);
for i = 1:N
    sim_RT = y_RT(:, i);
    sim_RT_dist = arrayfun(@(x) mean(sim_RT(sub_data.Outcome_p_sad==x, 1)), 0:20:100);
    all_RT_dists(i, :) = sim_RT_dist;
    plot(0:20:100, sim_RT_dist, 'linewidth', 1, 'color', [.5,.5,.5]);
end
mean_RT = mean(all_RT_dists, 1);
plot(0:20:100, mean_RT, 'linewidth', 3, 'Color', [0    0.4470    0.7410]);
set(gca, 'Xtick', 0:20:100)

