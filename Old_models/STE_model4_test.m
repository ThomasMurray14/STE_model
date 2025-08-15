% Perceptual = VKF
% response = Psychometric (prediction of PSE)


close all; clear;
addpath('custom_hgf_2');
addpath('vkf')

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

prc_model_config = prc_vkf_binary_config(); % perceptual model
obs_model_config = obs2_psychometric_config(); % response model
optim_config     = tapas_quasinewton_optim_config(); % optimisation algorithm
optim_config.nRandInit = 10;

lambda = 0.1;
v0 = .15;
omega = .1;
w0 = omega;
prc_params = [lambda, v0, omega, w0];

b0 = .5;
b1 = 0;
zeta = 10;
obs_params = [b0, b1, zeta];

sim = tapas_simModel(u_sub,...
    'prc_vkf_binary',...
    prc_params,...
    'obs2_psychometric',...
    obs_params,...
    123456789);

% visualise psychometric
% sim_psychometric = arrayfun(@(x) mean(sim.y(sub_data.p_sad==x), 'omitnan'), 0:.2:1);
% figure('name', 'simulated psychometric'); 
% plot(0:.2:1, sim_psychometric, 'linewidth', 3);
% set(gca, 'Ylim', [0,1], 'Xtick', 0:.2:1)

sim.traj.contingency = sub_data.Contingency;
plot_vkf(sim);

%% recover parameters

optim_config.nRandInit = 5;

% prc_model_config.logv0sa = 0;
% prc_model_config.logw0sa = 0;
prc_model_config = tapas_align_priors(prc_model_config);

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

% tapas_fit_plotCorr(est)
% plot_vkf(est);



%% full parameter recovery (OLD)
% 
% lambda_range = [0 1];
% omega_range = [0 8];
% b0_range = [.1 .9];
% b1_range = [-.5, .5];
% zeta_range = [.1, 40];
% 
% 
% N=300;
% [lambda_sim, lambda_est, omega_sim, omega_est, b0_sim, b0_est, b1_sim, b1_est, zeta_sim, zeta_est] = deal(nan(N,1));
% 
% for i=1:N
%     % get parameters to simulate
%     lambda  = lambda_range(1)   + (lambda_range(2)  -lambda_range(1))   *rand;
%     omega   = omega_range(1)    + (omega_range(2)   -omega_range(1))    *rand;
%     b0      = b0_range(1)       + (b0_range(2)      -b0_range(1))       *rand;
%     b1      = b1_range(1)       + (b1_range(2)      -b1_range(1))       *rand;
%     zeta    = zeta_range(1)     + (zeta_range(2)    -zeta_range(1))     *rand;
% 
%     % parameter vectors
%     v0 = .2; % fix
%     w0 = .2; % fix
%     prc_params = [lambda, v0, omega, w0];
%     obs_params = [b0, b1, zeta];
% 
%     % simulate
%     try
%     sim = tapas_simModel(u,...
%         'prc_vkf_binary',...
%         prc_params,...
%         'obs2_psychometric',...
%         obs_params);
% 
%         if ~any(isnan(sim.y))
% 
%             % fit to simulated
%             est = tapas_fitModel(...
%                 sim.y,...
%                 sim.u,...
%                 prc_model_config,...
%                 obs_model_config,...
%                 optim_config);
% 
%             % get estimated parameters
%             lambda_sim(i)   = lambda;
%             omega_sim(i)    = omega;
%             b0_sim(i) = b0;
%             b1_sim(i) = b1;
%             zeta_sim(i) = zeta;
% 
%             lambda_est(i)   = est.p_prc.lambda;
%             omega_est(i)    = est.p_prc.omega;
%             b0_est(i) = est.p_obs.b0;
%             b1_est(i) = est.p_obs.b1;
%             zeta_est(i) = est.p_obs.zeta;
%         end
%     catch
%     end
% end
% 
% 
% figure('name', 'lambda');   hold on; scatter(lambda_sim, lambda_est);       refline(1,0); xlabel('Simulated'); ylabel('Estimated'); title('Lambda');
% figure('name', 'omega');   hold on; scatter(omega_sim, omega_est);       refline(1,0); xlabel('Simulated'); ylabel('Estimated'); title('Omega');
% figure('name', 'b0');   hold on; scatter(b0_sim, b0_est);       refline(1,0); xlabel('Simulated'); ylabel('Estimated'); title('B0');
% figure('name', 'b1');   hold on; scatter(b1_sim, b1_est);       refline(1,0); xlabel('Simulated'); ylabel('Estimated'); title('B1');
% figure('name', 'zeta'); hold on; scatter(zeta_sim, zeta_est);   refline(1,0); xlabel('Simulated'); ylabel('Estimated'); title('Zeta');
% 
% save('STE_model4_recovery.mat', 'lambda_sim', 'lambda_est', 'omega_sim', 'omega_est', 'b0_sim', 'b0_est', 'b1_sim', 'b1_est', 'zeta_sim', 'zeta_est');
% 
% 
% 





