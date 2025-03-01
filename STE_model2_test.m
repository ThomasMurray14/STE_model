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



%% full parameter recovery

close all;

% reload model configs
prc_model_config = prc2_ehgf_binary_config(); % perceptual model
obs_model_config = obs2_psychometric_config(); % response model
optim_config     = tapas_quasinewton_optim_config(); % optimisation algorithm
optim_config.nRandInit = 10;


% set priors
prc_model_config.omsa(2) = 4;
prc_model_config.omsa(3) = 4;
prc_model_config = tapas_align_priors(prc_model_config);

obs_model_config.b0sa=2;
obs_model_config.b1sa=2;
obs_model_config.logzetasa=4;

obs_model_config = tapas_align_priors(obs_model_config);

% set N iterations
N = 300;

% preallocate sim and est parameters
[om2_sim, om2_est, om3_sim, om3_est, b0_sim, b0_est, b1_sim, b1_est, zeta_sim, zeta_est] = deal(nan(N,1));

for i = 1:N
    sim = tapas_sampleModel(u, prc_model_config, obs_model_config);
    om2_sim(i)  = sim.p_prc.om(2);
    om3_sim(i)  = sim.p_prc.om(3);
    b0_sim(i)   = sim.p_obs.b0;
    b1_sim(i)   = sim.p_obs.b1;
    zeta_sim(i) = sim.p_obs.zeta;

    if ~any(isnan(sim.y)) % only no missing trials
        try
            est = tapas_fitModel(...
                sim.y,...
                sim.u,...
                prc_model_config,...
                obs_model_config,...
                optim_config);
        
            % get estimated parameters
            if ~isinf(est.optim.LME)
                om2_est(i)  = est.p_prc.om(2);
                om3_est(i)  = est.p_prc.om(3);
                b0_est(i)   = est.p_obs.b0;
                b1_est(i)   = est.p_obs.b1;
                zeta_est(i) = est.p_obs.zeta;
            end
        catch
            fprintf('\nCannot fit')
        end
    end
end

% figure('name', 'om2');  hold on; scatter(om2_sim, om2_est);     refline(1,0); xlabel('Simulated'); ylabel('Estimated'); title('Omega2');
% figure('name', 'om3');  hold on; scatter(om3_sim, om3_est);     refline(1,0); xlabel('Simulated'); ylabel('Estimated'); title('Omega3');
% figure('name', 'b0');   hold on; scatter(b0_sim, b0_est);       refline(1,0); xlabel('Simulated'); ylabel('Estimated'); title('B0');
% figure('name', 'b1');   hold on; scatter(b1_sim, b1_est);       refline(1,0); xlabel('Simulated'); ylabel('Estimated'); title('B1');
% figure('name', 'zeta'); hold on; scatter(zeta_sim, zeta_est);   refline(1,0); xlabel('Simulated'); ylabel('Estimated'); title('Zeta');

save('STE_model2_recovery.mat', 'om2_sim', 'om2_est', 'om3_sim', 'om3_est', 'b0_sim', 'b0_est', 'b1_sim', 'b1_est', 'zeta_sim', 'zeta_est');




%% full parameter recovery (OLD)
% 
% om2_range = [-5 1];
% om3_range = [-5 1];
% b0_range = [.1 .9];
% b1_range = [-.5, .5];
% zeta_range = [.1, 40];
% 
% N=300;
% [om2_sim, om2_est, om3_sim, om3_est, b0_sim, b0_est, b1_sim, b1_est, zeta_sim, zeta_est] = deal(nan(N,1));
% 
% 
% for i=1:N
%     % get parameters to simulate
%     om2 = om2_range(1) + (om2_range(2)-om2_range(1))*rand;
%     om3 = om3_range(1) + (om3_range(2)-om3_range(1))*rand;
%     b0 = b0_range(1) + (b0_range(2)-b0_range(1))*rand;
%     b1 = b1_range(1) + (b1_range(2)-b1_range(1))*rand;
%     zeta = zeta_range(1) + (zeta_range(2)-zeta_range(1))*rand;
% 
%     % parameter vectors
%     prc_params(13) = om2;
%     prc_params(14) = om3;
%     obs_params = [b0, b1, zeta];
% 
%     % simulate
%     try
%     sim = tapas_simModel(u,...
%         'prc2_ehgf_binary',...
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
%             om2_sim(i) = om2;
%             om3_sim(i) = om3;
%             b0_sim(i) = b0;
%             b1_sim(i) = b1;
%             zeta_sim(i) = zeta;
%             om2_est(i) = est.p_prc.om(2);
%             om3_est(i) = est.p_prc.om(3);
%             b0_est(i) = est.p_obs.b0;
%             b1_est(i) = est.p_obs.b1;
%             zeta_est(i) = est.p_obs.zeta;
%         end
%     catch
%     end
% end
% 
% 
% 
% figure('name', 'om2');  hold on; scatter(om2_sim, om2_est);     refline(1,0); xlabel('Simulated'); ylabel('Estimated'); title('Omega2');
% figure('name', 'om3');  hold on; scatter(om3_sim, om3_est);     refline(1,0); xlabel('Simulated'); ylabel('Estimated'); title('Omega3');
% figure('name', 'b0');   hold on; scatter(b0_sim, b0_est);       refline(1,0); xlabel('Simulated'); ylabel('Estimated'); title('B0');
% figure('name', 'b1');   hold on; scatter(b1_sim, b1_est);       refline(1,0); xlabel('Simulated'); ylabel('Estimated'); title('B1');
% figure('name', 'zeta'); hold on; scatter(zeta_sim, zeta_est);   refline(1,0); xlabel('Simulated'); ylabel('Estimated'); title('Zeta');
% 
% 
% save('STE_model2_recovery.mat', 'om2_sim', 'om2_est', 'om3_sim', 'om3_est', 'b0_sim', 'b0_est', 'b1_sim', 'b1_est', 'zeta_sim', 'zeta_est');


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

