%% TO DO

% Change the response model:
% try merging correct/incorrect + logrt + confidence

% can we add mu3 to explain RTs?
% split corr/incorr or 
% tapas_ehgf_binary_combObs_plotTraj
% bias of happy vs sad (perception or response)

%%
clear;
close all;

addpath('custom_hgf');

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
y = [sub_data.correct,sub_data.logRT];

u_sub = u(~isnan(sub_data.Response_idx),:);
y_sub = y(~isnan(sub_data.Response_idx),:);


%% Get configuration structures
prc_model_config = prc1_ehgf_binary_pu_tbt_config(); % perceptual model
obs_model_config = obs1_comb_obs_config(); % response model
optim_config     = tapas_quasinewton_optim_config(); % optimisation algorithm
optim_config.nRandInit = 5;


%% simulate responses

prc_model_config.logalmu = log(0.1);
prc_model_config.rhomu(2) = 0;
prc_model_config.ommu(2)=-2;
prc_model_config = tapas_align_priors(prc_model_config);

r_temp = [];
r_temp.c_prc.n_levels = 3;
prc_params = prc1_ehgf_binary_pu_tbt_transp(r_temp, prc_model_config.priormus);

obs_params = obs_model_config.priormus;
obs_params(1) = exp(obs_params(1));
obs_params(7) = exp(obs_params(7));

sim = tapas_simModel(u_sub,...
    'prc1_ehgf_binary_pu_tbt',...
    prc_params,...
    'obs1_comb_obs',...
    obs_params,...
    123456789);


% Visualise Psychometric
% y is "correct"
% sim_sad = (sub_data.Cue_idx == 1 & sim.y(:,1) == 1) + (sub_data.Cue_idx == 0 & sim.y(:,1) == 0);
sim_sad = (sub_data.p_sad>.5 & sim.y(:,1)==1) + (sub_data.p_sad<.5 & sim.y(:,1)==0);

figure('name', 'simulated psychometric'); hold on;
sim_psychometric = arrayfun(@(x) mean(sim_sad(sub_data.Outcome_p_sad==x, 1)), 0:20:100);
plot(0:20:100, sim_psychometric, 'linewidth', 3);
set(gca, 'Ylim', [0,1], 'Xtick', 0:20:100)



%% Plot trajectory
% prc1_ehgf_binary_tbt_plotTraj(sim);


%% recover parameters

prc_model_config = prc1_ehgf_binary_pu_tbt_config(); % perceptual model

prc_model_config.omsa(2)=3;
prc_model_config = tapas_align_priors(prc_model_config);

obs_model_config.logzesa=2;
obs_model_config.beta0sa=2;
obs_model_config.beta1sa=2;
obs_model_config.beta2sa=2;
obs_model_config.beta3sa=2;
obs_model_config.beta4sa=0;
obs_model_config.logsasa=2;
obs_model_config = tapas_align_priors(obs_model_config);


est = tapas_fitModel(...
    sim.y,...
    sim.u,...
    prc_model_config,...
    obs_model_config,...
    optim_config);

% Check parameter identifiability
% tapas_fit_plotCorr(est)
% prc1_ehgf_binary_tbt_plotTraj(est)



%%
% %% full parameter recovery (OLD)
% 
% N=300;
% 
% rho2_range = [-2 2];
% om2_range = [-5 -1];
% al_range = [.005 2]; % not sure about ranges...
% 
% om2_sim=nan(N,1);
% om2_est=nan(N,1);
% rho2_sim=nan(N,1);
% rho2_est=nan(N,1);
% al_sim=nan(N,1);
% al_est=nan(N,1);
% 
% 
% for i=1:N
%     % get parameters to simulate
%     om2 = om2_range(1) + (om2_range(2)-om2_range(1))*rand;
%     rho2 = rho2_range(1) + (rho2_range(2)-rho2_range(1))*rand;
%     al = al_range(1) + (al_range(2)-al_range(1))*rand;
% 
%     % parameter vectors
%     prc_params(8) = rho2;
%     prc_params(13) = om2;
%     prc_params(15) = al; % specify in native space
% 
%     % simulate
%     sim = tapas_simModel(u,...
%         'prc1_ehgf_binary_pu_tbt',...
%         prc_params,...
%         'obs1_comb_obs',...
%         obs_params);
% 
%     if ~any(isnan(sim.y))
%         try
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
%             rho2_sim(i) = rho2;
%             al_sim(i) = al;
%             om2_est(i) = est.p_prc.om(2);
%             rho2_est(i) = est.p_prc.rho(2);
%             al_est(i) = est.p_prc.al;
%         catch
% 
%         end
%     end
% end
% 
% 
% figure('name', 'om2');  hold on; scatter(om2_sim, om2_est);     refline(1,0); xlabel('Simulated'); ylabel('Estimated'); title('Omega2');
% figure('name', 'rho2'); hold on; scatter(rho2_sim, rho2_est);   refline(1,0); xlabel('Simulated'); ylabel('Estimated'); title('Rho');
% figure('name', 'al');   hold on; scatter(al_sim, al_est);       refline(1,0); xlabel('Simulated'); ylabel('Estimated'); title('Alpha');
% 
% save('STE_model1_recovery.mat', 'om2_sim', 'om2_est', 'rho2_sim', 'rho2_est', 'al_sim', 'al_est');


%% Simulate psychometric functions
% config structures
prc_model_config = prc1_ehgf_binary_pu_tbt_config(); % perceptual model
obs_model_config = obs1_comb_obs_config(); % response model
optim_config     = tapas_quasinewton_optim_config(); % optimisation algorithm
optim_config.nRandInit = 5;


figure('name', 'simulated psychometric'); hold on;

for rho = [-5, 0, 5]

    prc_model_config.logalmu = log(.05);
    prc_model_config.rhomu(2) = rho;
    prc_model_config.ommu(2)=-2;
    prc_model_config = tapas_align_priors(prc_model_config);

    r_temp = [];
    r_temp.c_prc.n_levels = 3;
    prc_params = prc1_ehgf_binary_pu_tbt_transp(r_temp, prc_model_config.priormus);

    obs_params = obs_model_config.priormus;
    obs_params(1) = exp(obs_params(1));
    obs_params(7) = exp(obs_params(7));

    sim = tapas_simModel(u_sub,...
        'prc1_ehgf_binary_pu_tbt',...
        prc_params,...
        'obs1_comb_obs',...
        obs_params,...
        123456789);


    % Visualise Psychometric
    % y is "correct"
    sim_sad = (sub_data.p_sad>.5 & sim.y(:,1)==1) + (sub_data.p_sad<.5 & sim.y(:,1)==0);

    sim_psychometric = arrayfun(@(x) mean(sim_sad(sub_data.Outcome_p_sad==x, 1)), 0:20:100);
    plot(0:20:100, sim_psychometric, 'linewidth', 3, 'DisplayName', sprintf('\\rho = %1.1f', rho));
end

set(gca, 'Ylim', [0,1], 'Xtick', 0:20:100)
ylabel('p(Sad)')
xlabel('%Sad')
legend()


