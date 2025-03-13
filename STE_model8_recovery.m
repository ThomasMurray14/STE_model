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


%% Load configs and set priors
prc_model_config = prc3_ehgf_binary_pu_tbt_config(); % perceptual model
obs_model_config = obs_logrt_linear_binary_config(); % response model

optim_config     = tapas_quasinewton_optim_config(); % optimisation algorithm
optim_config.nRandInit = 5;

% SEt priors
prc_model_config.logalmu = log(0.1);
prc_model_config.logalsa = 2;
prc_model_config.rhomu(2) = 0;
prc_model_config.rhosa(2) = 2;
prc_model_config.ommu(2) = -2;
prc_model_config.omsa(2) = 2;
prc_model_config = tapas_align_priors(prc_model_config);

r_temp = [];
r_temp.c_prc.n_levels = 3;
prc_params = prc3_ehgf_binary_pu_tbt_transp(r_temp, prc_model_config.priormus);

obs_model_config.be0sa = 2;
obs_model_config.be1sa = 2;
obs_model_config.be2sa = 2;
obs_model_config.be3sa = 2;
obs_model_config = tapas_align_priors(obs_model_config);
obs_logrt_linear_binary_transp([], obs_model_config.priormus);

%% Run recovery
% set N iterations
N = 300;

% preallocate sim and est parameters
[om2_sim, om2_est, ...
    rho2_sim, rho2_est, ...
    al_sim, al_est, ...
    b0_sim, b0_est, ...
    b1_sim, b1_est, ...
    b2_sim, b2_est, ...
    b3_sim, b3_est, ...
    ze_sim, ze_est] = deal(nan(N,1));

for i = 1:N
    sim = tapas_sampleModel(u, prc_model_config, obs_model_config);
    om2_sim(i)  = sim.p_prc.om(2);
    rho2_sim(i) = sim.p_prc.rho(2);
    al_sim(i)   = sim.p_prc.al;
    b0_sim(i)   = sim.p_obs.be0;
    b1_sim(i)   = sim.p_obs.be1;
    b2_sim(i)   = sim.p_obs.be2;
    b3_sim(i)   = sim.p_obs.be3;
    ze_sim(i)   = sim.p_obs.ze;

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
                rho2_est(i) = est.p_prc.rho(2);
                al_est(i)   = est.p_prc.al;
                b0_est(i)       = est.p_obs.be0;
                b1_est(i)       = est.p_obs.be1;
                b2_est(i)       = est.p_obs.be2;
                b3_est(i)       = est.p_obs.be3;
                ze_est(i)       = est.p_obs.ze;
            end
        catch
            fprintf('\nCannot fit')
        end
    end
end

figure('name', 'om2');  hold on; scatter(om2_sim, om2_est);   refline(1,0); xlabel('Simulated'); ylabel('Estimated'); title('Omega2');
figure('name', 'rho2'); hold on; scatter(rho2_sim, rho2_est); refline(1,0); xlabel('Simulated'); ylabel('Estimated'); title('Rho2');
figure('name', 'al');   hold on; scatter(al_sim, al_est);     refline(1,0); xlabel('Simulated'); ylabel('Estimated'); title('Alpha');
figure('name', 'b0');   hold on; scatter(b0_sim, b0_est);     refline(1,0); xlabel('Simulated'); ylabel('Estimated'); title('B0');
figure('name', 'b1');   hold on; scatter(b1_sim, b1_est);     refline(1,0); xlabel('Simulated'); ylabel('Estimated'); title('B1');
figure('name', 'b2');   hold on; scatter(b2_sim, b2_est);     refline(1,0); xlabel('Simulated'); ylabel('Estimated'); title('B2');
figure('name', 'b3');   hold on; scatter(b3_sim, b3_est);     refline(1,0); xlabel('Simulated'); ylabel('Estimated'); title('B3');
figure('name', 'ze');   hold on; scatter(ze_sim, ze_est);     refline(1,0); xlabel('Simulated'); ylabel('Estimated'); title('Zeta');


save('STE_model8_recovery.mat', ...
    'om2_sim', 'om2_est', ...
    'rho2_sim', 'rho2_est', ...
    'al_sim', 'al_est', ...
    'b0_sim', 'b0_est',...
    'b1_sim', 'b1_est',...
    'b2_sim', 'b1_est',...
    'b3_sim', 'b3_est',...
    'ze_sim', 'ze_est');


