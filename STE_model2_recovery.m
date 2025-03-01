% Script to run parameter recovery on model 2

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



