% Model 6 recovery script
% 
% Perceptual = eHGF (Nace remix)
% Response = sigmoid (using posterior)
% 

%%
clear;
close all;

addpath('custom_hgf_3');

%% Load model inputs

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

% model input
u = [sub_data.u_al, cue];




%% full parameter recovery
close all;

% reload model configs
prc_model_config = prc3_ehgf_binary_pu_tbt_config(); % perceptual model
obs_model_config = tapas_unitsq_sgm_config(); % response model
obs_model_config.predorpost = 2; % use posterior
optim_config     = tapas_quasinewton_optim_config(); % optimisation algorithm
optim_config.nRandInit = 5;


% set priors
prc_model_config.logalmu = log(0.1);
prc_model_config.logalsa = 2;
prc_model_config.rhomu(2) = 0;
prc_model_config.rhosa(2) = 2;
prc_model_config.ommu(2)=-2;
prc_model_config.omsa(2)=16;
prc_model_config = tapas_align_priors(prc_model_config);

obs_model_config.logzemu = log(5);
obs_model_config.logzesa = 2;
obs_model_config = tapas_align_priors(obs_model_config);


N=300;

% preallocate sim and est parameters
[om2_sim, om2_est, rho2_sim, rho2_est, al_sim, al_est, ze_sim, ze_est] = deal(nan(N,1));

for i = 1:N
    sim = tapas_sampleModel(u, prc_model_config, obs_model_config);
    om2_sim(i)  = sim.p_prc.om(2);
    rho2_sim(i) = sim.p_prc.rho(2);
    al_sim(i)   = sim.p_prc.al;
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
                ze_est(i)   = est.p_obs.ze;
            end
        catch
            fprintf('\nCannot fit')
        end
    end
end


save('STE_model6_recovery.mat', 'om2_sim', 'om2_est', 'rho2_sim', 'rho2_est', 'al_sim', 'al_est', 'ze_sim', 'ze_est');


figure('name', 'om2');  hold on; scatter(om2_sim, om2_est);     refline(1,0); xlabel('Simulated'); ylabel('Estimated'); title('Omega2');
figure('name', 'rho2'); hold on; scatter(rho2_sim, rho2_est);   refline(1,0); xlabel('Simulated'); ylabel('Estimated'); title('Rho');
figure('name', 'al');   hold on; scatter(log(al_sim), log(al_est));       refline(1,0); xlabel('Simulated'); ylabel('Estimated'); title('Alpha');
figure('name', 'al');   hold on; scatter(log(ze_sim), log(ze_est));       refline(1,0); xlabel('Simulated'); ylabel('Estimated'); title('Zeta');

% log alpha minimum at -3
% omega mimimum at -6



