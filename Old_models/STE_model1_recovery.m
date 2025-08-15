% Script to run parameter recovery for Model 1

%%
clear;
close all;

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


%% full parameter recovery
close all;

% reload model configs
prc_model_config = prc1_ehgf_binary_pu_tbt_config(); % perceptual model
obs_model_config = obs1_comb_obs_config(); % response model
optim_config     = tapas_quasinewton_optim_config(); % optimisation algorithm
optim_config.nRandInit = 10;

% set priors
prc_model_config.ommu(2)    =-3;
prc_model_config.omsa(2)    =4;
prc_model_config.rhosa(2)   =2;
prc_model_config.logalsa    =4;
prc_model_config = tapas_align_priors(prc_model_config);

obs_model_config.logzesa=2;
obs_model_config.beta0sa=2;
obs_model_config.beta1sa=2;
obs_model_config.beta2sa=2;
obs_model_config.beta3sa=2;
obs_model_config.beta4sa=0;
obs_model_config.logsasa=2;
obs_model_config = tapas_align_priors(obs_model_config);

% set N iterations
N = 300;

% preallocate sim and est parameters
[om2_sim, om2_est, rho2_sim, rho2_est, al_sim, al_est] = deal(nan(N,1));

for i = 1:N
    sim = tapas_sampleModel(u, prc_model_config, obs_model_config);
    om2_sim(i)  = sim.p_prc.om(2);
    rho2_sim(i) = sim.p_prc.rho(2);
    al_sim(i)   = sim.p_prc.al;
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
            end
        catch
            fprintf('\nCannot fit')
        end
    end
end

% figure('name', 'om2');  hold on; scatter(om2_sim, om2_est);     refline(1,0); xlabel('Simulated'); ylabel('Estimated'); title('Omega2');
% figure('name', 'rho2'); hold on; scatter(rho2_sim, rho2_est);   refline(1,0); xlabel('Simulated'); ylabel('Estimated'); title('Rho');
% figure('name', 'al');   hold on; scatter(al_sim, al_est);       refline(1,0); xlabel('Simulated'); ylabel('Estimated'); title('Alpha');

save('STE_model1_recovery.mat', 'om2_sim', 'om2_est', 'rho2_sim', 'rho2_est', 'al_sim', 'al_est');



