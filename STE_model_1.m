%% TODO

% Try including confidence



%%
clear;
close all;

addpath('STE_model1_with_conf');

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
sub_data.resp_state = double(sub_data.Cue_idx == sub_data.Response_idx);


% model input
u = [sub_data.u_al, cue];
y = [sub_data.resp_state, sub_data.logRT, sub_data.Confidence_idx];

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
obs_params(8) = 3;%exp(obs_params(8));

sim = tapas_simModel(u_sub,...
    'prc1_ehgf_binary_pu_tbt',...
    prc_params,...
    'obs1_comb_obs',...
    obs_params,...
    123456789);

%%
% Visualise Psychometric
sim_sad = (sub_data.Cue_idx == 1 & sim.y(:,1) == 1) + (sub_data.Cue_idx == 0 & sim.y(:,1) == 0);

figure('name', 'simulated psychometric'); hold on;
sim_psychometric = arrayfun(@(x) mean(sim_sad(sub_data.Outcome_p_sad==x, 1)), 0:20:100);
plot(0:20:100, sim_psychometric, 'linewidth', 3);
set(gca, 'Ylim', [0,1], 'Xtick', 0:20:100)

% Visualise confidence
figure('name', 'simulated confidence'); hold on;
sim_conf = sim.y(:,3);
sim_gaussian = arrayfun(@(x) mean(sim_conf(sub_data.Outcome_p_sad==x, 1)), 0:20:100);
plot(0:20:100, sim_gaussian, 'linewidth', 3);
set(gca, 'Ylim', [0,1], 'Xtick', 0:20:100)

%% Plot trajectory
prc1_ehgf_binary_tbt_plotTraj(sim);


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

