% Script to model the STE data, using Nace custom HGF but just using binary
% responses (no RT)




close all; clear;


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

% response_idx    = sub_data.Response_idx;
sub_data.logRT              = log(sub_data.Response_RT);
sub_data = sub_data(~isnan(sub_data.logRT),:); % remove nans 

sub_data.correct = sub_data.Outcome_idx == sub_data.Response_idx;

%happy_face - 0.5

u = [sub_data.u_al, sub_data.Cue_idx];
y = [sub_data.correct];


%% models
prc_model_config = prc1_ehgf_binary_pu_tbt_config(); % perceptual model
obs_model_config = tapas_unitsq_sgm_config();%tapas_logrt_linear_binary_config(); % response model
optim_config     = tapas_quasinewton_optim_config(); % optimisation algorithm
optim_config.nRandInit = 5;


%% simulate responses
r_temp = [];
r_temp.c_prc.n_levels = 3;
prc_params = prc1_ehgf_binary_pu_tbt_transp(r_temp, prc_model_config.priormus);

prc_params(1); %mu_0mu(1)
prc_params(2); %mu_0mu(2)
prc_params(3); %mu_0mu(3)
prc_params(4); %logsa_0mu(1);
prc_params(5); %logsa_0mu(2);
prc_params(6); %logsa_0mu(3);
prc_params(7); %rho(1);
prc_params(8) = 100; %rho(2);
prc_params(9); %rho(3);
prc_params(10)=.001; %logkamu(1);
prc_params(11); %logkamu(2);
prc_params(12); %ommu(1);
prc_params(13)=-1; %ommu(2);
prc_params(14)=-1; %ommu(3);
prc_params(15) = 0.002; %logalmu;
prc_params(16); %eta0mu;
prc_params(17); %eta1mu;



obs_params = obs_model_config.priormus;

sim = tapas_simModel(u,...%[u;u;u;u],...
    'prc1_ehgf_binary_pu_tbt',...
    prc_params,...
    'tapas_unitsq_sgm',...
    obs_params, ...
    123456789);

% visualise psychometric
sim_psychometric = arrayfun(@(x) mean(sim.y(sub_data.Outcome_p_sad==x, 1)), 0:20:100);
% figure('name', 'simulated psychometric'); hold on;
plot(0:20:100, sim_psychometric, 'linewidth', 3);
set(gca, 'Ylim', [0,1], 'Xtick', 0:20:100)



