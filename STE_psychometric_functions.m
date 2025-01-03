
% Script to play around with modelling the TOSSTE data
% close all; clear; clc;

addpath([pwd, '\custom hgf']);

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

u = [sub_data.u_al, sub_data.Cue_idx];
y = [sub_data.correct,sub_data.logRT];


%% Get configuration structures
prc_model_config = tapas_ehgf_binary_pu_tbt_config(); % perceptual model
obs_model_config = m1_comb_obs_config();%tapas_logrt_linear_binary_config(); % response model
optim_config     = tapas_quasinewton_optim_config(); % optimisation algorithm

optim_config.nRandInit = 5;



r_temp = [];
r_temp.c_prc.n_levels = 3;
prc_params = tapas_ehgf_binary_pu_tbt_transp(r_temp, prc_model_config.priormus);

trials = sub_data(:, ["Trial_N", "Contingency", "Cue_idx", "Outcome_idx", "Outcome_p_sad"]);

rhos = 10;

sim_ys = zeros(192, numel(rhos));
ps = [0,20,40,60,80,100];

p_sads=zeros(numel(rhos), 6);

for i=1:numel(rhos)

    prc_params(1); %mu_0mu(1)
    prc_params(2); %mu_0mu(2)
    prc_params(3); %mu_0mu(3)
    prc_params(4); %logsa_0mu(1);
    prc_params(5); %logsa_0mu(2);
    prc_params(6); %logsa_0mu(3);
    prc_params(7); %rho(1);
    prc_params(8) = rho; %rho(2);
    prc_params(9); %rho(3);
    prc_params(10); %logkamu(1);
    prc_params(11); %logkamu(2);
    prc_params(12); %ommu(1);
    prc_params(13); %ommu(2);
    prc_params(14); %ommu(3);
    prc_params(15) =.02; %logalmu;
    prc_params(16); %eta0mu;
    prc_params(17); %eta1mu;
    
    obs_params = obs_model_config.priormus;
    obs_params(1) = exp(obs_params(1));
    obs_params(7) = exp(obs_params(7));
    
    sim = tapas_simModel(u,...
        'tapas_ehgf_binary_pu_tbt',...
        prc_params,...
        'm1_comb_obs',...
        obs_params,...
        123456789);
    
    sim_y = sim.y(:,1);
    sim_ys(:, i) = sim_y;

    
    for ip=1:6
        p_sads(i, ip) = sum(sim_y(trials.Outcome_p_sad==ps(ip)))/32;
    end


end


%%
close all;
plot(1:6, p_sads)
set(gca, 'Ylim', [0, 1])


