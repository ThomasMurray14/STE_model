%% MODEL 2
% Perceptual model = 2 level eHGF
% Response model = combined psychometric and logRT
% 
% 
% TODO
%   parameters/config for obs2
% 

#####
%  SCRIPTS ARE NOT FINISHED
#####


%%
clear;
close all;
addpath('STE_model_psychometric');

%% 

% example data (to get contingencies etc)
sub_data = readtable('STE_data\10369536_A_Threat.csv');

% Contingency space
cue = sub_data.Cue_idx;
cue(cue==0) = -1; % balanced contrast coding for the "happy bias" (not needed in this model)
outcome = sub_data.Outcome_p_sad/100;
u_al = 1 - (0.5*(1+cue)-cue.*outcome);
state = double(sub_data.Cue_idx == sub_data.Outcome_idx);

sub_data.u_al = u_al;
sub_data.state=state;
sub_data.p_sad=outcome;

% Responses
sub_data.logRT = log(sub_data.Response_RT);
sub_data.resp_state = double(sub_data.Cue_idx == sub_data.Response_idx);

% model input
u = sub_data.u_al;
y = [sub_data.resp_state, sub_data.logRT];


%% Get configuration structures
prc_model_config = prc2_ehgf_binary_config(); % perceptual model
obs_model_config = obs2_comb_obs_config(); % response model
optim_config     = tapas_quasinewton_optim_config(); % optimisation algorithm
optim_config.nRandInit = 5;












