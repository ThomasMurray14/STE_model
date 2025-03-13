% SCript to fit model 6 to the data

%%
clear;
close all;

addpath('STE_model6');

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

% model input
u = [sub_data.u_al, cue];


%% Get model configs and set priors
% Ok so I'm using the obs model that Nace set up to use in the combined
% model, which uses 'correct' as the response. This isn't the same as I
% used for parameter recovery, but should (in theory) be the same. However
% I should probably run that again...

prc_model_config = prc3_ehgf_binary_pu_tbt_config(); % perceptual model
obs_model_config = obs1_unitsq_sgm_tbt_config(); % response model
obs_model_config.predorpost = 2; % use posterior

optim_config     = tapas_quasinewton_optim_config(); % optimisation algorithm
optim_config.nRandInit = 5;



prc_model_config.logalmu = log(0.4);
prc_model_config.logalsa = 2;
prc_model_config.rhomu(2) = 0;
prc_model_config.rhosa(2) = 2;
prc_model_config.ommu(2)=-2;
prc_model_config.omsa(2)=2;
prc_model_config.ommu(3)=0; % turn off third level
prc_model_config.omsa(3)=0; % turn off third level
prc_model_config.logkamu(2) = -inf; % turn off third level
prc_model_config.logkasa(2) = 0; % turn off third level
prc_model_config = tapas_align_priors(prc_model_config);

obs_model_config.logzemu = log(5);
obs_model_config.logzesa = 1;



%% Loop through data and fit model to each

STE_dir = dir('STE_data\*.csv');
N_files = numel(STE_dir);
model_fits(N_files) = struct('ID', [], 'group', '', 'condition', '', 'est', struct());

for i = 1:N_files
    % get ID, group, condition
    f_name = fullfile(STE_dir(i).folder, STE_dir(i).name);
    [~,name,~] = fileparts(f_name);
    tokens = split(name, '_');
    model_fits(i).ID           = str2double(tokens{1});
    model_fits(i).group        = tokens{2};
    model_fits(i).condition    = tokens{3};
    
    % load data
    sub_data = readtable(f_name);

    % model input
    u_sub = u;

    % subject responses
    % Use correct - turn's out Nace was on to something here
    sub_data.correct = sub_data.Outcome_idx == sub_data.Response_idx;
    y = sub_data.correct;

    % remove missing
    missed = isnan(sub_data.Response_idx);
    u_sub = u_sub(~missed, :);
    y = y(~missed, :);

    % fit model
    try
        model_fits(i).est = tapas_fitModel(...
            y,...
            u_sub,...
            prc_model_config,...
            obs_model_config,...
            optim_config);
    catch
    end
end

save('STE_model6_fits.mat', 'model_fits');


