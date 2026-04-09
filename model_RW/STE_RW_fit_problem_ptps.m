% Function to run parameter recovery for ResponseBias model

clear; close all; clc
addpath('../helper_functions/')


%% Get inputs
% example data (to get contingencies etc)
sub_data = readtable('..\STE_data\10369536_A_Threat.csv');
[u,y] = data_prep(sub_data);


%% Get configuration structures
[prc_config, obs_config] = STE_RW_config;
optim_config     = tapas_quasinewton_optim_config(); % optimisation algorithm
optim_config.nRandInit = 10;


%% Fit data

f_names = {'12641097_C_Safe.csv'};%, '12643580_C_Safe.csv', '12643580_C_Threat.csv', '12655060_D_Safe.csv'};

STE_dir = dir('..\STE_data\*.csv');
N_files = numel(f_names);
model_fits(N_files) = struct('ID', [], 'group', '', 'condition', '', 'est', struct());

% Just to be fancy
completion_times = zeros(N_files, 1);

for i = 1:N_files
    tic;
    fprintf('\nFitting dataset %i\n', i);
    if i>1
        avg_iter_time = mean(completion_times(1:i));
        fprintf('\tAverage iteration time = %1.2fs', avg_iter_time)
        estimated_total_time = avg_iter_time * ((N_files-i) + 1);
        fprintf('\n\tEstimated completion time = %im, %1.2fs\n\n', floor(estimated_total_time/60), rem(estimated_total_time,60));
    end

    % get ID, group, condition
    f_name = fullfile(STE_dir(i).folder, f_names{i});
    [~,name,~] = fileparts(f_name);
    tokens = split(name, '_');
    model_fits(i).ID           = str2double(tokens{1});
    model_fits(i).group        = tokens{2};
    model_fits(i).condition    = tokens{3};
    
    % load data
    sub_data = readtable(f_name);
    [~, y] = data_prep(sub_data);

    % remove missing
    missed = isnan(sub_data.Response_idx);
    y(missed, :) = NaN;

    % filter RT
    RT_filter = sub_data.Response_RT < 100;
    y(RT_filter, :) = NaN;

    % fit model
    try
        model_fits(i).est = tapas_fitModel(...
            y,...
            u,...
            prc_config,...
            obs_config,...
            optim_config);
    catch
    end
    completion_times(i) = toc;
end
