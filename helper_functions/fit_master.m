function model_fits = fit_master(u, prc_config, obs_config, optim_config)
% master function to fit models to actual data
% 
% input:
%   u             - model input
%   prc_config    - perceptual model config
%   obs_config    - response model config
%   optim_config  - optimisation algorithm config
% 


STE_dir = dir('..\STE_data\*.csv');
N_files = numel(STE_dir);
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
    f_name = fullfile(STE_dir(i).folder, STE_dir(i).name);
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


end