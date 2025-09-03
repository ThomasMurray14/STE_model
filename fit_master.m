function model_fits = fit_master(u, prc_model_config, obs_model_config, optim_config)
% master function to fit models to actual data
% 
% input:
%   u                   - model input. Missed responses for each subject
%                           are removed
%   prc_model_config
%   obs_model_config
%   optim_config

% 
% NOTE - this only works for models using resp_state as y (model 1 so far).
% If I make other models that fit to response space (not state) then this
% won't work.

STE_dir = dir('STE_data\*.csv');
N_files = numel(STE_dir);
model_fits(N_files) = struct('ID', [], 'group', '', 'condition', '', 'est', struct());

for i = 1:N_files
    fprintf('\nFitting dataset %i\n', i)

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
    sub_data.logRT = log(sub_data.Response_RT);
    sub_data.resp_state = double(sub_data.Cue_idx == sub_data.Response_idx);
    y = [sub_data.resp_state, sub_data.logRT];%, sub_data.Confidence_idx];

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



end