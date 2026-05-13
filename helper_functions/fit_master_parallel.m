function model_fits = fit_master_parallel(u, prc_config, obs_config, optim_config, n_workers)
% master function to fit models to actual data
%
% input:
%   u             - model input
%   prc_config    - perceptual model config
%   obs_config    - response model config
%   optim_config  - optimisation algorithm config
%   n_workers     - number of parallel workers
% 
% 
% Would be nice to figure out a way of visualising completion.....

STE_dir = dir('..\STE_data\*.csv');
N_files = numel(STE_dir);

% Pre-allocate as a struct array (required for parfor indexing)
model_fits(N_files) = struct('ID', [], 'group', '', 'condition', '', 'est', struct());

% Start parallel pool with 8 workers (if not already running)
pool = gcp('nocreate');
if isempty(pool)
    parpool('local', n_workers);
elseif pool.NumWorkers ~= n_workers
    delete(pool);
    parpool('local', n_workers);
end

completion_times = zeros(N_files, 1);
start_time       = tic;

parfor i = 1:N_files
    iter_start = tic;

    % get ID, group, condition
    f_name = fullfile(STE_dir(i).folder, STE_dir(i).name);
    [~, name, ~] = fileparts(f_name);
    tokens = split(name, '_');

    % Build a local struct for this iteration (parfor requirement)
    result        = struct('ID', [], 'group', '', 'condition', '', 'est', struct());
    result.ID        = str2double(tokens{1});
    result.group     = tokens{2};
    result.condition = tokens{3};

    % load data
    sub_data = readtable(f_name);
    [~, y]   = data_prep(sub_data);

    % remove missing
    missed       = isnan(sub_data.Response_idx);
    y(missed, :) = NaN;

    % filter RT
    RT_filter       = sub_data.Response_RT < 100;
    y(RT_filter, :) = NaN;

    % fit model
    try
        result.est = fitModel( ...
            y, ...
            u, ...
            prc_config, ...
            obs_config, ...
            optim_config);
    catch err
        warning('fit_master:fitFailed', ...
            'Dataset %i failed: %s', i, err.message);
    end

    completion_times(i) = toc(iter_start);
    model_fits(i)       = result;
end

% Summary timing report (printed after all workers finish)
total_time = toc(start_time);
fprintf('\nAll %i datasets fitted.\n', N_files);
fprintf('Total wall time : %im %1.2fs\n', floor(total_time/60), rem(total_time,60));
fprintf('Mean per-dataset: %1.2fs\n',     mean(completion_times));
fprintf('Slowest dataset : %1.2fs\n',     max(completion_times));

end