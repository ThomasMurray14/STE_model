function param_table = print_fits_results(m, prc_param_names, prc_param_idx, obs_param_names, obs_param_idx)
% function to print summary of fit parameters

% e.g.
% prc_param_names = {'om2', 'rho2', 'al'};
% prc_param_idx   = [13, 8, 15];
% obs_param_names = {'ze', 'beta0', 'beta1', 'beta2', 'beta3', 'beta4', 'sa'};
% obs_param_idx   = [1, 2, 3, 4, 5, 6, 7];

% get list of IDs
IDs = unique([m.ID])';
N = numel(IDs);

% Set ID in table
param_table = table;
param_table.ID = IDs;

% preallocate space in table
all_params = [prc_param_names, obs_param_names];
for iParam = 1:numel(all_params)
    param_table.([all_params{iParam}, '_S']) = nan(N, 1);
    param_table.([all_params{iParam}, '_T']) = nan(N, 1);
end

% add parameters to table
for i = 1:numel(m)
    idx = param_table.ID == m(i).ID;
    c = m(i).condition(1); %S or T
    
    % if estimated
    if ~isempty(m(i).est)
        % loop through prc params
        for iPrc = 1:numel(prc_param_names)
            param_table(idx, [prc_param_names{iPrc}, '_', c]) = {m(i).est.p_prc.p(prc_param_idx(iPrc))};
        end
        
        % loop through obs params
        for iObs = 1:numel(obs_param_names)
            param_table(idx, [obs_param_names{iObs}, '_', c]) = {m(i).est.p_obs.p(obs_param_idx(iObs))};
        end
    end
end

% print summary
for iP = 1:numel(all_params)
    fprintf('\nParam = %s', all_params{iP});
    S = param_table{:, [all_params{iP}, '_S']};
    T = param_table{:, [all_params{iP}, '_T']};

    fprintf('\n\tMean (S.D.):')
    fprintf('\n\t\tSafe=%1.3f (%1.3f)', mean(S, 'omitmissing'), std(S, 'omitmissing'));
    fprintf('\n\t\tThreat=%1.3f (%1.3f)', mean(T, 'omitmissing'), std(T, 'omitmissing'));

    
    [h,p,ci,stats] = ttest(S, T);
    fprintf('\n\tT-test:')
    fprintf('\n\t\tt(%i)=%1.3f, p=%1.3f', stats.df, stats.tstat, p);
    if h == 1; fprintf('**'); end

    fprintf('\n')
end



end