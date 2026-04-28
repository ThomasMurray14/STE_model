% script to write parameter estimates to csv

model = 'ResponseBias';
fits = importdata(['..\model_', model, '\model_', model, '_fit.mat']);

% Specify parameters
p_names = {'om2', 'al', 'zeta0', 'zeta1', 'beta0', 'beta1', 'beta2', 'beta3', 'beta4', 'sa'};
p_idx = [13, 15, 1, 2, 3, 4, 5, 6, 7, 8];
p_mod = {'prc', 'prc', 'obs', 'obs', 'obs', 'obs', 'obs', 'obs', 'obs', 'obs'};
p_space = {'native', 'log', 'log', 'native', 'native', 'native', 'native', 'native', 'native', 'log'};

IDs = unique([fits(:).ID]);

S = struct();
conds = {'Safe', 'Threat'};

% Loop subjects
for iSub = 1:numel(IDs)
    % Add ID
    S(iSub).ID = IDs(iSub);
    
    % Loop through safe/threat
    for iCond = 1:2

        % Get fits for this ID and condition
        i_fit = fits([fits.ID] == IDs(iSub) & strcmp({fits.condition}, conds{iCond}));

        % if this ID/condition has been estimated
        if ~isempty(i_fit.est)

            % Loop parameters
            for iP = 1:numel(p_names)
                S(iSub).([conds{iCond}, '_', p_names{iP}]) = i_fit.est.(['p_', p_mod{iP}]).p(p_idx(iP));

                if strcmp(p_space{iP}, 'log')
                    S(iSub).([conds{iCond}, '_log', p_names{iP}]) = log(i_fit.est.(['p_', p_mod{iP}]).p(p_idx(iP)));
                end
            end

        else % fill with NaN
            for iP = 1:numel(p_names)
                S(iSub).([conds{iCond}, '_', p_names{iP}]) = NaN;
                if strcmp(p_space{iP}, 'log')
                    S(iSub).([conds{iCond}, '_log', p_names{iP}]) = NaN;
                end
            end

        end

    end
end

% Convert to table
T = struct2table(S);

% Remove bad fits - some settle on priors
for iCond = 1:2
    bad_idx = find(T.([conds{iCond} '_om2']) == -2);
    for i = 1:numel(bad_idx)
        if iCond == 1
            T(bad_idx(i), 2:11) = {NaN};
        elseif iCond == 2
            T(bad_idx(i), 12:21) = {NaN};
        end
    end
end


% Get change scores (Threat-Safe)
for iP = 1:numel(p_names)
    T.(['Delta_', p_names{iP}]) = T.(['Threat_', p_names{iP}]) - T.(['Safe_', p_names{iP}]);
    if strcmp(p_space{iP}, 'log')
        T.(['Delta_log', p_names{iP}]) = T.(['Threat_log', p_names{iP}]) - T.(['Safe_log', p_names{iP}]);
    end
end

% Write table
writetable(T, ['../full_analysis/', model, '_params.csv']);


