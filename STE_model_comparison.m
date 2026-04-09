% Model comparison using VBA toolbox


% Name models as folder names
model_names = {'PerceptualBias', 'PerceptualBias3L', 'PredictionBias', 'PredictionBias3L',...
    'ResponseBias', 'ResponseBias3L', 'RW', 'SuttonK1'};
N_models = numel(model_names);

%% Load model fits
for iM = 1:N_models
    models(iM).model_fits = importdata(['model_', model_names{iM}, '\model_', model_names{iM}, '_fit.mat']);
end

%% Load excluded participants
fits = readtable('full_analysis\fits.csv');
IDs = fits.ID;
include = logical(fits.Include);
excluded_IDs = IDs(~include);

% remove IDs from model_fits
excluded_idx = ismember([models(1).model_fits.ID], excluded_IDs);
for i = 1:N_models
    models(i).model_fits(excluded_idx) = [];
end

N_inc = numel(models(1).model_fits);

%% Get LME
LMEs = nan(N_inc, N_models);
for iP = 1:N_inc
    for iM = 1:N_models
        if ~isempty(models(iM).model_fits(iP).est)
            LMEs(iP, iM) = models(iM).model_fits(iP).est.optim.LME;
        end
    end
end
valid = ~(isnan(LMEs) + isinf(LMEs));
LMEs_valid = LMEs(~any(~valid, 2), :);


%%% Check N ~valid for each model, exclude if necessary


%% Model comparison
options.modelNames = model_names;
[posterior,out] = VBA_groupBMC(LMEs', options);


