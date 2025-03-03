% Script to get model comparison metrics

%% load model fits
models(1).model_fits = importdata('STE_model1_fits.mat');
models(2).model_fits = importdata('STE_model2_fits.mat');
models(3).model_fits = importdata('STE_model3_fits.mat');
models(4).model_fits = importdata('STE_model4_fits.mat');
models(5).model_fits = importdata('STE_model5_fits.mat');

model_names = {
    'eHGF (Nace remix)';
    'eHGF-psychometric (PSE)';
    'eHGF-psychometric (slope)';
    'VKF-psychometric (PSE)';
    'VKF-psychometric (slope)'};

%% Safe/Threat Concatenated
LMEs = nan(400, 5);
BICs = nan(400, 5);
AICs = nan(400, 5);
for iP = 1:400
    for iM = 1:5
        if ~isempty(models(iM).model_fits(iP).est)
            LMEs(iP, iM) = models(iM).model_fits(iP).est.optim.LME;
            BICs(iP, iM) = models(iM).model_fits(iP).est.optim.BIC;
            AICs(iP, iM) = models(iM).model_fits(iP).est.optim.AIC;
        end
    end
end


valid = ~(isnan(LMEs) + isinf(LMEs));
LMEs_valid = LMEs(~any(~valid, 2), :);
BICs_valid = BICs(~any(~valid, 2), :);
AICs_valid = AICs(~any(~valid, 2), :);
BICs_sum = sum(BICs_valid, 1);
AICs_sum = sum(AICs_valid, 1);

[alpha,exp_r,xp,pxp,bor] = spm_BMS(LMEs_valid);

figure;
bar(1:5, alpha, 'EdgeColor', 'none');
set(gca, 'xticklabels', model_names)
title('BMS')

N_fit = 400 - (sum(isnan(LMEs), 1) + sum(isinf(LMEs), 1));
fprintf('\nN fit:')
for i = 1:5
    fprintf('\t\n%s: %i', model_names{i}, N_fit(i));
end
fprintf('\nBayes Ombinus Risk: %1.3f', bor);
fprintf('\n\nAIC (sum):')
for i = 1:5
    fprintf('\t\n%s: %1.1f', model_names{i}, AICs_sum(i));
end
fprintf('\n\nBIC (sum):')
for i = 1:5
    fprintf('\t\n%s: %1.1f', model_names{i}, BICs_sum(i));
end




%% Safe only
safe_idx = find(arrayfun(@(x) strcmp(models(1).model_fits(x).condition, 'Safe'), 1:400));

LMEs = nan(200, 5);
BICs = nan(200, 5);
AICs = nan(200, 5);
for i = 1:200
    iP = safe_idx(i);
    for iM = 1:5
        if ~isempty(models(iM).model_fits(iP).est)
            LMEs(i, iM) = models(iM).model_fits(iP).est.optim.LME;
            BICs(i, iM) = models(iM).model_fits(iP).est.optim.BIC;
            AICs(i, iM) = models(iM).model_fits(iP).est.optim.AIC;
        end
    end
end

valid = ~(isnan(LMEs) + isinf(LMEs));
LMEs_valid = LMEs(~any(~valid, 2), :);
BICs_valid = BICs(~any(~valid, 2), :);
AICs_valid = AICs(~any(~valid, 2), :);
BICs_sum = sum(BICs_valid, 1);
AICs_sum = sum(AICs_valid, 1);
[alpha,exp_r,xp,pxp,bor] = spm_BMS(LMEs_valid);

figure;
bar(1:5, alpha, 'EdgeColor', 'none');
set(gca, 'xticklabels', model_names)
title('BMS (Safe only)')

N_fit = 200 - (sum(isnan(LMEs), 1) + sum(isinf(LMEs), 1));
fprintf('\nN fit:')
for i = 1:5
    fprintf('\t\n%s: %i', model_names{i}, N_fit(i));
end
fprintf('\nBayes Ombinus Risk: %1.3f', bor);
fprintf('\n\nAIC (sum):')
for i = 1:5
    fprintf('\t\n%s: %1.1f', model_names{i}, AICs_sum(i));
end
fprintf('\n\nBIC (sum):')
for i = 1:5
    fprintf('\t\n%s: %1.1f', model_names{i}, BICs_sum(i));
end



%% Threat only
threat_idx = find(arrayfun(@(x) strcmp(models(1).model_fits(x).condition, 'Threat'), 1:400));

LMEs = nan(200, 5);
BICs = nan(200, 5);
AICs = nan(200, 5);
for i = 1:200
    iP = threat_idx(i);
    for iM = 1:5
        if ~isempty(models(iM).model_fits(iP).est)
            LMEs(i, iM) = models(iM).model_fits(iP).est.optim.LME;
            BICs(i, iM) = models(iM).model_fits(iP).est.optim.BIC;
            AICs(i, iM) = models(iM).model_fits(iP).est.optim.AIC;
        end
    end
end

valid = ~(isnan(LMEs) + isinf(LMEs));
LMEs_valid = LMEs(~any(~valid, 2), :);
BICs_valid = BICs(~any(~valid, 2), :);
AICs_valid = AICs(~any(~valid, 2), :);
BICs_sum = sum(BICs_valid, 1);
AICs_sum = sum(AICs_valid, 1);
[alpha,exp_r,xp,pxp,bor] = spm_BMS(LMEs_valid);

figure;
bar(1:5, alpha, 'EdgeColor', 'none');
set(gca, 'xticklabels', model_names)
title('BMS (Threat only)')

N_fit = 200 - (sum(isnan(LMEs), 1) + sum(isinf(LMEs), 1));
fprintf('\nN fit:')
for i = 1:5
    fprintf('\t\n%s: %i', model_names{i}, N_fit(i));
end
fprintf('\n\nBayes Ombinus Risk: %1.3f', bor);
fprintf('\n\nAIC (sum):')
for i = 1:5
    fprintf('\t\n%s: %1.1f', model_names{i}, AICs_sum(i));
end
fprintf('\n\nBIC (sum):')
for i = 1:5
    fprintf('\t\n%s: %1.1f', model_names{i}, BICs_sum(i));
end








