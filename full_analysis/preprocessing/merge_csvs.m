% Script to merge csvs from different outputs

psychometric_functions = readtable('psychometric_function_fits.csv');
attention_check = readtable('attention_check.csv');
demographics = readtable('demographics.csv');
transdiagnostic_factors = readtable('predictedFactorScores.csv');
VAS = readtable('state_anxiety_VAS.csv');
model_params = readtable('..\ResponseBias_params2.csv');

% Make sure all have 'ID'
psychometric_functions = renamevars(psychometric_functions, 'Var1', 'ID');
demographics = renamevars(demographics, 'subjIDs', 'ID');
transdiagnostic_factors = renamevars(transdiagnostic_factors, 'subjIDs', 'ID');

% Remove group from psychometric (demo contains group, will come first)
psychometric_functions = removevars(psychometric_functions, 'Group');

% Binarise transdiagnostic (for visualisations)
transdiagnostic_factors = removevars(transdiagnostic_factors, 'Group');
transdiagnostic_factors.AD_binary = transdiagnostic_factors.AD > median(transdiagnostic_factors.AD);
transdiagnostic_factors.Compul_binary = transdiagnostic_factors.Compul > median(transdiagnostic_factors.Compul);
transdiagnostic_factors.SW_binary = transdiagnostic_factors.SW > median(transdiagnostic_factors.SW);

% VAS effects
VAS = removevars(VAS, 'Group');
VAS = renamevars(VAS, 'expt_start', 'VAS_expt_start');
VAS = renamevars(VAS, 'post_ratings', 'VAS_post_ratings');
VAS = renamevars(VAS, 'pre_safe', 'VAS_pre_safe');
VAS = renamevars(VAS, 'pre_threat', 'VAS_pre_threat');
VAS = renamevars(VAS, 'post_threat', 'VAS_post_threat');
VAS.VAS_safe_vs_baseline = VAS.VAS_pre_safe - VAS.VAS_expt_start;
VAS.VAS_threat_vs_baseline = VAS.VAS_pre_threat - VAS.VAS_expt_start;
VAS.VAS_threat_vs_safe = VAS.VAS_pre_threat - VAS.VAS_pre_safe;

% Attention check
attention_check = removevars(attention_check, 'Group');
attention_check = renamevars(attention_check, 'Safe_Check_1', 'Attention_check_safe1');
attention_check = renamevars(attention_check, 'Safe_Check_2', 'Attention_check_safe2');
attention_check = renamevars(attention_check, 'Safe_Check_3', 'Attention_check_safe3');
attention_check = renamevars(attention_check, 'Threat_Check_1', 'Attention_check_threat1');
attention_check = renamevars(attention_check, 'Threat_Check_2', 'Attention_check_threat2');
attention_check = renamevars(attention_check, 'Threat_Check_3', 'Attention_check_threat3');
attention_check.Attention_check_sum = ...
    attention_check.Attention_check_safe1 + ...
    attention_check.Attention_check_safe2 + ...
    attention_check.Attention_check_safe3 + ...
    attention_check.Attention_check_threat1 + ...
    attention_check.Attention_check_threat2 + ...
    attention_check.Attention_check_threat3;



% Merge
tables = {demographics, psychometric_functions, transdiagnostic_factors, VAS, attention_check, model_params};
T = tables{1};
for i = 2:numel(tables)
    T = outerjoin(T, tables{i}, 'Keys', 'ID', 'MergeKeys', true);
end

% Add inclusion criteria
T.Attention_check_include = T.Attention_check_sum > 4; % i.e. missed no more than 2
T.N_missing_safe_include = T.Safe_N_missed < 20;
T.N_missing_threat_include = T.Threat_N_missed < 20;
T.Include = T.Attention_check_include & T.N_missing_safe_include & T.N_missing_threat_include;

% save
writetable(T, '../fits2.csv');




