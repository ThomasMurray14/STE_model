% Script to simulate behavioural responses and RTs from a given model,
% using best fitting parameters


% Specify model
model = 'ResponseBias';
addpath(['../model_', model, '/']);
addpath('../helper_functions');

% Load model fits
fits = importdata(['../model_', model, '/model_', model, '_fit.mat']);

% Excluded bad participants
T = readtable('full_analysis\fits.csv');
excluded_IDs = T.ID(~logical(T.Include));
fits(ismember([fits.ID], excluded_IDs)) = [];

% Split by cond
safe_fits = fits(strcmp({fits.condition}, 'Safe'));
threat_fits = fits(strcmp({fits.condition}, 'Threat'));

% Get parameters
safe_prc_p = cell2mat(arrayfun(@(x) x.est.p_prc.p, safe_fits, 'UniformOutput', false)');
safe_obs_p = cell2mat(arrayfun(@(x) x.est.p_obs.p, safe_fits, 'UniformOutput', false)');
threat_prc_p = cell2mat(arrayfun(@(x) x.est.p_prc.p, threat_fits, 'UniformOutput', false)');
threat_obs_p = cell2mat(arrayfun(@(x) x.est.p_obs.p, threat_fits, 'UniformOutput', false)');

% Get mean parameters
safe_prc_p_mean = mean(safe_prc_p, 1);
safe_obs_p_mean = mean(safe_obs_p, 1);
threat_prc_p_mean = mean(threat_prc_p, 1);
threat_obs_p_mean = mean(threat_obs_p, 1);

% Get functions for this model
f = str2func(['STE_', model, '_config']);
[prc_config, obs_config] = f();

% Get inputs etc
sub_data = readtable('..\STE_data\10369536_A_Threat.csv');
[u,~] = data_prep(sub_data);

% Simulate
N=200;
[safe_resp, safe_logRT, threat_resp, threat_logRT] = deal(zeros(size(u,1), N));

for i = 1:N
    % Simulare safe
    sim_safe = tapas_simModel(u,...
        prc_config.model,...
        safe_prc_p_mean,...
        obs_config.model,...
        safe_obs_p_mean);
    resp_state_safe = sim_safe.y(:,1);
    safe_resp(:, i) = (u(:,2) == 1 & resp_state_safe == 1) + (u(:,2) == -1 & resp_state_safe == 0); % 1=Sad,0=Happy
    safe_logRT(:, i) = sim_safe.y(:,2);
    
    % Simulate threat
    sim_threat = tapas_simModel(u,...
        prc_config.model,...
        threat_prc_p_mean,...
        obs_config.model,...
        threat_obs_p_mean);
    resp_state_threat = sim_threat.y(:,1);
    threat_resp(:, i) = (u(:,2) == 1 & resp_state_threat == 1) + (u(:,2) == -1 & resp_state_threat == 0); % 1=Sad,0=Happy
    threat_logRT(:, i) = sim_threat.y(:,2);
end


% Get data indices
expected_idx = strcmp(sub_data.Expectedness, 'E');
unexpected_idx = strcmp(sub_data.Expectedness, 'UE');
stable_idx = sub_data.Block_N == 1;
volatile_idx = sub_data.Block_N > 1;


%% logRT, expectedness

% Get simulated logRTs
safe_logRT_expected_sim     = mean(safe_logRT(expected_idx, :), 1)';
safe_logRT_unexpected_sim   = mean(safe_logRT(unexpected_idx, :), 1)';
threat_logRT_expected_sim   = mean(threat_logRT(expected_idx, :), 1)';
threat_logRT_unexpected_sim = mean(threat_logRT(unexpected_idx, :), 1)';

% get real logRTs
safe_logRT_expected_real = T.Safe_expected_respLogRT;
safe_logRT_unexpected_real = T.Safe_unexpected_respLogRT;
threat_logRT_expected_real = T.Threat_expected_respLogRT;
threat_logRT_unexpected_real = T.Threat_unexpected_respLogRT;

% Plot
cols = colororder;
safe_col = cols(1,:);
threat_col = cols(2,:);
figure; hold on;
plot([0, 1], [mean(safe_logRT_expected_real, 'omitnan'), mean(safe_logRT_unexpected_real, 'omitnan')], 'color', safe_col, 'DisplayName', 'Real Safe');
plot([0, 1], [mean(threat_logRT_expected_real, 'omitnan'), mean(threat_logRT_unexpected_real, 'omitnan')], 'color', threat_col, 'DisplayName', 'Real Threat');
plot([0, 1], [mean(safe_logRT_expected_sim), mean(safe_logRT_unexpected_sim)], 'color', safe_col, 'linestyle', '--', 'DisplayName', 'Sim Safe');
plot([0, 1], [mean(threat_logRT_expected_sim), mean(threat_logRT_unexpected_sim)], 'color', threat_col, 'linestyle', '--', 'DisplayName', 'Sim Threat');
set(gca, ...
    'XLim', [0, 1], ...
    'XTick', [0, 1], ...
    'XTickLabels', {'Expected', 'Unexpected'});
ylabel('logRT');
legend;

%% logRT, volatility

% Get simulated logRTs
safe_logRT_stable_sim     = mean(safe_logRT(stable_idx, :), 1)';
safe_logRT_volatile_sim   = mean(safe_logRT(volatile_idx, :), 1)';
threat_logRT_stable_sim   = mean(threat_logRT(stable_idx, :), 1)';
threat_logRT_volatile_sim = mean(threat_logRT(volatile_idx, :), 1)';

% get real logRTs
safe_logRT_stable_real = T.Safe_stable_respLogRT;
safe_logRT_volatile_real = T.Safe_volatile_respLogRT;
threat_logRT_stable_real = T.Threat_stable_respLogRT;
threat_logRT_volatile_real = T.Threat_volatile_respLogRT;

% Plot
cols = colororder;
safe_col = cols(1,:);
threat_col = cols(2,:);
figure; hold on;
plot([0, 1], [mean(safe_logRT_stable_real, 'omitnan'), mean(safe_logRT_volatile_real, 'omitnan')], 'color', safe_col, 'DisplayName', 'Real Safe');
plot([0, 1], [mean(threat_logRT_stable_real, 'omitnan'), mean(threat_logRT_volatile_real, 'omitnan')], 'color', threat_col, 'DisplayName', 'Real Threat');
plot([0, 1], [mean(safe_logRT_stable_sim), mean(safe_logRT_volatile_sim)], 'color', safe_col, 'linestyle', '--', 'DisplayName', 'Sim Safe');
plot([0, 1], [mean(threat_logRT_stable_sim), mean(threat_logRT_volatile_sim)], 'color', threat_col, 'linestyle', '--', 'DisplayName', 'Sim Threat');
set(gca, ...
    'XLim', [0, 1], ...
    'XTick', [0, 1], ...
    'XTickLabels', {'Stable', 'Volatile'});
ylabel('logRT');

