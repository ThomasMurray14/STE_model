% Function to run parameter recovery for PerceptualBias model


% TO DO:
%   Add est and sim structures to recov (see SuttonK1 parameter recovery)



clear; close all; clc
addpath('..');

%% Get inputs
% example data (to get contingencies etc)
sub_data = readtable('..\STE_data\10369536_A_Threat.csv');

% Contingency space
cue = sub_data.Cue_idx;
cue(cue==0) = -1; % balanced contrast coding for the "happy bias"
outcome = sub_data.Outcome_p_sad/100;
u_al = 1 - (0.5*(1+cue)-cue.*outcome);
state = double(sub_data.Cue_idx == sub_data.Outcome_idx);

sub_data.u_al = u_al;
sub_data.state=state;
sub_data.p_sad=outcome;

% Responses
sub_data.logRT = log(sub_data.Response_RT);
sub_data.resp_state = double(sub_data.Cue_idx == sub_data.Response_idx);

% model input
u = [sub_data.u_al, cue];
y = [sub_data.resp_state, sub_data.logRT];


%% Get configuration structures
[prc_config, obs_config] = STE_PerceptualBias3L_config;
optim_config     = tapas_quasinewton_optim_config(); % optimisation algorithm
optim_config.nRandInit = 5;


%% Run parameter recovery

% run recovery
N=200;

% Parameters to recover
prc_param_names = {'om2', 'om3', 'rho2', 'al'};
prc_param_idx   = [13, 14, 8, 15];
prc_param_space = {'native', 'native', 'native', 'log'};
obs_param_names = {'ze', 'beta0', 'beta1', 'beta2', 'beta3', 'beta4', 'beta5', 'sa'};
obs_param_idx   = [1, 2, 3, 4, 5, 6, 7, 8];
obs_param_space = {'log', 'native', 'native', 'native', 'native', 'native', 'native', 'log'};

% preallocate sim and est parameters
all_params = [prc_param_names, obs_param_names];
for iP = 1:numel(all_params)
    recov.(all_params{iP}).sim = nan(N, 1);
    recov.(all_params{iP}).est = nan(N, 1);
end
recov.LME = nan(N, 1); % store LME
recov.AIC = nan(N, 1); % store AIC
recov.BIC = nan(N, 1); % store BIC
recov.est = cell(N,1); % store full est structure

% Just to be fancy
completion_times = zeros(N, 1);

% Main loop
for i = 1:N
    tic;
    fprintf('\n--------------------------------\n')
    fprintf('\nParameter recovery iteration %i\n', i);
    if i>1
        avg_iter_time = mean(completion_times(1:i));
        fprintf('\tAverage iteration time = %1.2fs', avg_iter_time)
        estimated_total_time = avg_iter_time * ((N-i) + 1);
        fprintf('\n\tEstimated completion time = %im, %1.2fs\n\n', floor(estimated_total_time/60), rem(estimated_total_time,60));
    end

    sim = tapas_sampleModel(u, prc_config, obs_config);
    recov.sim{i} = sim;
    
    % Store simulated prc params
    for iP = 1:numel(prc_param_names)
        param_name=prc_param_names{iP};
        recov.(param_name).sim(i) = sim.p_prc.p(prc_param_idx(iP));
        recov.(param_name).space = prc_param_space{iP};
    end
    
    % store simulated obs params
    for iP = 1:numel(obs_param_names)
        param_name=obs_param_names{iP};
        recov.(param_name).sim(i) = sim.p_obs.p(obs_param_idx(iP));
        recov.(param_name).space = obs_param_space{iP};
    end


    % estimate parameters from simulated responses
    if ~any(isnan(sim.y)) % only if no missing trials
        try
            est = tapas_fitModel(...
                sim.y,...
                sim.u,...
                prc_config,...
                obs_config,...
                optim_config);
            recov.est{i} = est; % store est (for debugging)

            % get estimated parameters
            if ~isinf(est.optim.LME)
                % store fit metrics
                recov.LME(i) = est.optim.LME;
                recov.AIC(i) = est.optim.AIC;
                recov.BIC(i) = est.optim.BIC;
                
                % Store estimated prc params
                for iP = 1:numel(prc_param_names)
                    param_name=prc_param_names{iP};
                    recov.(param_name).est(i) = est.p_prc.p(prc_param_idx(iP));
                end
                
                % store simulated obs params
                for iP = 1:numel(obs_param_names)
                    param_name=obs_param_names{iP};
                    recov.(param_name).est(i) = est.p_obs.p(obs_param_idx(iP));
                end

            end
        catch
            fprintf('\nCannot fit')
        end
        completion_times(i) = toc;
    end
end

save('PerceptualBias3L_recovery.mat', 'recov');
recovery_figures(recov);



