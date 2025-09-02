function recov = parameter_recovery_master(u, prc_model_config, obs_model_config, optim_config, N, prc_params, obs_params, figs)
% master function to perform parameter recovery
% 
% input:
%   u                   - model input
%   prc_model_config
%   obs_model_config
%   optim_config
%   N                   - N simulations
%   prc_params          - cell with names of perceptual model params
%   obs_params          - cell with names of observation model params
%   figs                - boolean to produce figs


% preallocate sim and est parameters
all_params = [prc_params, obs_params];
for iP = 1:numel(all_params)
    recov.(all_params{iP}).sim = zeros(N, 1);
    recov.(all_params{iP}).est = zeros(N, 1);
end


% Main loop
for i = 1:N
    fprintf('\nParameter recovery iteration %i\n', i);
    sim = tapas_sampleModel(u, prc_model_config, obs_model_config);
    
    % Store simulated prc params
    for iP = 1:numel(prc_params)
        param_name=prc_params{iP};
        param_idx=str2double(param_name(end)); % clunky method for dealing with some parameters ending in idx (e.g. om2)
        if isnan(param_idx)
            recov.(param_name).sim(i) = sim.p_prc.(param_name);
        else
            recov.(param_name).sim(i) = sim.p_prc.(param_name(1:end-1))(param_idx);
        end
    end
    
    % store simulated obs params
    for iP = 1:numel(obs_params)
        param_name=obs_params{iP};
        param_idx=str2double(param_name(end)); % clunky method for dealing with some parameters ending in idx (e.g. om2)
        if isnan(param_idx)
            recov.(param_name).sim(i) = sim.p_obs.(param_name);
        else
            recov.(param_name).sim(i) = sim.p_obs.(param_name(1:end-1))(param_idx);
        end
    end


    % estimate parameters from simulated responses
    if ~any(isnan(sim.y)) % only if no missing trials
        try
            est = tapas_fitModel(...
                sim.y,...
                sim.u,...
                prc_model_config,...
                obs_model_config,...
                optim_config);

            % get estimated parameters
            if ~isinf(est.optim.LME)

                % Store estimated prc params
                for iP = 1:numel(prc_params)
                    param_name=prc_params{iP};
                    param_idx=str2double(param_name(end)); % clunky method for dealing with some parameters ending in idx (e.g. om2)
                    if isnan(param_idx)
                        recov.(param_name).est(i) = est.p_prc.(param_name);
                    else
                        recov.(param_name).est(i) = est.p_prc.(param_name(1:end-1))(param_idx);
                    end
                end
                
                % store simulated obs params
                for iP = 1:numel(obs_params)
                    param_name=obs_params{iP};
                    param_idx=str2double(param_name(end)); % clunky method for dealing with some parameters ending in idx (e.g. om2)
                    if isnan(param_idx)
                        recov.(param_name).est(i) = est.p_obs.(param_name);
                    else
                        recov.(param_name).est(i) = est.p_obs.(param_name(1:end-1))(param_idx);
                    end
                end

            end
        catch
            fprintf('\nCannot fit')
        end
    end
end


if figs
    for iP = 1:numel(all_params)
        figure('name', all_params{iP});  
        hold on; 
        scatter(recov.(all_params{iP}).sim, recov.(all_params{iP}).est);
        refline(1,0); 
        xlabel('Simulated');
        ylabel('Estimated');
        title(all_params{iP})
    end
end


end