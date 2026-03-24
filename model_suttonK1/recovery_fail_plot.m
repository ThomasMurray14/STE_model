% script to plot which simulated parameters fail to recover


recov = importdata('model_suttonK1_recovery.mat');


failed = isnan(recov.LME);

param_names = fieldnames(recov);
param_names = param_names(1:end-5);

for i = 1:numel(param_names)
    figure('name', param_names{i}); hold on;
    S = recov.(param_names{i}).sim(~failed);
    F = recov.(param_names{i}).sim(failed);
    
    if strcmp(recov.(param_names{i}).space, 'log')
        S = log(S);
        F = log(F);
    end

    scatter(S, ones(size(S)) + .1, 'filled','blue');
    scatter(F, ones(size(F)), 'filled','red');
end



figure;
i = 1;
for iR = 1:numel(param_names)
    for iC = 1:numel(param_names)
        subplot(numel(param_names), numel(param_names), i);
        hold on;
        
        s_r = recov.(param_names{iR}).sim(~failed);
        f_r = recov.(param_names{iR}).sim(failed);
        s_c = recov.(param_names{iC}).sim(~failed);
        f_c = recov.(param_names{iC}).sim(failed);
        
        if strcmp(recov.(param_names{iR}).space, 'log')
            s_r = log(s_r);
            f_r = log(f_r);
        end

        if strcmp(recov.(param_names{iC}).space, 'log')
            s_c = log(s_c);
            f_c = log(f_c);
        end
        
        scatter(s_r, s_c, 'filled', 'blue');
        scatter(f_r, f_c, 'filled', 'red');

        i = i+1;
    end
end


