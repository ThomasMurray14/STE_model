function recovery_figures(recov)

all_params = fieldnames(recov);
all_params(ismember(all_params, {'LME', 'AIC', 'BIC'})) = [];

for iP = 1:numel(all_params)
    sim = recov.(all_params{iP}).sim;
    est = recov.(all_params{iP}).est;
    valid = ~isnan(est);
    sim = sim(valid);
    est = est(valid);
    figname = all_params{iP};
    if strcmp(recov.(all_params{iP}).space, 'log')
        sim = log(sim);
        est = log(est);
        figname = ['log(', figname, ')'];
    elseif strcmp(recov.(all_params{iP}).space, 'logit')
        % figure this out...
    end
    figure('name', figname);  
    hold on;
    scatter(sim, est);
    refline(1,0); 
    xlabel('Simulated');
    ylabel('Estimated');
    corr_string = sprintf("%s\nPearson's r = %.3f; Spearman's \\rho = %.3f; Kendall's \\tau = %.3f", ...
        figname,...
        corr(sim, est), ...
        corr(sim, est, 'type', 'Spearman'), ...
        corr(sim, est, 'type', 'Kendall'));
    title(corr_string)
end






end