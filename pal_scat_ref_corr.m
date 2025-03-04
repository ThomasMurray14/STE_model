function f = pal_scat_ref_corr(sim, est, name, varargin)

% varargin = fig name (with path)

valid = ~isnan(est);
est = est(valid);
sim = sim(valid);

scatter(sim,est, [], 'blue', 'filled', 'MarkerFaceAlpha', 0.5)

corr_string = sprintf("%s\nPearson's r = %.3f; Spearman's \\rho = %.3f; Kendall's \\tau = %.3f", ...
    name,...
    corr(sim, est), ...
    corr(sim, est, 'type', 'Spearman'), ...
    corr(sim, est, 'type', 'Kendall'));

set_title = corr_string;

title(set_title)
ylabel('Estimated');
xlabel('Simulated');

h=refline(1,0);
h.LineStyle="--";
h.LineWidth=1.2;
h.Color="r";

f = gcf;
if ~isempty(varargin)
    saveas(f, varargin{1})
end

end