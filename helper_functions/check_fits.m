% script to check parameter fits between fit and fit2 (high ze)

model_name = 'PerceptualBias';

fits = importdata(['..\model_', model_name, '\model_', model_name, '_fit.mat']);
fits2 = importdata(['..\model_', model_name, '\model_', model_name, '_fit2.mat']);


params = {
    'om2',   @(s) s.est.p_prc.om(2)
    'rho',   @(s) s.est.p_prc.rho(2)
    'al',    @(s) s.est.p_prc.al
    'zeta',  @(s) s.est.p_obs.ze
    'beta0', @(s) s.est.p_obs.beta0
    'beta1', @(s) s.est.p_obs.beta1
    'beta2', @(s) s.est.p_obs.beta2
    'beta3', @(s) s.est.p_obs.beta3
    'beta4', @(s) s.est.p_obs.beta4
    'sa',    @(s) s.est.p_obs.sa
};

for i = 1:size(params,1)
    name = params{i,1};
    fn   = params{i,2};

    a = arrayfun(@(x) safeget(fits(x), fn), 1:400);
    b = arrayfun(@(x) safeget(fits2(x), fn), 1:400);

    plot_corr(a,b,name)
end

function val = safeget(s, path)
    if isempty(s.est)
        val = NaN;
    else
        val = path(s);
    end
end

function plot_corr(a,b,name)
    idx = ~isnan(a) & ~isnan(b);
    figure('name', name)
    scatter(a(idx), b(idx))
    [r,p] = corr(a(idx)', b(idx)');
    title(sprintf('%s: r = %.2f, p = %.3f', name, r, p))
end



