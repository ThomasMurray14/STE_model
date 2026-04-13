% Script to visualise parameter correlations


model = 'ResponseBias';
fits = importdata(['..\model_', model, '\model_', model, '_fit.mat']);

C = nan([size(fits(1).est.optim.Corr), numel(fits)]);
for i = 1:numel(fits)
    C(:,:,i) = fits(i).est.optim.Corr;
end

mean_C = mean(C, 3, 'omitmissing');

r = fits(1).est;

% Nicked from hgf toolbox:
% Determine indices of parameters to optimize (i.e., those that are not fixed and not NaN)
prc_ind = r.c_prc.priorsas;
prc_ind(isnan(prc_ind)) = 0;
prc_ind = find(prc_ind);

obs_ind = r.c_obs.priorsas;
obs_ind(isnan(obs_ind)) = 0;
obs_ind = find(obs_ind);

n_par   = length(prc_ind) + length(obs_ind);

% Find names of optimized parameters to use them as tick labels 
names_prc = fieldnames(r.p_prc);
fields = struct2cell(r.p_prc);
expnms_prc = [];
for k = 1:length(names_prc)
    for l= 1:length(fields{k})
    expnms_prc = [expnms_prc, names_prc(k)];
    end
end
expnms_prc = expnms_prc(1:length(r.p_prc.p))';

names_obs = fieldnames(r.p_obs);
fields = struct2cell(r.p_obs);
expnms_obs = [];
for k = 1:length(names_obs)
    for l= 1:length(fields{k})
    expnms_obs = [expnms_obs, names_obs(k)];
    end
end
expnms_obs = expnms_obs(1:length(r.p_obs.p))';

ticklabels = {[expnms_prc(prc_ind); expnms_obs(obs_ind)]};

% Plot matrix
imagesc(mean_C, [-1 1]);
axis('square');
set(gca,'xtick',1:n_par);
set(gca,'ytick',1:n_par);
set(gca,'xticklabel',ticklabels{:});
set(gca,'yticklabel',ticklabels{:});
colorbar;
title('Average parameter correlation', 'FontSize', 15, 'FontWeight', 'bold');
