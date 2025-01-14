% Perceptual = binary HGF; response = binary (Tom's psychometric VIP)

% prc model is normal binary HGF, but edited to ignore column 2 of u


close all; clear;
addpath('custom_hgf_2');

% example data (to get contingencies etc)
sub_data = readtable('STE_data\10369536_A_Threat.csv');

% input in contingency space
state = double(sub_data.Cue_idx == sub_data.Outcome_idx);

% stimulus intensity
sub_data.p_sad = sub_data.Outcome_p_sad/100;

% model input
u = [state, sub_data.p_sad];

u_sub = u;



%% simulate responses

prc_model_config = prc2_hgf_binary_config(); % perceptual model
obs_model_config = obs2_psychometric_config(); % response model
% obs_model_config = tapas_unitsq_sgm_config(); % response model
optim_config     = tapas_quasinewton_optim_config(); % optimisation algorithm
optim_config.nRandInit = 5;

r_temp = [];
r_temp.c_prc.n_levels = 3;
prc_params = prc2_hgf_binary_transp(r_temp, prc_model_config.priormus);
obs_params = obs_model_config.priormus;


sim = tapas_simModel(u_sub,...
    'prc2_hgf_binary',...
    prc_params,...
    'obs2_psychometric',...
    obs_params,...
    123456789);

% visualise psychometric
sim_psychometric = arrayfun(@(x) mean(sim.y(sub_data.p_sad==x, 1)), 0:.2:1);
% figure('name', 'simulated psychometric'); hold on;
plot(0:.2:1, sim_psychometric, 'linewidth', 3);
set(gca, 'Ylim', [0,1], 'Xtick', 0:.2:1)


%% recover parameters
prc_model_config = prc2_hgf_binary_config(); % perceptual model
obs_model_config = obs2_psychometric_config(); % response model
% obs_model_config = tapas_unitsq_sgm_config(); % response model

optim_config.nRandInit = 5;

% adjust prior variance
prc_model_config.omsa(2)=1;
prc_model_config.priorsas(13)=prc_model_config.omsa(2);
prc_model_config.omsa(3)=1;
prc_model_config.priorsas(14)=prc_model_config.omsa(3);


% do estimate regression coefficients
obs_model_config.b0sa=2;
obs_model_config.b1sa=2;
obs_model_config.priorsas(1)=obs_model_config.b0sa;
obs_model_config.priorsas(2)=obs_model_config.b1sa;

est = tapas_fitModel(...
    sim.y,...
    sim.u,...
    prc_model_config,...
    obs_model_config,...
    optim_config);

tapas_fit_plotCorr(est)
tapas_hgf_binary_plotTraj(est)


%% full parameter recovery

om2_range = [-4 -1];
om3_range = [-8 -3];
b0_range = [.1 .9];
b1_range = [-.5, .5];

N=100;

om2_sim=nan(N,1);
om2_est=nan(N,1);
om3_sim=nan(N,1);
om3_est=nan(N,1);
b0_sim=nan(N,1);
b0_est=nan(N,1);
b1_sim=nan(N,1);
b1_est=nan(N,1);

for i=1:N
    % get parameters to simulate
    om2 = om2_range(1) + (om2_range(2)-om2_range(1))*rand;
    om3 = om3_range(1) + (om3_range(2)-om3_range(1))*rand;
    b0 = b0_range(1) + (b0_range(2)-b0_range(1))*rand;
    b1 = b1_range(1) + (b1_range(2)-b1_range(1))*rand;
    
    % parameter vectors
    prc_params(13) = om2;
    prc_params(14) = om3;
    obs_params = [b0, b1];
    
    % simulate
    sim = tapas_simModel(u,...
        'prc2_hgf_binary',...
        prc_params,...
        'obs2_psychometric',...
        obs_params);
    
    % fit to simulated
    est = tapas_fitModel(...
        sim.y,...
        sim.u,...
        prc_model_config,...
        obs_model_config,...
        optim_config);
    
    % get estimated parameters
    om2_sim(i) = om2;
    om3_sim(i) = om3;
    b0_sim(i) = b0;
    b1_sim(i) = b1;
    om2_est(i) = est.p_prc.om(2);
    om3_est(i) = est.p_prc.om(3);
    b0_est(i) = est.p_obs.b0;
    b1_est(i) = est.p_obs.b1;

end




%% fit real data 

close all;

prc_model_config.omsa(3)=0;
prc_model_config.priorsas(14)=prc_model_config.omsa(3);

obs_model_config.b0sa=1;
obs_model_config.b1sa=1;
obs_model_config.priorsas(1)=obs_model_config.b0sa;
obs_model_config.priorsas(2)=obs_model_config.b1sa;

optim_config     = tapas_quasinewton_optim_config(); % optimisation algorithm
optim_config.nRandInit = 20;



data_dir=[pwd,'\STE_data\'];
STE_files = dir(fullfile(data_dir, '*.csv'));
model_fits = repmat(struct(), numel(STE_files), 1);
LMEs = nan(numel(STE_files), 1);

for i_file=1:numel(STE_files)
    fname = STE_files(i_file).name;
    fprintf('\n%s\n', fname);
    sub_data = readtable(fullfile(data_dir, fname));
    fname_tokens = regexp(fname, '(\d+)_(\w)_(\w+)\.csv', 'tokens');
    [model_fits(i_file).ID, model_fits(i_file).group, model_fits(i_file).condition] = fname_tokens{:}{:};

    % Model input
    state = double(sub_data.Cue_idx == sub_data.Outcome_idx); % input in contingency space
    sub_data.p_sad = sub_data.Outcome_p_sad/100; % stimulus intensity
    u_sub = [state, sub_data.p_sad]; % model input
    
    % response
    y_sub = sub_data.Response_idx; % binary response
    
    % remove missing trials ----- Need to sort this out when finished -
    % better to have ignored trials for averaging
    nan_trials = isnan(y_sub);
    u_sub = u_sub(~nan_trials, :);
    y_sub = y_sub(~nan_trials);

    % add to model_fits
    model_fits(i_file).y = y_sub;

    % fit data
    try
        est = tapas_fitModel(...
            y_sub,... 
            u_sub,...
            prc_model_config,...
            obs_model_config,...
            optim_config);
        
        model_fits(i_file).est = est;
        LMEs(i_file) = est.optim.LME;
    catch
        model_fits(i_file).est = [];
    end
end

save('STE_model4_LME.mat', 'LMEs');









