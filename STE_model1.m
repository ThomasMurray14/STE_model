%% TO DO

% Change the response model:
% try merging correct/incorrect + logrt + confidence

% can we add mu3 to explain RTs?
% split corr/incorr or 
% tapas_ehgf_binary_combObs_plotTraj
% bias of happy vs sad (perception or response)

%%

close all; clear;
% addpath([pwd, '\custom_hgf']);


%% 

% example data (to get contingencies etc)
sub_data = readtable('STE_data\10369536_A_Safe.csv');

% Contingency space
sub_data.p_sad = sub_data.Outcome_p_sad/100;
state = double(sub_data.Cue_idx == sub_data.Outcome_idx);
u_al = sub_data.p_sad;
tr_temp = or((state ==1 & sub_data.Outcome_idx== 0), (state == 0 & sub_data.Outcome_idx == 1));
u_al(tr_temp) = 1 - sub_data.p_sad(tr_temp);

sub_data.u_al = u_al;
sub_data.state=state;


% Responses
sub_data.logRT = log(sub_data.Response_RT);
sub_data.correct = sub_data.Outcome_idx == sub_data.Response_idx;

% sub_data = sub_data(~isnan(sub_data.logRT),:); % remove nans 
%happy_face - 0.5

u = [sub_data.u_al, sub_data.Cue_idx];
y = [sub_data.correct,sub_data.logRT];

u_sub = u(~isnan(sub_data.Response_idx),:);
y_sub = y(~isnan(sub_data.Response_idx),:);


%% Get configuration structures
prc_model_config = prc1_ehgf_binary_pu_tbt_config(); % perceptual model
obs_model_config = obs1_comb_obs_config(); % response model
optim_config     = tapas_quasinewton_optim_config(); % optimisation algorithm
optim_config.nRandInit = 5;


%% simulate responses

r_temp = [];
r_temp.c_prc.n_levels = 3;
prc_params = prc1_ehgf_binary_pu_tbt_transp(r_temp, prc_model_config.priormus);

prc_params(1); %mu_0mu(1)
prc_params(2); %mu_0mu(2)
prc_params(3); %mu_0mu(3)
prc_params(4); %logsa_0mu(1);
prc_params(5); %logsa_0mu(2);
prc_params(6); %logsa_0mu(3);
prc_params(7); %rho(1);
prc_params(8) = 0; %rho(2);
prc_params(9); %rho(3);
prc_params(10); %logkamu(1);
prc_params(11); %logkamu(2);
prc_params(12); %ommu(1);
prc_params(13); %ommu(2);
prc_params(14); %ommu(3);
prc_params(15) = 0.002; %logalmu;
prc_params(16); %eta0mu;
prc_params(17); %eta1mu;


obs_params = obs_model_config.priormus;
obs_params(1) = exp(obs_params(1)); %%%%%%%%%%%% Why do I do this bit?
obs_params(7) = exp(obs_params(7));

sim = tapas_simModel(u_sub,...
    'prc1_ehgf_binary_pu_tbt',...
    prc_params,...
    'obs1_comb_obs',...
    obs_params,...
    123456789);



%% plot
% prc1_ehgf_binary_tbt_plotTraj(sim);

% Visualise Psychometric
% sim_psychometric = arrayfun(@(x) mean(sim.y(sub_data.Outcome_p_sad==x, 1)), 0:20:100);
% figure('name', 'simulated psychometric'); hold on;
% plot(0:20:100, sim_psychometric, 'linewidth', 3);
% set(gca, 'Ylim', [0,1], 'Xtick', 0:20:100)



%% recover parameters

prc_model_config = prc1_ehgf_binary_pu_tbt_config(); % perceptual model

est = tapas_fitModel(...
    sim.y,...
    sim.u,...
    prc_model_config,...
    obs_model_config,...
    optim_config);

% Check parameter identifiability
tapas_fit_plotCorr(est)
prc1_ehgf_binary_tbt_plotTraj(est)


%% fit real data 

close all;

prc_model_config = prc1_ehgf_binary_pu_tbt_config(); % perceptual model %%%%% Do I need to keep loading this?
obs_model_config = obs1_comb_obs_config();%tapas_logrt_linear_binary_config(); % response model
optim_config     = tapas_quasinewton_optim_config(); % optimisation algorithm
optim_config.nRandInit = 5;

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

    % prepare data
    sub_data.p_sad = sub_data.Outcome_p_sad/100;
    sub_data.u_al = u_al; % these are the same across everyone
    sub_data.state = state;

    sub_data.logRT = log(sub_data.Response_RT);
    sub_data = sub_data(~isnan(sub_data.logRT),:);    % remove nans 
    sub_data.correct = sub_data.Outcome_idx == sub_data.Response_idx;
    
    y_sub = [sub_data.correct, sub_data.logRT];
    u_sub = u(~isnan(sub_data.logRT),:);
    model_fits(i_file).y = y_sub;
    
    % fit data
    try
        est = tapas_fitModel(...
            y_sub,... %,confidence],...
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

save('STE_model1_LME.mat', 'LMEs');




%%
sum(cellfun(@isempty, {model_fits.est})) % N with no fit


%% Visualise psychometric

i_sub=1;
sub_data = readtable(fullfile(data_dir, STE_files(i_sub).name));
sub_fit = model_fits(i_sub).est;

p_sad = [0, 20, 40, 60, 80, 100];
p_sad_resp = arrayfun(@(x) mean(sub_data.Response_idx(sub_data.Outcome_p_sad == x), 'omitnan'), p_sad);

% simulate using fit parameters
N_sim=30;
sim_resp = zeros(N_sim, 6);
for i_sim = 1:30
    fit_prc_params = sub_fit.p_prc.p;
    fit_obs_params = sub_fit.p_obs.p;
    
    sim = tapas_simModel(u,...
        'prc1_ehgf_binary_pu_tbt',...
        fit_prc_params,...
        'obs1_comb_obs',...
        fit_obs_params);

    % simulated responses
    sim_resp(i_sim, :) = arrayfun(@(x) mean(sim.y(sub_data.Outcome_p_sad==x, 1)), p_sad);

end


figure; hold on;
plot(p_sad, p_sad_resp, 'linewidth', 3);
plot(p_sad, sim_resp, '--')
set(gca, 'Ylim', [0, 1])











%% Write fit parameters to csv file
% IDs = unique(cellfun(@str2double, {model_fits.ID}));
% conds = {'Safe', 'Threat'};
% var_names = {'ID', 'group', 'S_om2', 'S_al', 'S_rho', 'T_om2', 'T_al', 'T_rho'};
% var_types = {'double', 'string', 'double','double','double','double','double','double'};
% parameter_fits = table(...
%     'Size', [numel(IDs), numel(var_names)],...
%     'VariableNames', var_names,...
%     'VariableTypes', var_types);
% 
% for i=1:numel(IDs)
%     sub_fits = model_fits(strcmp(num2str(IDs(i)), {model_fits.ID}));
%     parameter_fits(i, 'ID') = {IDs(i)};
%     parameter_fits(i, 'group') = {sub_fits(1).group};
%     for c=1:2
%         cond_fit = sub_fits(strcmp({sub_fits.condition}, conds{c}));
%         parameter_fits(i, [conds{c}(1), '_om2']) = {cond_fit.est.p_prc.om(2)};
%         parameter_fits(i, [conds{c}(1), '_al']) = {cond_fit.est.p_prc.al};
%         parameter_fits(i, [conds{c}(1), '_rho']) = {cond_fit.est.p_prc.rho(2)};
%     end
% end
% 
% % writetable(parameter_fits, 'parameter_fits.csv');





%%


% Check parameter identifiability
% tapas_fit_plotCorr(est)

% tapas_ehgf_binary_tbt_plotTraj(est)
% 
% figure;histogram(est.optim.yhat(:,1),50)
% 
% figure;plot(1:length(est.optim.yhat(:,1)), est.optim.yhat(:,1))
% 
% % figure;histogram(log(RT(2:end))) %% Doesn't work
% figure;plot(sub_data.logRT, est.optim.yhat(:,2), '.') 
% close all
% 
%% check behavior
% sub_data_selected = sub_data;
% sub_data_selected.muhat1_state = est.traj.muhat(:,1);
% sub_data_selected.muhat2_state = est.traj.muhat(:,2);
% sub_data_selected.muhat2_corr = sub_data_selected.muhat2_state;
% sub_data_selected.muhat2_corr(sub_data_selected.state == 0) = -sub_data_selected.muhat2_state(sub_data_selected.state == 0);
% sub_data_selected.sahat1 = est.traj.sahat(:,1);
% sub_data_selected.sahat2 = est.traj.sahat(:,2);
% 
% sub_data_selected.logRT = sub_data_selected.logRT - nanmean(sub_data_selected.logRT);
% 
% % corr(table2array(sub_data_selected(2:end,[13,9,10,11,12,14]))) %% Doesn't work
% sub_data_selected.Expectedness01 = strcmp(sub_data_selected.Expectedness,'E');
% lm = fitlm(sub_data_selected, 'logRT ~muhat2_corr');
% 
% % Display model summary
% disp(lm);
% figure;
% hold on
% x = sub_data_selected.muhat2_corr;%muhat2_state;
% plot(strcmp(sub_data_selected.Expectedness,'UE'), sub_data_selected.muhat2_corr, '.')
% plot(x, sub_data_selected.logRT, '.')
% plot(x(strcmp(sub_data_selected.Expectedness,'UE')), sub_data_selected.logRT(strcmp(sub_data_selected.Expectedness,'UE')), 'o')
% 
% %% recover pars
% prc_params = est.p_prc.p;
% obs_params = est.p_obs.p;
% 
% sim = tapas_simModel(u,...
%     'tapas_ehgf_binary_pu_tbt',...
%     prc_params,...
%     'm1_comb_obs',...
%     obs_params,...
%     123456789);
% est_sim = tapas_fitModel(...
%     sim.y,...
%     sim.u,...
%     prc_model_config,...
%     obs_model_config,...
%     optim_config);
% 
% prc_params = est.p_prc.p;
% obs_params = est.p_obs.p;
% 
% prc_params_sim = est_sim.p_prc.p;
% obs_params_sim = est_sim.p_obs.p;
% 

