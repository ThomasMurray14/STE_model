% Perceptual = eHGF (Nace version), Response = logRT

close all; clear;

% example data (to get contingencies etc)
sub_data = readtable('STE_data\10369536_A_Threat.csv');

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

u = [sub_data.u_al, sub_data.Cue_idx];
y = [sub_data.logRT];

u_sub = u(~isnan(sub_data.Response_idx),:);
y_sub = y(~isnan(sub_data.Response_idx),:);




%% models
prc_model_config = prc1_ehgf_binary_pu_tbt_config(); % perceptual model
obs_model_config = tapas_logrt_linear_binary_config(); % response model
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
prc_params(8) = 100; %rho(2);
prc_params(9); %rho(3);
prc_params(10)=.001; %logkamu(1);
prc_params(11); %logkamu(2);
prc_params(12); %ommu(1);
prc_params(13)=-1; %ommu(2);
prc_params(14)=-1; %ommu(3);
prc_params(15) = 0.002; %logalmu;
prc_params(16); %eta0mu;
prc_params(17); %eta1mu;


obs_model_config.be4sa = 0;
obs_model_config.priorsas(5) = obs_model_config.be4sa;

obs_params = obs_model_config.priormus;



sim = tapas_simModel(u_sub,...%[u;u;u;u],...
    'prc1_ehgf_binary_pu_tbt',...
    prc_params,...
    'tapas_logrt_linear_binary',...
    obs_params,...
    123456789);


%% plot
% prc1_ehgf_binary_tbt_plotTraj(sim);


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
tapas_hgf_binary_plotTraj(est)




%% fit real data 

close all;

prc_model_config = prc1_ehgf_binary_pu_tbt_config(); % perceptual model
obs_model_config = tapas_logrt_linear_binary_config(); % response model
obs_model_config.be4sa = 0; % don't estimate beta 4 (influence of top level)
obs_model_config.priorsas(5) = obs_model_config.be4sa;

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
    sub_data.correct = sub_data.Outcome_idx == sub_data.Response_idx;
    sub_data.logRT = log(sub_data.Response_RT);
    
    u_sub = u(~isnan(sub_data.Response_idx),:);
    y_sub = sub_data.logRT(~isnan(sub_data.Response_idx));
    
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

save('STE_model3_LME.mat', 'LMEs');




