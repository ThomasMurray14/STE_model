function [prc_config, obs_config] = STE_SuttonK1_config


% Perceptual model
prc_config = prc_sutton_k1_binary_config();
prc_config.logmumu        = log(3); % set mu prior mean
prc_config.logmusa        = 1; % set mu prior sa
prc_config.logRhatmu      = log(1); % set prior mean Rhat
prc_config.logRhatsa      = 0; % set prior mean Rhat
prc_config.logitvhat_1mu  = tapas_logit(0.5, 1); % set vhat1 to .5
prc_config.logitvhat_1sa  = 0; % fix
prc_config.logh_1mu       = log(.005);
prc_config.logh_1sa       = 0; % fix
prc_config = tapas_align_priors(prc_config);


% Response model
obs_config = obs_suttonK1_comb_obs_config();
obs_config.logzemu = log(1);
obs_config.logzesa = 2;
obs_config.beta0mu = 6.5;
obs_config.beta0sa = 4;
obs_config.beta1mu = -1;
obs_config.beta1sa = 4;
obs_config.beta2mu = 4;
obs_config.beta2sa = 2;
obs_config.logsamu = log(0.1);
obs_config.logsasa = 2;
obs_config = tapas_align_priors(obs_config);


end