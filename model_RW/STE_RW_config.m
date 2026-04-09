function [prc_config, obs_config] = STE_RW_config

prc_config = prc_RW_binary_config();

prc_config.logitv_0mu = tapas_logit(.5, 1);
prc_config.logitv_0sa = 0;
prc_config.logitalmu = tapas_logit(.5, 1);
prc_config.logitalsa = 1;

obs_config = obs_RW_comb_obs_config();

obs_config.logzemu = log(1);
obs_config.logzesa = 2;
obs_config.beta0mu = 6.5;
obs_config.beta0sa = 4;
obs_config.beta1mu = 4;
obs_config.beta1sa = 4;
obs_config.beta2mu = 4;
obs_config.beta2sa = 4;
obs_config.logsamu = log(0.1);
obs_config.logsasa = 2;




end