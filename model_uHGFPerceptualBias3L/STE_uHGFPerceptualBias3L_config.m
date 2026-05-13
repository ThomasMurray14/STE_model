function [prc_config, obs_config] = STE_uHGFPerceptualBias3L_config

% Perceptual model
% prc_config = prc_PerceptualBias3L_ehgf_binary_pu_tbt_config;
prc_config = uhgf_binary_pu_config;

prc_config.model = 'prc_uHGFPerceptualBias3L_binary_pu_tbt';
prc_config.ommu(2)    = -2;
prc_config.omsa(2)    = 4;
prc_config.ommu(3)    = -2;
prc_config.omsa(3)    = 4;
prc_config.rhomu(2)   = 0; % bias towards sad
prc_config.rhosa(2)   = 4;
prc_config.logkamu(1) = log(1);
prc_config.logkasa(1) = 0;
prc_config.logkamu(2) = log(1); % make sure 3 level
prc_config.logkasa(2) = 0;
prc_config.logalmu    = log(.05); % perceptual uncertainty
prc_config.logalsa    = 2;
prc_config = align_priors(prc_config);
prc_config.prc_fun = @prc_uHGFPerceptualBias3L_binary_pu_tbt;
prc_config.transp_prc_fun = @prc_uHGFPerceptualBias3L_binary_pu_tbt_transp;


% Response model
obs_config = obs_uHGFPerceptualBias3L_comb_obs_config;

obs_config.logzemu = log(5);%log(1);
obs_config.logzesa = 2;
obs_config.beta0mu = 6.5000;
obs_config.beta0sa = 4;

obs_config.beta1mu = 0;
obs_config.beta1sa = 4;

obs_config.beta2mu = 0;
obs_config.beta2sa = 4;

obs_config.beta3mu = 0;
obs_config.beta3sa = 4;

obs_config.beta4mu = 4;
obs_config.beta4sa = 4;

obs_config.beta5mu = 0;
obs_config.beta5sa = 4;

obs_config.logsasa = log(.1);
obs_config.logsasa = 2;

obs_config = align_priors(obs_config);


end