function c = obs4_psychometric_config

% parameters for linear regression of volatility onto PSE. B0 describes
% baseline position of PSE, B1 describes influence of precision weighted
% surprise

% Config structure
c = struct;

% Is the decision based on predictions or posteriors? Comment as appropriate.
% c.predorpost = 1; % Predictions (I think)?
c.predorpost = 2; % Posteriors

% Model name
c.model = 'obs4_psychometric';

% Sufficient statistics of Gaussian parameter priors
% B0
c.logitb0mu = tapas_logit(0.5, 1);
c.logitb0sa = 2;

% B1
c.b1mu = 0;
c.b1sa = 2;

% Zeta (psychometric function slope)
c.logitzetamu = tapas_logit(0.5, 1); % scaled to 0-100 in model
c.logitzetasa = 4;


% Gather prior settings in vectors
c.priormus = [
    c.logitb0mu,...
    c.b1mu,...
    c.logitzetamu,...
         ];

c.priorsas = [
    c.logitb0sa,...
    c.b1sa,...
    c.logitzetasa,...
         ];


% Model filehandle
c.obs_fun = @obs4_psychometric;

% Handle to function that transforms observation parameters to their native space
% from the space they are estimated in
c.transp_obs_fun = @obs4_psychometric_transp;

end