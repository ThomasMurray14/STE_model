function c = obs2_psychometric_config

% parameters for linear regression of volatility onto PSE. B0 describes
% baseline position of PSE, B1 describes influence of volatility

% Config structure
c = struct;

% Is the decision based on predictions or posteriors? Comment as appropriate.
% c.predorpost = 1; % Predictions (I think)?
c.predorpost = 2; % Posteriors

% Model name
c.model = 'obs2_psychometric_config';

% Sufficient statistics of Gaussian parameter priors
% B0
c.b0mu = .5;
c.b0sa = 1;

% B1
c.b1mu = 0;
c.b1sa = 1;

% Beta (psychometric function slope)
c.logbetamu = log(3);
c.logbetasa = 2;

% Gather prior settings in vectors
c.priormus = [
    c.b0mu,...
    c.b1mu,...
    c.logbetamu,...
         ];

c.priorsas = [
    c.b0sa,...
    c.b1sa,...
    c.logbetasa,...
         ];


% Model filehandle
c.obs_fun = @obs2_psychometric;

% Handle to function that transforms observation parameters to their native space
% from the space they are estimated in
c.transp_obs_fun = @obs2_psychometric_transp;

end