function c = obs3_psychometric_config

% parameters for linear regression of volatility onto PSE. B0 describes
% baseline position of PSE, B1 describes influence of volatility

% Config structure
c = struct;

% Is the decision based on predictions or posteriors? Comment as appropriate.
% c.predorpost = 1; % Predictions (I think)?
c.predorpost = 2; % Posteriors

% Model name
c.model = 'obs3_psychometric_config';

% Sufficient statistics of Gaussian parameter priors
% B0
c.b0mu = .5;
c.b0sa = .2;

% B1
c.b1mu = 0;
c.b1sa = 0;

% Beta_2
c.b2mu = 0; 
c.b2sa = 0;

% Beta_3
c.b3mu = 0; 
c.b3sa = 2;

% Beta_4
c.b4mu = 0; 
c.b4sa = 2;

% Gather prior settings in vectors
c.priormus = [
    c.b0mu,...
    c.b1mu,...
    c.b2mu,...
    c.b3mu,...
    c.b4mu,...
         ];

c.priorsas = [
    c.b0sa,...
    c.b1sa,...
    c.b2sa,...
    c.b3sa,...
    c.b4sa,...
         ];

% Model filehandle
c.obs_fun = @obs3_psychometric;

% Handle to function that transforms observation parameters to their native space
% from the space they are estimated in
c.transp_obs_fun = @obs3_psychometric_transp;

end