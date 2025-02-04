function c = obs3_psychometric_config

% parameters for linear regression of volatility onto slope parameter. 
% B0 describes baseline slope, B1 describes influence of volatility

% Config structure
c = struct;

% Is the decision based on predictions or posteriors? Comment as appropriate.
% c.predorpost = 1; % Predictions (I think)?
c.predorpost = 2; % Posteriors

% Model name
c.model = 'obs3_psychometric_config';

% Sufficient statistics of Gaussian parameter priors
% B0
c.b0mu = 2.5;
c.b0sa = 2;

% B1
c.b1mu = 0;
c.b1sa = 2;

% Alpha (PSE)
c.alphamu = .5;
c.alphasa = 1;


% Gather prior settings in vectors
c.priormus = [
    c.b0mu,...
    c.b1mu,...
    c.alphamu,...
         ];

c.priorsas = [
    c.b0sa,...
    c.b1sa,...
    c.alphasa,...
         ];


% Model filehandle
c.obs_fun = @obs3_psychometric;

% Handle to function that transforms observation parameters to their native space
% from the space they are estimated in
c.transp_obs_fun = @obs3_psychometric_transp;

end