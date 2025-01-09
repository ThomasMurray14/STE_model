function obs2_psychometric(r, infStates, ptrans)
% Tom's attempt at transforming the inferred perceptual states into
% psychometric function parameters

be0  = ptrans(1);
be1  = ptrans(2);

% Initialize returned log-probabilities, predictions,
% and residuals as NaNs so that NaN is returned for all
% irregualar trials
n = size(infStates,1);
logp = NaN(n,1);
yhat = NaN(n,1);
res  = NaN(n,1);

% % Weed irregular trials out from responses and inputs
% y = r.y(:,1);
% y(r.irr) = [];
% u = r.u(:,1);
% u(r.irr) = [];




end