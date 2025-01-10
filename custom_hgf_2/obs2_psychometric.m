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

% Get x2
% Assumed structure of infStates:
% dim 1: time (ie, input sequence number)
% dim 2: HGF level
% dim 3: 1: muhat, 2: sahat, 3: mu, 4: sa
% I'm using mu (not muhat), following logrt model.
x2 = infStates(:,2,3); % use second level

% Predict response using psychometric function
gamma = 0;
lambda = 0;
beta = 3; %%%% Check if this is sensible
for i = 1:numel(x2)
    alpha = b0 + (b1*x2(i)); % regression for PSE
    morph = # % %sad

    % probability of response = 1 (using Weibull function)
    p_1 = gamma + (1 - gamma - lambda) * (1 - exp(-(morph/alpha).^beta));
    
    % predicted response
    resp = rand < p_1;
end




end