function [logp, yhat, res] = obs2_psychometric(r, infStates, ptrans)
% Tom's attempt at transforming the inferred perceptual states into
% psychometric function parameters

% u(:,1) = contingency space
% u(:,2) = intensity (%stim morph)
% y = binary response (not in contingency space)

% get parameters
b0  = ptrans(1);
b1  = ptrans(2);

% Initialize returned log-probabilities, predictions,
% and residuals as NaNs so that NaN is returned for all
% irregualar trials
n = size(infStates,1);
logp = NaN(n,1);
yhat = NaN(n,1);

% Get x2
% Assumed structure of infStates:
% dim 1: time (ie, input sequence number)
% dim 2: HGF level
% dim 3: 1: muhat, 2: sahat, 3: mu, 4: sa
% I'm using mu (not muhat), following logrt model.
x3 = infStates(:,3,3); % use third level

% Get stimulus intensity
intensity = r.u(:,2);

% Predict response using psychometric function
gamma = 0;
lambda = 0;
beta = 3; % Shape parameter (fixed for now)

for i = 1:numel(x3)
    alpha = b0 + (b1*x3(i)); % regression for PSE

    stim_intensity = intensity(i); % get %sad

    % probability of response = 1 (using Weibull function)
    p_1 = gamma + (1 - gamma - lambda) * (1 - exp(-(stim_intensity/alpha).^beta));

    % Avoid numerical issues with probabilities near 0 or 1
    p_1 = max(min(p_1, 1 - 1e-16), 1e-16); % Might need to fix this....................

    % store probability
    yhat(i)=p_1;

    % Compute log probability of observed response
    if ~isnan(r.y(i)) % Ignore irregular trials
        logp(i) = r.y(i) * log(p_1) + (1 - r.y(i)) * log(1 - p_1);
    end

end

% Calculate residuals (standardized difference between observed and predicted)
res = (r.y - yhat) ./ sqrt(yhat .* (1 - yhat));



end