function y = obs2_psychometric_sim(r, infStates, p)


% Get parameters
b0  = p(1);
b1  = p(2);

% Number of trials
n = size(infStates,1);

% Stimulus intensity
intensity = r.u(:,2);

% get volatility
mu3 = infStates(:,3,3); % use third level

% estimate trialwise PSE
PSE = b0 + b1.*exp(mu3);

% Predict response using psychometric function
gamma = 0;
lambda = 0;
beta = 3; % Shape parameter (fixed for now)

prob = nan(n, 1);
for i = 1:n
    alpha = PSE(i); % regression for PSE

    stim_intensity = intensity(i); % get %sad

    % probability of response = 1 (using Weibull function)
    prob(i) = gamma + (1 - gamma - lambda) * (1 - exp(-(stim_intensity/alpha).^beta));

end

% Simulate
y = binornd(1, prob);


end