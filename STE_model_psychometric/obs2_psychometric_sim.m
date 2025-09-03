function y = obs2_psychometric_sim(r, infStates, p)


% Get parameters
b0  = p(1);
b1  = p(2);
zeta = p(3); % shape parameter

% Number of trials
n = size(infStates,1);

% Stimulus intensity
intensity = r.u(:,2);

% get volatility
mu3 = infStates(:,3,3); % use third level

% estimate trialwise PSE
PSE = b0 + b1.*exp(mu3);

% set gamma/lambda
gamma = 0;
lambda = 0;

prob = nan(n, 1);
for i = 1:n

    alpha = PSE(i); % regression for PSE
    x = intensity(i); % get %sad

    % probability of response = 1 (using logistic)
    prob(i) = gamma + (1 - gamma - lambda).*(1./(1+exp(-1*(zeta).*(x-alpha))));

end

% Simulate
y = binornd(1, prob);



end