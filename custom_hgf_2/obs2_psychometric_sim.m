function y = obs2_psychometric_sim(r, infStates, p)


% Get parameters
b0  = p(1);
b1  = p(2);

% Number of trials
n = size(infStates,1);

% Stimulus intensity
intensity = r.u(:,2);

% Get x3
x3 = infStates(:,3,3); % use second level

% Predict response using psychometric function
gamma = 0;
lambda = 0;
beta = 3; % Shape parameter (fixed for now)

prob = nan(n, 1);
for i = 1:numel(x3)
    alpha = b0 + (b1*x3(i)); % regression for PSE

    stim_intensity = intensity(i); % get %sad

    % probability of response = 1 (using Weibull function)
    prob(i) = gamma + (1 - gamma - lambda) * (1 - exp(-(stim_intensity/alpha).^beta));

end

% Simulate
y = binornd(1, prob);


end