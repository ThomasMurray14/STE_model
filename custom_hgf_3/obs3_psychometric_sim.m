function y = obs3_psychometric_sim(r, infStates, p)

u = r.u(:,1);

% get parameters
b0 = p(1);
b1 = p(2);
b2 = p(3);
b3 = p(4);
b4 = p(5);

% Number of trials
n = size(infStates,1);

% Stimulus intensity
intensity = r.u(:,2);

% Extract trajectories of interest from infStates
mu1hat = infStates(:,1,1);
sa1hat = infStates(:,1,2);
mu2    = infStates(:,2,3);
sa2    = infStates(:,2,4);
mu3    = infStates(:,3,3);

% Surprise
% ~~~~~~~~
m1hreg = mu1hat;
poo = m1hreg.^u.*(1-m1hreg).^(1-u); % probability of observed outcome
surp = -log2(poo);

% Bernoulli variance (aka irreducible uncertainty, risk) 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bernv = sa1hat;

% Inferential variance (aka informational or estimation uncertainty, ambiguity)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
inferv = tapas_sgm(mu2, 1).*(1 -tapas_sgm(mu2, 1)).*sa2; % transform down to 1st level

% Phasic volatility (aka environmental or unexpected uncertainty)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pv = tapas_sgm(mu2, 1).*(1-tapas_sgm(mu2, 1)).*exp(mu3); % transform down to 1st level

% Calculate predicted PSE
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PSE = b0 +b1.*surp +b2.*bernv +b3.*inferv +b4.*pv;

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