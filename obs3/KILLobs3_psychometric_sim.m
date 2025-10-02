function [y, prob] = obs3_psychometric_sim(r, infStates, p)
% simulates responses from psychometric function


% Transform parameters into native space
alpha = p(1);
beta = p(2);
gamma = p(3);
lambda = p(4);

% psychometric function
F = @(x,alpha,beta,gamma,lambda) gamma + (1 - gamma - lambda) ./ (1 + exp(-beta * (x - alpha)));

% inputs
u = r.u;

% loop through trials
prob = nan(size(u));
for i=1:numel(u)
    prob(i) = F(u(i), alpha, beta, gamma, lambda);
end

% Initialize random number generator
if isnan(r.c_sim.seed)
    rng('shuffle');
else
    rng(r.c_sim.seed);
end

% Simulate responses
y = binornd(1, prob);

end
