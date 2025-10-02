function [y, logrt] = obs3_gaussianRT_sim(r, infStates, p)
% [y, logrt] = m1_logrt_linear_binary_sim(r, infStates, p)
%
% Simulates logRTs with Gaussian noise.
% (Designed to be compatible with the HGF Toolbox as part of TAPAS).
%
% INPUT
%   r             struct      Struct obtained from tapas_simModel.m
%   infStates     tensor      Tensor containing inferred states from the
%                             perceptual model    
%   p             vector      1xP vector with free param values (nat space)
%
%   OPTIONAL:
%
% OUTPUT    
%   y             vector       Nx1 vector with simulated logRTs
%   logrt         vector       Nx1 vector containing noise-free predictions


% Parameters
mu = p(1);
sigma0 = p(2);
gamma = p(3);
lambda = p(4);
sigma1 = p(5);


% gaussian function
G = @(x,mu,sigma0,gamma,lambda)gamma+(lambda*exp(-((x-mu)^2)/(2*(sigma0^2))));

% inputs
u = r.u;

n = numel(u);


% Predict logRT
logrt=nan(size(u));
for i = 1:numel(u)
    logrt(i) = G(u(i), mu, sigma0, gamma, lambda);
end


% Initialize random number generator
if isnan(r.c_sim.seed)
    rng('shuffle');
else
    rng(r.c_sim.seed);
end

% Simulate (with noise)
y = logrt+sqrt(sigma1)*randn(n, 1);

end
