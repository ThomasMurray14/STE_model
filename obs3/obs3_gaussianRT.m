function [logp, yhat, res] = obs3_gaussianRT(r, infStates, ptrans)
% [logp, yhat, res] = m1_logrt_linear_binary(r, infStates, ptrans)
%
% Calculates the log-probability of log-reaction times y 
%
% INPUT
%   r             struct      Struct obtained from tapas_fitModel.m fct
%   infStates     tensor      Tensor containing inferred states from the
%                             perceptual model    
%   ptrans        vector      1xP vector with free param values (est space)
%
%   OPTIONAL:
%
% OUTPUT    
%   logp          vector       1xN vector containing trialwise log
%                              probabilities of logRTs
%   yhat          vector       1xN vector containing noise-free predictions
%   res           vector       1xN vector containing residuals
%



% Transform parameters to their native space
mu = tapas_sgm(ptrans(1), 1);
sigma0 = exp(ptrans(2));
gamma = exp(ptrans(3));
lambda = exp(ptrans(4));
sigma1 = exp(ptrans(5));



% gaussian function
G = @(x,mu,sigma0,gamma,lambda)gamma+(lambda*exp(-((x-mu)^2)/(2*(sigma0^2))));

% inputs
u = r.u;

% responses
y = r.y(:,2);

% probabilities
yhat=nan(size(u));
logp=nan(size(u));
res =nan(size(u));

% Predict logRT
logrt=nan(size(u));
for i = 1:numel(u)
    logrt(i) = G(u(i), mu, sigma0, gamma, lambda);
end


% Calculate log-probabilities (taken straight from original tapas model)
% 
% Note: 8*atan(1) == 2*pi (this is used to guard against
% errors resulting from having used pi as a variable).
reg = ~ismember(1:n,r.irr);
logp(reg) = -1/2.*log(8*atan(1).*sigma1) -(y-logrt).^2./(2.*sigma1);
yhat(reg) = logrt;
res(reg) = y-logrt;



return;
