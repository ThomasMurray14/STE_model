function [logp, yhat, res] = obs3_psychometric(r, infStates, ptrans)
% Calculates the log-probability of response y=1

% Transform parameters into native space
alpha = tapas_sgm(ptrans(1), 1); % PSE
beta = exp(ptrans(2)); % slope
gamma = tapas_sgm(ptrans(3), 1); % upper asymptote
lambda = tapas_sgm(ptrans(4), 1); % lower asymptote

% psychometric function
F = @(x,alpha,beta,gamma,lambda) gamma + (1 - gamma - lambda) ./ (1 + exp(-beta * (x - alpha)));

% inputs
u = r.u;

% responses
y = r.y(:,1);

% probabilities
yhat=nan(size(u));
logp=nan(size(u));

for i = 1:numel(u)
    % p(response==1)
    p_1 = F(u(i),alpha,beta,gamma,lambda);

    % predicted outcome
    yhat(i)=p_1;

    % probability of observed outcome
    logp(i) = r.y(i) * log(p_1) + (1 - r.y(i)) * log(1 - p_1);

end

% Calculate residuals (standardized difference between observed and predicted)
res = (y - yhat) ./ sqrt(yhat .* (1 - yhat));

return;
