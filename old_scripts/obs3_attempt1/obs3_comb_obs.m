function [logp, yhat, res, logp_split] = obs3_comb_obs(r, infStates, ptrans)
% [logp, yhat, res, logp_split] = obs3_comb_obs(r, infStates, ptrans)
%
% Calculates the combined log-probability of binary and continuous
% behavioural responses (NO LEARNING)
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
%                              probabilities of responses
%   yhat          matrix       Nx2 matrix containing noise-free predictions
%   res           matrix       Nx2 matrix containing responses (bin + cont)
%   logp_split    matrix       Nx2 matrix containing trialwise logll values
%


%% Separate parameters
ptrans_psychometric = ptrans(1:4); % parameters for the psychometric model
ptrans_gaussian = ptrans(5:9); % parameters for the gaussianRT model


%% binary part of the response model

% compute log likelihood (binary responses)
[logp_binary, yhat_binary, res_binary] = ...
    obs3_psychometric(r, infStates, ptrans_psychometric);

%% continuous part of the response model

% prepare inputs for logRT GLM assuming that the input is a matrix with 2
% columns (1: binary predictions, 2: log response times).
% rt = r;
% rt.y = r.y(:,2);

% Compute the log likelihood (logRTs)
[logp_reactionTime, yhat_reactionTime, res_reactionTime] = ...
    obs3_gaussianRT(r, infStates, ptrans_gaussian);


%% confidence part of the response model
% confidence01 = r.y(:,4);
% 
% [logp_confidence, yhat_confidence, res_confidence] = ...
%     tapas_softmax_binary(confidence01, infStates, ptrans);

%% get combined log likelihood of two response data modalities
logp = logp_binary + logp_reactionTime;
logp_split = [logp_binary logp_reactionTime];

%% return predictions and responses for each response data modality
yhat = [yhat_binary yhat_reactionTime];
res = [res_binary res_reactionTime];

end
