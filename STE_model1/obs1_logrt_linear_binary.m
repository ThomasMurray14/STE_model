function [logp, yhat, res] = obs1_logrt_linear_binary(r, infStates, ptrans)
% [logp, yhat, res] = m1_logrt_linear_binary(r, infStates, ptrans)
%
% Calculates the log-probability of log-reaction times y (in units of
% log-ms) according to the linear log-RT model developed for the SPIRL
% study based on the Lawson logRT GLM.
% (Designed to be compatible with the HGF Toolbox as part of TAPAS).
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
% _________________________________________________________________________
% Author: Alex Hess
%
% Copyright (C) 2023 Translational Neuromodeling Unit
%                    Institute for Biomedical Engineering
%                    University of Zurich & ETH Zurich
%
% This file is released under the terms of the GNU General Public Licence
% (GPL), version 3. You can redistribute it and/or modify it under the
% terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any
% later version.
%
% This file is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
% more details.
% 
% You should have received a copy of the GNU General Public License along
% with this program. If not, see <https://www.gnu.org/licenses/>.
% _________________________________________________________________________

% Transform parameters to their native space
be0  = ptrans(1);
be1  = ptrans(2);
be2  = ptrans(3);
be3  = ptrans(4);
be4  = ptrans(5);
sa   = exp(ptrans(6));

% Initialize returned log-probabilities, predictions,
% and residuals as NaNs so that NaN is returned for all
% irregualar trials
n = size(infStates,1);
logp = NaN(n,1);
yhat = NaN(n,1);
res  = NaN(n,1);

% Weed irregular trials out from responses and inputs
y = r.y(:,2);
y(r.irr) = [];

u_al = r.u(:,1);
state = u_al>0.5;

u_al(r.irr) = [];
state(r.irr) = [];

% Extract trajectories of interest from infStates
mu1hat = infStates(:,1,1);
mu2hat = infStates(:,2,1);
sa1hat = infStates(:,1,2);
sa2hat = infStates(:,2,2);
mu3hat = infStates(:,3,1);


% move variables from state (contingency space) to response space
% mu2hat_resp = mu2hat;
% mu2hat_resp(state ==0) = -mu2hat(state ==0);

% mu1hat_resp = mu1hat;
% mu1hat_resp(state ==0) = 1-mu1hat(state ==0);

% sahat1_resp = mu1hat_resp.*(1-mu1hat_resp);

% % Surprise
% % ~~~~~~~~
% m1hreg = mu1hat;
% m1hreg(r.irr) = [];
% poo = m1hreg.^u.*(1-m1hreg).^(1-u); % probability of observed outcome
% surp = -log2(poo);
% surp_shifted = [1; surp(1:(length(surp)-1))];

% Calculate predicted log-reaction time
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

mu1hat(r.irr) = [];
sa1hat(r.irr) = [];
sa2hat(r.irr) = [];
mu2hat(r.irr) = [];
logrt = be0 +be1.*mu2hat + be2.*sa1hat +be3.*sa2hat +be4.*mu3hat;

% Calculate log-probabilities for non-irregular trials
% Note: 8*atan(1) == 2*pi (this is used to guard against
% errors resulting from having used pi as a variable).
reg = ~ismember(1:n,r.irr);
logp(reg) = -1/2.*log(8*atan(1).*sa) -(y-logrt).^2./(2.*sa);
yhat(reg) = logrt;
res(reg) = y-logrt;

return;
