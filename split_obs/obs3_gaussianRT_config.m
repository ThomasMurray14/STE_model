function [c] = obs3_gaussianRT_config()
% [c] = obs3_comb_obs_config()
%
% Contains the prior configurations of the combined response model M1.
% (Designed to be compatible with the HGF Toolbox as part of TAPAS).
%
% INPUT
%   argin         type  
%
%   OPTIONAL:
%
% OUTPUT    
%   c             struct       Struct containing obs model prior configs
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


%% Config structure
c = struct;

% Is the decision based on predictions or posteriors? Comment as appropriate.
% c.predorpost = 1; % Predictions
c.predorpost = 2; % Posteriors

% Model name
c.model = 'obs3_gaussianRT';


%% Sufficient statistics of Gaussian parameter priors


%-----------------------------------
% Model for Reaction Time fit
% Gaussian function:
% gamma+(lambda*exp(-((x-mu)^2)/(2*(sigma^2))))

% mu (logit space)
c.logitmumu=tapas_logit(.5, 1);
c.logitmusa=4;

% sigma0 (log space)
c.logsigma0mu=log(0.2);
c.logsigma0sa=4;

% gamma1 (log space)
c.loggamma1mu=log(0.2);
c.loggamma1sa=4;

% lambda1 (log space)
c.loglambda1mu=log(0.2);
c.loglambda1sa=4;

% sigma1 (log space) - noise term
c.logsigma1mu = log(2);
c.logsigma1sa = 4;


%% Gather prior settings in vectors
c.priormus = [
    c.logitmumu,...
    c.logsigma0mu,...
    c.loggamma1mu,...
    c.loglambda1mu,...
    c.logsigma1mu,...
    ];


c.priorsas = [
    c.logitmusa,...
    c.logsigma0sa,...
    c.loggamma1sa,...
    c.loglambda1sa,...
    c.logsigma1sa,...
    ];

% Model filehandle
c.obs_fun = @obs3_gaussianRT;

% Handle to function that transforms perceptual parameters to their native
% space from the space they are estimated in
c.transp_obs_fun = @obs3_gaussianRT_transp;

return;
