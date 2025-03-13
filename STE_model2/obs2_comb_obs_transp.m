function [pvec, pstruct] = obs2_comb_obs_transp(r, ptrans)
% [pvect, pstruct] = comb_obs_transp(r, ptrans)
%
% Transforms parameter values from estimation into native space.
% (Designed to be compatible with the HGF Toolbox as part of TAPAS).
%
% INPUT
%   r            struct      Struct obtained from tapas_fitModel.m fct
%   ptrans       vector      1xP vector containing model params (est space)
%
%   OPTIONAL:
%
% OUTPUT    
%   pvect        vector      1xP vector containing model params (nat space)
%   pstruct      struct      Empty struct
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



% empty array
pstruct = struct();

% vector with parameter values transformed back into native space
pvec = ptrans;

% Psychometric
pvec(1)     = ptrans(1); % b0 (baseline PSE)
pstruct.b0  = pvec(1);

pvec(2)     = ptrans(2); % b1 (volatility-on-PSE)
pstruct.b0  = pvec(2);

pvec(3)     = exp(ptrans(3)); % zeta (slope)
pstruct.zeta = pvec(3);

% logRT
pvec(4)     = ptrans(4);         % beta0
pstruct.beta0 = pvec(4);

pvec(5)     = ptrans(5);         % beta1
pstruct.beta1 = pvec(5);

pvec(6)     = ptrans(6);         % beta2
pstruct.beta2 = pvec(6);

pvec(7)     = ptrans(7);         % beta3
pstruct.beta3 = pvec(7);

pvec(8)     = ptrans(8);         % beta4
pstruct.beta4 = pvec(8);

pvec(9)     = exp(ptrans(9));    % sa (logRT model noise parameter)
pstruct.sa  = pvec(9);


end
