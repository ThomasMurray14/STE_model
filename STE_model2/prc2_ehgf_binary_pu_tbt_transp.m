function [pvec, pstruct] = prc2_ehgf_binary_pu_tbt_transp(r, ptrans)
% --------------------------------------------------------------------------------------------------
% Copyright (C) 2012-2015 Christoph Mathys, TNU, UZH & ETHZ
%
% This file is part of the HGF toolbox, which is released under the terms of the GNU General Public
% Licence (GPL), version 3. You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version). For further details, see the file
% COPYING or <http://www.gnu.org/licenses/>.


pvec    = NaN(1,length(ptrans));
pstruct = struct;

% l = r.c_prc.n_levels;

pvec(1:3)         = ptrans(1:3);                           % mu_0
pstruct.mu_0      = pvec(1:3);

pvec(4:6)         = exp(ptrans(4:6));                  % sa_0
pstruct.sa_0      = pvec(4:6);

pvec(7)           = tapas_sgm(ptrans(7), 1);                     % rho
pstruct.rho       = pvec(7);

pvec(8:9)       = exp(ptrans(8:9));              % ka
pstruct.ka        = pvec(8:9);

pvec(10:12)   = ptrans(10:12);                     % om
pstruct.om        = pvec(10:12);

pvec(13)         = exp(ptrans(13));                      % al
pstruct.al        = pvec(13);

pvec(14)       = ptrans(14);                         % eta0
pstruct.eta0      = pvec(14);

pvec(15)       = ptrans(15);                         % eta1
pstruct.eta1      = pvec(15);

return;
