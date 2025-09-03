function pstruct = prc2_ehgf_binary_pu_tbt_namep(pvec)
% --------------------------------------------------------------------------------------------------
% Copyright (C) 2012-2015 Christoph Mathys, TNU, UZH & ETHZ
%
% This file is part of the HGF toolbox, which is released under the terms of the GNU General Public
% Licence (GPL), version 3. You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version). For further details, see the file
% COPYING or <http://www.gnu.org/licenses/>.


pstruct = struct;

% l = (length(pvec)-2)/5;
% 
% if l ~= floor(l)
%     error('tapas:hgf:UndetNumLevels', 'Cannot determine number of levels');
% end

pstruct.mu_0      = pvec(1:3);
pstruct.sa_0      = pvec(4:6);
pstruct.rho       = pvec(7);
pstruct.ka        = pvec(8:9);
pstruct.om        = pvec(10:12);
pstruct.al        = pvec(13);
pstruct.eta0      = pvec(14);
pstruct.eta1      = pvec(15);




return;
