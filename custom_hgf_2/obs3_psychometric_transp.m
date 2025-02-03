function [pvec, pstruct] = obs3_psychometric_transp(r, ptrans)

pvec    = NaN(1,length(ptrans));
pstruct = struct;

pvec(1)    = ptrans(1); % b0
pstruct.b0 = pvec(1);

pvec(2)    = ptrans(2); % b1
pstruct.b1 = pvec(2);

pvec(3)    = exp(ptrans(3)); % zeta
pstruct.zeta = pvec(3);


end

