function [pvec, pstruct] = obs3_psychometric_transp(r, ptrans)

pvec    = NaN(1,length(ptrans));
pstruct = struct;

pvec(1)    = exp(ptrans(1)); % b0
pstruct.b0 = pvec(1);

pvec(2)    = exp(ptrans(2)); % b1
pstruct.b1 = pvec(2);

pvec(3)    = ptrans(3); % alpha
pstruct.alpha = pvec(3);


end

