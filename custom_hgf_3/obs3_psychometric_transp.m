function [pvec, pstruct] = obs3_psychometric_transp(r, ptrans)

pvec    = NaN(1,length(ptrans));
pstruct = struct;

pvec(1)    = ptrans(1); % b0
pstruct.b0 = pvec(1);

pvec(2)    = ptrans(2); % b1
pstruct.b1 = pvec(2);

pvec(3)    = ptrans(3); % b2
pstruct.b2 = pvec(3);

pvec(4)    = ptrans(4); % b3
pstruct.b3 = pvec(4);

pvec(5)    = ptrans(5); % b4
pstruct.b4 = pvec(5);

end