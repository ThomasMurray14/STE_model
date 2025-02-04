function pstruct = obs3_psychometric_namep(pvec)

pstruct = struct;

pstruct.b0 = pvec(1); % native space
pstruct.b1 = pvec(2); % native space
pstruct.alpha = pvec(3);

end