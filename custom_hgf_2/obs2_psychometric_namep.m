function pstruct = obs2_psychometric_namep(pvec)

pstruct = struct;

pstruct.b0 = pvec(1);
pstruct.b1 = pvec(2);
pstruct.beta = pvec(3); % native space

end