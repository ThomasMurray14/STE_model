function y = obs_RW_logrt_linear_binary_sim(r, infStates, p)
% Simulates logRTs with Gaussian noise
%
% --------------------------------------------------------------------------------------------------
% Copyright (C) 2016 Christoph Mathys, TNU, UZH & ETHZ
%
% This file is part of the HGF toolbox, which is released under the terms of the GNU General Public
% Licence (GPL), version 3. You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version). For further details, see the file
% COPYING or <http://www.gnu.org/licenses/>.

% Get parameters
be0  = p(1);
be1  = p(2);
be2  = p(3);
sa   = p(4);

% Number of trials
n = size(infStates,1);

% inputs and stim noise
u_al = r.u(:,1);
u = u_al>0.5;
stim_noise = 0.5-abs(u_al-.5); 



% Extract trajectories of interest from infStates
da = r.traj.da; % prediction error
v = r.traj.v; % posterior
vhat = r.traj.vhat; % prediction




% Calculate predicted log-reaction time
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
logrt = be0 +be1.*abs(da) + be2.*stim_noise;

% Initialize random number generator
if isnan(r.c_sim.seed)
    rng('shuffle');
else
    rng(r.c_sim.seed);
end

% Simulate
y = logrt+sqrt(sa)*randn(n, 1);

end
