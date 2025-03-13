function y = obs_logrt_linear_binary_sim(r, infStates, p)
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
be3  = p(4);
% be4  = p(5);
sa   = p(5);

% Number of trials
n = size(infStates,1);

% Inputs
% y = r.y(:,2);
% y(r.irr) = [];

u_al = r.u(:,1);
state = u_al>0.5;

% u_al(r.irr) = [];
% state(r.irr) = [];

% Extract trajectories of interest from infStates
mu1hat = infStates(:,1,1);
mu2hat = infStates(:,2,1);
sa1hat = infStates(:,1,2);
sa2hat = infStates(:,2,2);
mu3hat = infStates(:,3,1);

% move variables from state (contingency space) to response space
mu2hat_resp = mu2hat;
mu2hat_resp(state ==0) = -mu2hat(state ==0);

mu1hat_resp = mu1hat;
mu1hat_resp(state ==0) = 1-mu1hat(state ==0);

sahat1_resp = mu1hat_resp.*(1-mu1hat_resp);


% Surprise
% ~~~~~~~~
% poo = mu1hat.^u.*(1-mu1hat).^(1-u); % probability of observed outcome
% surp = -log2(poo);
% surp_shifted = [1; surp(1:(length(surp)-1))];

% Calculate predicted log-reaction time
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 
% mu1hat_resp(r.irr) = [];
% sahat1_resp(r.irr) = [];
% sa2hat(r.irr) = [];
% mu3hat(r.irr) = [];
logrt = be0 +be1.*mu2hat + be2.*sa1hat +be3.*sa2hat;% +be4.*mu3hat;

% 
% logrt = be0 + be1.*surp_shifted + be2.*sa1hat + be3.*sa2hat + ...
%     be4.*exp(mu3hat);

% Initialize random number generator
if isnan(r.c_sim.seed)
    rng('shuffle');
else
    rng(r.c_sim.seed);
end

% Simulate
y = logrt+sqrt(sa)*randn(n, 1);

end
