function visualise_psychometric(u, sub_data, prc_name, prc_params, obs_name, obs_params, N)
% Function to simulate responses N times and visualise the p(sad) and RT
% distributions

all_y = nan(size(u,1), 2, N);
for i = 1:N
    sim = tapas_simModel(u,...
        prc_name,...
        prc_params,...
        obs_name,...
        obs_params);
    all_y(:,:,i) = sim.y;
end

y_resp = squeeze(all_y(:,1,:)); % uxN
y_RT = squeeze(all_y(:,2,:));

figure;

% visualise psychometric
subplot(2,1,1); hold on;
all_resp_dists = zeros(N, 6);
for i = 1:N
    sim_sad = (sub_data.Cue_idx == 1 & y_resp(:,i) == 1) + (sub_data.Cue_idx == 0 & y_resp(:,i) == 0);
    sim_resp_dist = arrayfun(@(x) mean(sim_sad(sub_data.Outcome_p_sad==x, 1)), 0:20:100);
    all_resp_dists(i, :) = sim_resp_dist;
    plot(0:20:100, sim_resp_dist, 'linewidth', 1, 'color', [.7,.7,.7]);
end
mean_resp = mean(all_resp_dists, 1);
plot(0:20:100, mean_resp, 'linewidth', 3, 'Color', [0    0.4470    0.7410]);
set(gca, 'Ylim', [0,1], 'Xtick', 0:20:100)
ylabel('p(Sad)')
title('simulated response distribution');


% visualise RT
subplot(2,1,2); hold on;
all_RT_dists = zeros(N, 6);
for i = 1:N
    sim_RT = y_RT(:, i);%(sub_data.Cue_idx == 1 & y_RT(:,i) == 1) + (sub_data.Cue_idx == 0 & y_RT(:,i) == 0);
    sim_RT_dist = arrayfun(@(x) mean(sim_RT(sub_data.Outcome_p_sad==x, 1)), 0:20:100);
    all_RT_dists(i, :) = sim_RT_dist;
    plot(0:20:100, sim_RT_dist, 'linewidth', 1, 'color', [.7,.7,.7]);
end
mean_RT = mean(all_RT_dists, 1);
plot(0:20:100, mean_RT, 'linewidth', 3, 'Color', [0    0.4470    0.7410]);
set(gca, 'Xtick', 0:20:100)
xlabel('%Sad');
ylabel('logRT');
title('simulated RT distribution');

end
