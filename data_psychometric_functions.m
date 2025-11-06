% script to visualise the psychometric functions in the actual data



STE_dir = dir('STE_data\*.csv');
N_files = numel(STE_dir);

[psad_all, rt_all] = deal(zeros(N_files, 6));

for i = 1:N_files

    % get ID, group, condition
    f_name = fullfile(STE_dir(i).folder, STE_dir(i).name);
    
    % load data
    sub_data = readtable(f_name);

    % get p(sad)
    sub_psad = arrayfun(@(x) mean(sub_data.Response_idx(sub_data.Outcome_p_sad==x), 'omitnan'), 0:20:100);

    % get RT
    sub_rt = log(arrayfun(@(x) mean(sub_data.Response_RT(sub_data.Outcome_p_sad==x), 'omitnan'), 0:20:100));

    % add to full matrix
    psad_all(i,:) = sub_psad;
    rt_all(i,:) = sub_rt;

end

% remove missed
valid = ~any(isnan(psad_all), 2);
psad_all = psad_all(valid, :);
rt_all = rt_all(valid, :);



figure;
subplot(2,1,1); 
hold on;
mean_psad = mean(psad_all, 1);
% plot with jitter
for i=1:size(psad_all)
    scatter(...
        arrayfun(@(x) x+((rand-0.5)/3), 1:6), ...
        psad_all(i,:), 5, [.7,.7,.7])
end
% plot mean
plot(1:6, mean_psad, '-o',...
    'color', [0,0.4470,0.7410],...
    'linewidth', 1.5, 'MarkerSize',8,...
    'MarkerEdgeColor',[0,0.4470,0.7410],...
    'MarkerFaceColor',[1,1,1])

set(gca, 'xlim', [1,6])


subplot(2,1,2); 
hold on;
% plot with jitter
for i=1:size(psad_all)
    scatter(arrayfun(@(x) x+((rand-0.5)/3), 1:6), rt_all(i,:), 5, [.7,.7,.7])
end
% plot mean
plot(1:6, mean(rt_all, 1), '-o',...
    'color', [0,0.4470,0.7410],...
    'linewidth', 1.5, 'MarkerSize',8,...
    'MarkerEdgeColor',[0,0.4470,0.7410],...
    'MarkerFaceColor',[1,1,1])

set(gca,'xlim', [1,6], 'xtick', 1:6)










