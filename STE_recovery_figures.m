% script to loop through parameter recoveries and generate figures

figdir = 'C:\Users\Tom\OneDrive - University of Cambridge\Cambridge\TOS_STE\Figures';

%%
m1 = importdata('STE_model1_recovery.mat');
m1_name_map.al = 'alpha';
m1_name_map.om2 = 'omega2';
m1_name_map.rho2 = 'rho2';
p_names = fieldnames(m1_name_map);
for i=1:numel(p_names)
    est = m1.([p_names{i}, '_est']);
    sim = m1.([p_names{i}, '_sim']);
    name = m1_name_map.(p_names{i});
    
    valid = ~isnan(est);
    f = pal_scat_ref_corr(sim(valid), est(valid), name);
    saveas(f, [figdir, '\model1_recovery_', name, '.bmp']);
    close all;
end

%%
m2 = importdata('STE_model2_recovery.mat');
m2_name_map.b0 = 'beta0';
m2_name_map.b1 = 'beta1';
m2_name_map.om2 = 'omega2';
m2_name_map.om3 = 'omega3';
m2_name_map.zeta = 'zeta';
p_names = fieldnames(m2_name_map);
for i=1:numel(p_names)
    est = m2.([p_names{i}, '_est']);
    sim = m2.([p_names{i}, '_sim']);
    name = m2_name_map.(p_names{i});
    
    valid = ~isnan(est);
    f = pal_scat_ref_corr(sim(valid), est(valid), name);
    saveas(f, [figdir, '\model2_recovery_', name, '.bmp']);
    close all;
end

%%
m3 = importdata('STE_model3_recovery.mat');
m3_name_map.b0 = 'beta0';
m3_name_map.b1 = 'beta1';
m3_name_map.om2 = 'omega2';
m3_name_map.om3 = 'omega3';
m3_name_map.alpha = 'alpha';
p_names = fieldnames(m3_name_map);
for i=1:numel(p_names)
    est = m3.([p_names{i}, '_est']);
    sim = m3.([p_names{i}, '_sim']);
    name = m3_name_map.(p_names{i});
    
    valid = ~isnan(est);
    f = pal_scat_ref_corr(sim(valid), est(valid), name);
    saveas(f, [figdir, '\model3_recovery_', name, '.bmp']);
    close all;
end

%%
m4 = importdata('STE_model4_recovery.mat');
m4_name_map.b0 = 'beta0';
m4_name_map.b1 = 'beta1';
m4_name_map.lambda = 'lambda';
m4_name_map.omega = 'omega';
m4_name_map.zeta = 'zeta';
p_names = fieldnames(m4_name_map);
for i=1:numel(p_names)
    est = m4.([p_names{i}, '_est']);
    sim = m4.([p_names{i}, '_sim']);
    name = m4_name_map.(p_names{i});
    
    valid = ~isnan(est);
    f = pal_scat_ref_corr(sim(valid), est(valid), name);
    saveas(f, [figdir, '\model4_recovery_', name, '.bmp']);
    close all;
end

%%
m5 = importdata('STE_model5_recovery.mat');
m5_name_map.b0 = 'beta0';
m5_name_map.b1 = 'beta1';
m5_name_map.alpha = 'alpha';
m5_name_map.lambda = 'lambda';
m5_name_map.omega = 'omega';
p_names = fieldnames(m5_name_map);
for i=1:numel(p_names)
    est = m5.([p_names{i}, '_est']);
    sim = m5.([p_names{i}, '_sim']);
    name = m5_name_map.(p_names{i});
    
    valid = ~isnan(est);
    f = pal_scat_ref_corr(sim(valid), est(valid), name);
    saveas(f, [figdir, '\model5_recovery_', name, '.bmp']);
    close all;
end