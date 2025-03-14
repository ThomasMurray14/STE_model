% script to loop through parameter recoveries and generate figures

close all;

figdir = 'C:\Users\Tom\OneDrive - University of Cambridge\Cambridge\TOS_STE\Figures';

%%
m1 = importdata('STE_model1_recovery.mat');
pal_scat_ref_corr(m1.al_sim,        m1.al_est,      'alpha',    [figdir, '\model1_recovery_alpha.bmp']);
pal_scat_ref_corr(m1.om2_sim,       m1.om2_est,     'omega',    [figdir, '\model1_recovery_omega.bmp']);
pal_scat_ref_corr(m1.rho2_sim,      m1.rho2_est,    'rho',      [figdir, '\model1_recovery_rho.bmp']);
pal_scat_ref_corr(log(m1.al_sim),   log(m1.al_est), 'log(alpha)', [figdir, '\model1_recovery_logalpha.bmp']);

%%
m2 = importdata('STE_model2_recovery.mat');
pal_scat_ref_corr(m2.b0_sim,        m2.b0_est,      'beta0',    [figdir, '\model2_recovery_beta0.bmp']);
pal_scat_ref_corr(m2.b1_sim,        m2.b1_est,      'beta1',    [figdir, '\model2_recovery_beta1.bmp']);
pal_scat_ref_corr(m2.om2_sim,       m2.om2_est,     'omega2',   [figdir, '\model2_recovery_omega2.bmp']);
pal_scat_ref_corr(m2.om3_sim,       m2.om3_est,     'omega3',   [figdir, '\model2_recovery_omega3.bmp']);
pal_scat_ref_corr(m2.zeta_sim,      m2.zeta_est,    'zeta',     [figdir, '\model2_recovery_zeta.bmp']);
pal_scat_ref_corr(log(m2.zeta_sim),      log(m2.zeta_est),    'log(zeta)',     [figdir, '\model2_recovery_logzeta.bmp']);

%%
m3 = importdata('STE_model3_recovery.mat');
pal_scat_ref_corr(m3.b0_sim,        m3.b0_est,      'beta0',    [figdir, '\model3_recovery_beta0.bmp']);
pal_scat_ref_corr(m3.b1_sim,        m3.b1_est,      'beta1',    [figdir, '\model3_recovery_beta1.bmp']);
pal_scat_ref_corr(m3.om2_sim,       m3.om2_est,     'omega2',   [figdir, '\model3_recovery_omega2.bmp']);
pal_scat_ref_corr(m3.om3_sim,       m3.om3_est,     'omega3',   [figdir, '\model3_recovery_omega3.bmp']);
pal_scat_ref_corr(m3.alpha_sim,     m3.alpha_est,   'alpha',    [figdir, '\model3_recovery_alpha.bmp']);


%%
m4 = importdata('STE_model4_recovery.mat');
pal_scat_ref_corr(m4.b0_sim,        m4.b0_est,      'beta0',    [figdir, '\model4_recovery_beta0.bmp']);
pal_scat_ref_corr(m4.b1_sim,        m4.b1_est,      'beta1',    [figdir, '\model4_recovery_beta1.bmp']);
pal_scat_ref_corr(m4.lambda_sim,    m4.lambda_est,  'lambda',   [figdir, '\model4_recovery_lambda.bmp']);
pal_scat_ref_corr(m4.omega_sim,     m4.omega_est,   'omega',    [figdir, '\model4_recovery_omega.bmp']);
pal_scat_ref_corr(m4.zeta_sim,      m4.zeta_est,    'zeta',     [figdir, '\model4_recovery_zeta.bmp']);


pal_scat_ref_corr(log(m4.omega_sim),     log(m4.omega_est),   'log(omega)',    [figdir, '\model4_recovery_logomega.bmp']);
pal_scat_ref_corr(log(m4.zeta_sim),      log(m4.zeta_est),    'log(zeta)',     [figdir, '\model4_recovery_logzeta.bmp']);


%%
m5 = importdata('STE_model5_recovery.mat');
pal_scat_ref_corr(m5.b0_sim,        m5.b0_est,      'beta0',    [figdir, '\model5_recovery_beta0.bmp']);
pal_scat_ref_corr(m5.b1_sim,        m5.b1_est,      'beta1',    [figdir, '\model5_recovery_beta1.bmp']);
pal_scat_ref_corr(m5.alpha_sim,     m5.alpha_est,   'alpha',    [figdir, '\model5_recovery_alpha.bmp']);
pal_scat_ref_corr(m5.lambda_sim,    m5.lambda_est,  'lambda',   [figdir, '\model5_recovery_lambda.bmp']);
pal_scat_ref_corr(m5.omega_sim,     m5.omega_est,   'omega',    [figdir, '\model5_recovery_omega.bmp']);

pal_scat_ref_corr(log(m5.omega_sim),     log(m5.omega_est),   'log(omega)',    [figdir, '\model5_recovery_logomega.bmp']);




%% 
m6 = importdata('STE_model6_recovery_rho_scaled.mat');
pal_scat_ref_corr(m6.al_sim,        m6.al_est,      'alpha',    [figdir, '\model6rhoscaled_recovery_alpha.bmp']);
pal_scat_ref_corr(m6.om2_sim,       m6.om2_est,     'omega',    [figdir, '\model6rhoscaled_recovery_omega.bmp']);
pal_scat_ref_corr(m6.rho2_sim,      m6.rho2_est,    'rho',      [figdir, '\model6rhoscaled_recovery_rho.bmp']);
pal_scat_ref_corr(m6.ze_sim,        m6.ze_est,      'zeta',     [figdir, '\model6rhoscaled_recovery_zeta.bmp']);
pal_scat_ref_corr(log(m6.al_sim),   log(m6.al_est), 'log(alpha)', [figdir, '\model6rhoscaled_recovery_logalpha.bmp']);
pal_scat_ref_corr(log(m6.ze_sim),   log(m6.ze_est), 'log(zeta)',  [figdir, '\model6rhoscaled_recovery_logzeta.bmp']);


%% 
m6 = importdata('STE_model6_recovery.mat');
pal_scat_ref_corr(m6.al_sim,        m6.al_est,      'alpha',    [figdir, '\model6_recovery_alpha.bmp']);
pal_scat_ref_corr(m6.om2_sim,       m6.om2_est,     'omega',    [figdir, '\model6_recovery_omega.bmp']);
pal_scat_ref_corr(m6.rho2_sim,      m6.rho2_est,    'rho',      [figdir, '\model6_recovery_rho.bmp']);
pal_scat_ref_corr(m6.ze_sim,        m6.ze_est,      'zeta',     [figdir, '\model6_recovery_zeta.bmp']);
pal_scat_ref_corr(log(m6.al_sim),   log(m6.al_est), 'log(alpha)', [figdir, '\model6_recovery_logalpha.bmp']);
pal_scat_ref_corr(log(m6.ze_sim),   log(m6.ze_est), 'log(zeta)',  [figdir, '\model6_recovery_logzeta.bmp']);









