% Script to write model parameters to csv

clear; close all; clc

% load fits
m = importdata('STE_model7_fits_b0fixed.mat');
IDs = unique([m.ID])';

% load behavioural data
behavioural_data_dir = 'C:\Users\Tom\OneDrive - University of Cambridge\Cambridge\TOS_STE\DATA';
VAS = readtable([behavioural_data_dir, '\state_anxiety_VAS.csv']);
factor_scores = readtable([behavioural_data_dir, '\predictedFactorScores.csv']);
factor_scores = renamevars(factor_scores,'subjIDs','ID');
factor_scores = removevars(factor_scores, 'Group');

% concatenate behavioural data
T = outerjoin(VAS, factor_scores, 'Keys', 'ID', 'MergeKeys', true);

% preallocate space for model parameters
[T.S_om, T.S_b1, T.S_ze, T.T_om, T.T_b1, T.T_ze] = deal(nan(size(IDs)));

% add parameters to table
for i = 1:numel(m)
    idx = T.ID == m(i).ID;
    c = m(i).condition(1); %S or T
    if ~isempty(m(i).est)
        T(idx, [c, '_om']) = {m(i).est.p_prc.om(2)};
        T(idx, [c, '_b1']) = {m(i).est.p_obs.b1};
        T(idx, [c, '_ze']) = {m(i).est.p_obs.zeta};
    end
end

% calculate parameter differences
T.om_diff = T.S_om - T.T_om;
T.b1_diff = T.S_b1 - T.T_b1;
T.ze_diff = T.S_ze - T.T_ze;

% write csv
writetable(T, [behavioural_data_dir, '\HGFpsychometric_b0fixed_params.csv']);


