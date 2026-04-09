% Script to organise Hopkins reduced set of questions (from Gorilla) into
% csv, to input into python script.

raw_files = {...
    'C:\Users\Tom\OneDrive - University of Cambridge\Cambridge\SpotTheEmotion_ThreatOfScream\DATA\2024\raw\Safe_first\data_exp_156254-v16\data_exp_156254-v16_questionnaire-fo9f.csv';
    'C:\Users\Tom\OneDrive - University of Cambridge\Cambridge\SpotTheEmotion_ThreatOfScream\DATA\2024\raw\Threat_first\data_exp_170385-v4\data_exp_170385-v4_questionnaire-fo9f.csv';
    'C:\Users\Tom\OneDrive - University of Cambridge\Cambridge\SpotTheEmotion_ThreatOfScream\DATA\2025\raw\safe_first\data_exp_156254-v17_questionnaire-fo9f.csv';
    'C:\Users\Tom\OneDrive - University of Cambridge\Cambridge\SpotTheEmotion_ThreatOfScream\DATA\2025\raw\threat_first\data_exp_170385-v5_questionnaire-fo9f.csv'};
processed_tables = cell(4,1);

group_labels = {'A','B','C','D'};

for iRaw = 1:4
    % read raw table
    raw = readtable(raw_files{iRaw});
    
    % get list of IDs
    IDs = unique(raw.ParticipantPrivateID(~isnan(raw.ParticipantPrivateID)));
    
    % preallocate table for reduced set
    items = unique(raw.ObjectName(~cellfun(@isempty, raw.ObjectName)));
    reduced_v_names = [{'subjIDs'; 'Group'}; items];
    reduced_v_types = repmat({'doublenan'}, size(reduced_v_names));
    reduced_v_types{2} = 'string';
    reduced_set = table(...
        'Size', [numel(IDs), numel(reduced_v_names)],...
        'VariableNames', reduced_v_names,...
        'VariableTypes', reduced_v_types);
    
    % loop through IDs
    for iID = 1:numel(IDs)
        ID = IDs(iID);
        reduced_set(iID, 'subjIDs') = {ID};
        reduced_set(iID, 'Group') = {group_labels{iRaw}};
        sub_raw = raw(raw.ParticipantPrivateID == ID, :);
        sub_q = sub_raw(strcmp(sub_raw.Key, 'quantised'), :);
        for iItem = 1:numel(items)
            resp = sub_q.Response(strcmp(sub_q.ObjectName, items{iItem}));
            reduced_set(iID, items{iItem}) = {resp};
        end
    end
    
    % remove 'Rating Scale' (example scale in instructions)
    reduced_set = removevars(reduced_set, 'Rating Scale');

    % Re scale EAT (called EA in Gorilla output)
    reduced_set(:, {'EA1','EA2'}) = reduced_set(:, {'EA1','EA2'}) - 3;
    reduced_set(reduced_set{:, {'EA1'}}<0, {'EA1'}) = {0};
    reduced_set(reduced_set{:, {'EA2'}}<0, {'EA2'}) = {0};
    
    % rename incorrect OCI-R item
    reduced_set = renamevars(reduced_set, 'OC4', 'OCI4');
    reduced_v_names = reduced_set.Properties.VariableNames;
    
    % Subtract 1 from LSAS (i.e. so scores start at 0)
    LSAS_items = reduced_v_names(contains(reduced_v_names, 'LSAS'));
    reduced_set(:, LSAS_items) = reduced_set(:, LSAS_items) -1;
    
    % Subtract 1 from OCI-R
    OCI_items = reduced_v_names(contains(reduced_v_names, 'OCI'));
    reduced_set(:, OCI_items) = reduced_set(:, OCI_items) -1;
    
    % find mean of LSAS and remove raw
    LSAS_items = [2,5,6,7,8,9,10,11,12,14,15,16,18,19,20,21,22,23,24];
    for i = LSAS_items
        raw_LSAS = {['LSAS', num2str(i) 'a'], ['LSAS', num2str(i), 'b']};
        reduced_set(:, ['LSAS', num2str(i)]) = mean(reduced_set(:, raw_LSAS), 2);
        reduced_set = removevars(reduced_set, raw_LSAS);
    end
    
    % Reverse score SDS items (ZDS in raw form)
    SDS_to_reverse = {'ZDS11', 'ZDS12', 'ZDS14', 'ZDS16', 'ZDS17', 'ZDS18', 'ZDS20'};
    reduced_set(:, SDS_to_reverse) = 5 - reduced_set(:, SDS_to_reverse);
    
    % Reverse score STAI items
    STAI_to_reverse = {'STAI1', 'STAI3', 'STAI7', 'STAI10', 'STAI13', 'STAI16', 'STAI19'};
    reduced_set(:, STAI_to_reverse) = 5 - reduced_set(:, STAI_to_reverse);
    
    % Reverse score BIS items
    BIS_to_reverse = {'BIS9', 'BIS13', 'BIS20'};
    reduced_set(:, BIS_to_reverse) = 5 - reduced_set(:, BIS_to_reverse);
    
    % Reverse score AES (all)
    AES_items = reduced_v_names(contains(reduced_v_names, 'AES'));
    reduced_set(:, AES_items) = 5 - reduced_set(:, AES_items);
    
    % Get lists of old and new item names
    old_scales = {'ZDS', 'STAI', 'OCI', 'LSAS', 'BIS', 'SCZ', 'AUDIT', 'EA', 'AES'}; %(as they come off gorilla)
    new_scales = {'SDS', 'STAI', 'OCI', 'LSAS', 'BIS', 'SCZ', 'AUDIT', 'EAT', 'AES'};
    scale_N = [20, 20, 18, 24, 30, 43, 10, 26, 18];
    old_v_names = {'subjIDs'; 'Group'};
    new_v_names = {'subjIDs'; 'Group'};
    for iScale = 1:numel(old_scales)
        for iN = 1:scale_N(iScale)
            old_v_names = [old_v_names; {[old_scales{iScale}, num2str(iN)]}];
            new_v_names = [new_v_names; {[new_scales{iScale}, '_', num2str(iN)]}];
        end
    end
    
    % Create full table (using old scale names) - python script takes table
    % containing all questions, with 0's in non-used questions
    full_v_types = repmat({'doublenan'}, size(old_v_names));
    full_v_types{2} = 'string';
    group_table = table(...
        'Size', [numel(IDs), numel(old_v_names)],...
        'VariableNames', old_v_names,...
        'VariableTypes', full_v_types);
    
    % populate all cells with 0
    group_table{:,:} = 0;
    
    % populate with responses to reduced set
    group_table(:, reduced_set.Properties.VariableNames) = reduced_set(:, reduced_set.Properties.VariableNames);
    
    % change old item names to new item names
    group_table = renamevars(group_table, old_v_names, new_v_names);
    
    % Add processed table to cell
    processed_tables{iRaw} = group_table;
end

% concatenate two tables (group is not needed here)
full_table = [processed_tables{1}; processed_tables{2}; processed_tables{3}; processed_tables{4}];

% save table
writetable(full_table, 'transDiagQs_for_python.csv');

