% Script to get csv for each subject for STE tasks (safe and threat)

data_dir = 'C:\Users\Tom\OneDrive - University of Cambridge\Cambridge\SpotTheEmotion_ThreatOfScream\DATA\';

raw_files(1).group = 'A';
raw_files(1).safe = '2024\raw\Safe_first\data_exp_156254-v16\data_exp_156254-v16_task-na2a.csv';
raw_files(1).threat = '2024\raw\Safe_first\data_exp_156254-v16\data_exp_156254-v16_task-p4ov.csv';

raw_files(2).group = 'B';
raw_files(2).safe = '2024\raw\Threat_first\data_exp_170385-v4\data_exp_170385-v4_task-na2a.csv';
raw_files(2).threat = '2024\raw\Threat_first\data_exp_170385-v4\data_exp_170385-v4_task-p4ov.csv';

raw_files(3).group = 'C';
raw_files(3).safe = '2025\raw\safe_first\data_exp_156254-v17_task-na2a.csv';
raw_files(3).threat = '2025\raw\safe_first\data_exp_156254-v17_task-p4ov.csv';

raw_files(4).group = 'D';
raw_files(4).safe = '2025\raw\threat_first\data_exp_170385-v5_task-na2a.csv';
raw_files(4).threat = '2025\raw\threat_first\data_exp_170385-v5_task-p4ov.csv';

% set variables in new tables
v_names = {...
    'ID',...
    'Group',...
    'Safe_Check_1',...
    'Safe_Check_2',...
    'Safe_Check_3',...
    'Threat_Check_1',...
    'Threat_Check_2',...
    'Threat_Check_3'};
v_types = {...
    'doublenan',...
    'cellstr',...
    'doublenan',...
    'doublenan',...
    'doublenan',...
    'doublenan',...
    'doublenan',...
    'doublenan',...
    };

processed_tables = cell(2,1);

for iRaw = 1:4

    data.Safe = readtable([data_dir, raw_files(iRaw).safe]);
    data.Threat = readtable([data_dir, raw_files(iRaw).threat]);
    
    % get list of IDs
    IDs = unique(data.Safe.ParticipantPrivateID(~isnan(data.Safe.ParticipantPrivateID)));
    
    % preallocate space in table
    group_table = table(...
        'Size', [numel(IDs), numel(v_names)],...
        'VariableNames', v_names,...
        'VariableTypes', v_types);
    
    conds = {'Safe', 'Threat'};
    for iC = 1:2
        d = data.(conds{iC});
        for iID = 1:numel(IDs)
            
            % current ID
            ID = IDs(iID);
            
            % add ID to table
            group_table(iID, 'ID') = {ID};
            
            % Add group
            group_table(iID, 'Group') = {raw_files(iRaw).group};
            
            % this subject's data
            sub_raw = d(d.ParticipantPrivateID == ID, :);
            
            % remove anything that isn't the attention check
            sub_raw = sub_raw(strcmp(sub_raw.display, 'Attention_check'), :);
    
            % remove anything that isn't their responses
            sub_raw = sub_raw(strcmp(sub_raw.ZoneType, 'response_text_entry'), :);
            
            for iA = 1:3
                try 
                    resp = str2double(sub_raw.Response{iA});
                catch
                    resp = 0;
                end
                correct = resp == sub_raw.Attention_number(iA);
                group_table(iID, [conds{iC}, '_Check_', num2str(iA)]) = {correct};
            end
        end
    end
    
    processed_tables{iRaw} = group_table;

end

full_table = [processed_tables{1}; processed_tables{2}; processed_tables{3}; processed_tables{4}];

% writetable(full_table, 'attention_check.csv');



