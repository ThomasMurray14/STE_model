% Script to preprocess demographic data from part 2 project

% read raw demographic data
A_raw = readtable('C:\Users\Tom\OneDrive - University of Cambridge\Cambridge\SpotTheEmotion_ThreatOfScream\DATA\2024\raw\Safe_first\data_exp_156254-v16\data_exp_156254-v16_questionnaire-o7yh');
B_raw = readtable('C:\Users\Tom\OneDrive - University of Cambridge\Cambridge\SpotTheEmotion_ThreatOfScream\DATA\2024\raw\Threat_first\data_exp_170385-v4\data_exp_170385-v4_questionnaire-o7yh');
C_raw = readtable('C:\Users\Tom\OneDrive - University of Cambridge\Cambridge\SpotTheEmotion_ThreatOfScream\DATA\2025\raw\safe_first\data_exp_156254-v17_questionnaire-o7yh');
D_raw = readtable('C:\Users\Tom\OneDrive - University of Cambridge\Cambridge\SpotTheEmotion_ThreatOfScream\DATA\2025\raw\threat_first\data_exp_170385-v5_questionnaire-o7yh');

% get list of IDs
A_IDs = unique(A_raw.ParticipantPrivateID(~isnan(A_raw.ParticipantPrivateID)));
B_IDs = unique(B_raw.ParticipantPrivateID(~isnan(B_raw.ParticipantPrivateID)));
C_IDs = unique(C_raw.ParticipantPrivateID(~isnan(C_raw.ParticipantPrivateID)));
D_IDs = unique(D_raw.ParticipantPrivateID(~isnan(D_raw.ParticipantPrivateID)));

% specify question mappings
qmap = {
    'resp_demo_1', 'Age';
    'resp_demo_2', 'Gender';
    'resp_demo_3', 'Ethnicity';
    'resp_demo_4', 'Country_of_residence';
    'resp_demo_5', 'Nationality';
    'response-3', 'Education_level';
    'resp_demo_6', 'Diagnoses';
    'resp_demo_7', 'Psychoactive_Medication';
    'resp_demo_8', 'Medication_name';
    'resp_demo_9', 'Medication_dosage_frequency';
    'resp_demo_10', 'Medication_duration'};

% Preallocate space for processed table
v_names = [{'subjIDs'; 'Group'}; qmap(:, 2)];
v_types = {...
    'doublenan';
    'string';
    'doublenan';
    'cellstr';
    'cellstr';
    'cellstr';
    'cellstr';
    'cellstr';
    'cellstr';
    'cellstr';
    'cellstr';
    'cellstr';
    'cellstr'
    };

%% A
A_processed = table(...
    'Size', [numel(A_IDs), numel(v_names)],...
    'VariableNames', v_names,...
    'VariableTypes', v_types);

% Loop through subjects
for iID = 1:numel(A_IDs)
    ID = A_IDs(iID);
    A_processed(iID, 'subjIDs') = {ID};
    A_processed(iID, 'Group') = {'A'};
    sub_raw = A_raw(A_raw.ParticipantPrivateID == ID, :);
    for iQ = 1:size(qmap, 1)
        q_resp = sub_raw.Response{...
            contains(sub_raw.QuestionKey, qmap{iQ, 1}) ...
            & ~contains(sub_raw.QuestionKey, 'quantised')};
        if iQ == 1; q_resp = str2double(q_resp); end
        A_processed(iID, qmap{iQ, 2}) = {q_resp};
    end
end

%% B
B_processed = table(...
    'Size', [numel(B_IDs), numel(v_names)],...
    'VariableNames', v_names,...
    'VariableTypes', v_types);

% Loop through subjects
for iID = 1:numel(B_IDs)
    ID = B_IDs(iID);
    B_processed(iID, 'subjIDs') = {ID};
    B_processed(iID, 'Group') = {'B'};
    sub_raw = B_raw(B_raw.ParticipantPrivateID == ID, :);
    for iQ = 1:size(qmap, 1)
        q_resp = sub_raw.Response{...
            contains(sub_raw.QuestionKey, qmap{iQ, 1}) ...
            & ~contains(sub_raw.QuestionKey, 'quantised')};
        if iQ == 1; q_resp = str2double(q_resp); end
        B_processed(iID, qmap{iQ, 2}) = {q_resp};
    end
end

%% C
C_processed = table(...
    'Size', [numel(C_IDs), numel(v_names)],...
    'VariableNames', v_names,...
    'VariableTypes', v_types);

% Loop through subjects
for iID = 1:numel(C_IDs)
    ID = C_IDs(iID);
    C_processed(iID, 'subjIDs') = {ID};
    C_processed(iID, 'Group') = {'C'};
    sub_raw = C_raw(C_raw.ParticipantPrivateID == ID, :);
    for iQ = 1:size(qmap, 1)
        q_resp = sub_raw.Response{...
            contains(sub_raw.QuestionKey, qmap{iQ, 1}) ...
            & ~contains(sub_raw.QuestionKey, 'quantised')};
        if iQ == 1; q_resp = str2double(q_resp); end
        C_processed(iID, qmap{iQ, 2}) = {q_resp};
    end
end

%% D
D_processed = table(...
    'Size', [numel(D_IDs), numel(v_names)],...
    'VariableNames', v_names,...
    'VariableTypes', v_types);

% Loop through subjects
for iID = 1:numel(D_IDs)
    ID = D_IDs(iID);
    D_processed(iID, 'subjIDs') = {ID};
    D_processed(iID, 'Group') = {'D'};
    sub_raw = D_raw(D_raw.ParticipantPrivateID == ID, :);
    for iQ = 1:size(qmap, 1)
        q_resp = sub_raw.Response{...
            contains(sub_raw.QuestionKey, qmap{iQ, 1}) ...
            & ~contains(sub_raw.QuestionKey, 'quantised')};
        if iQ == 1; q_resp = str2double(q_resp); end
        D_processed(iID, qmap{iQ, 2}) = {q_resp};
    end
end

% Concatenate
full_demographics  = [A_processed; B_processed; C_processed; D_processed];


%save
% writetable(full_demographics, 'demographics.csv');
