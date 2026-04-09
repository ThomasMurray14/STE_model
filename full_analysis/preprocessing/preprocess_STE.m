 % Script to get csv for each subject for STE tasks (safe and threat)

% Root directory
data_dir = 'C:\Users\Tom\OneDrive - University of Cambridge\Cambridge\SpotTheEmotion_ThreatOfScream\DATA\2025\raw\';

% set variables in new table per subject
v_names = {...
    'Trial_N',...
    'Block_N',...
    'Contingency',...
    'Expectedness',...
    'Cue_idx',...
    'Cue_name',...
    'Outcome_idx',...
    'Outcome_name',...
    'Outcome_p_sad',...
    'Response_idx',...
    'Response_name',...
    'Response_RT',...
    'Confidence_idx',...
    'Confidence_name',...
    'Confidence_RT'};
v_types = {...
    'doublenan',...
    'doublenan',...
    'doublenan',...
    'cellstr',...
    'doublenan',...
    'cellstr',...
    'doublenan',...
    'cellstr',...
    'doublenan',...
    'doublenan',...
    'cellstr',...
    'doublenan',...
    'doublenan',...
    'cellstr',...
    'doublenan'};

% Name raw data files
raw_files(1).name = 'safe_first\data_exp_156254-v17_task-na2a.csv';
raw_files(1).condition = 'Safe';
raw_files(1).group = 'C'; % C = safe first

raw_files(2).name = 'safe_first\data_exp_156254-v17_task-p4ov.csv';
raw_files(2).condition = 'Threat';
raw_files(2).group = 'C';

raw_files(3).name = 'threat_first\data_exp_170385-v5_task-na2a.csv';
raw_files(3).condition = 'Safe';
raw_files(3).group = 'D'; % D = threat first

raw_files(4).name = 'threat_first\data_exp_170385-v5_task-p4ov.csv';
raw_files(4).condition = 'Threat';
raw_files(4).group = 'D';


% Loop through these files
for iRaw = 1:4

    % Read table
    T = readtable([data_dir, raw_files(iRaw).name]);

    % get list of IDs
    IDs = unique(T.ParticipantPrivateID(~isnan(T.ParticipantPrivateID)));
    
    for iID = 1:numel(IDs)
        ID = IDs(iID);
        fprintf('\n\tID: %i', ID);
    
        sub_raw = T(T.ParticipantPrivateID == ID, :);
        sub_raw = sub_raw(strcmp(sub_raw.display, 'Trials'), :); % remove anything that's not a trial
    
        sub_table = table(...
            'Size', [192, numel(v_names)],...
            'VariableNames', v_names,...
            'VariableTypes', v_types);
    
        for trialN = 1:192
            % get raw data for this trial only
            trial_raw = sub_raw(sub_raw.TrialNumber == trialN, :);
    
            % get trial info
            blockN = trial_raw.SubBlockNumber(1);
            expectedness = trial_raw.Expectedness{1};
            contingency = trial_raw.ContingencyI_e_P_cue1_Happy_(1);
        
            % get cue
            cue_file = trial_raw.Cue{1}; % e.g. cue2_330Hz_500ms.mp3
            cue_split = split(cue_file, '_');
            if strcmp(cue_split{1}, 'cue1')
                cue_i = 0;
                cue_n = 'High';
            else
                cue_i = 1;
                cue_n = 'Low';
            end
    
            % get outcome
            outcome_file = trial_raw.Outcome{1}; % e.g. 'sad3.png'
            if strcmp(outcome_file(1:end-5), 'happy')
                outcome_i = 0;
                outcome_n = 'Happy';
                if strcmp(outcome_file(end-4), '1')
                    outcome_p_sad = 0;
                elseif strcmp(outcome_file(end-4), '2')
                    outcome_p_sad = 20;
                elseif strcmp(outcome_file(end-4), '3')
                    outcome_p_sad = 40;
                end
            else
                outcome_i = 1;
                outcome_n = 'Sad';
                if strcmp(outcome_file(end-4), '1')
                    outcome_p_sad = 60;
                elseif strcmp(outcome_file(end-4), '2')
                    outcome_p_sad = 80;
                elseif strcmp(outcome_file(end-4), '3')
                    outcome_p_sad = 100;
                end
            end
    
            % get response
            response_screen_idx = strcmp(trial_raw.ScreenName, 'Respond');
            if any(response_screen_idx)
                response_n = trial_raw(response_screen_idx, 'Response').Response{1};
                if strcmp(response_n, 'Happy')
                    response_i = 0;
                elseif strcmp(response_n, 'Sad')
                    response_i = 1;
                else
                    response_i = NaN;
                end
            
                % get response RT
                response_RT = trial_raw(strcmp(trial_raw.ScreenName, 'Respond'), 'ReactionTime').ReactionTime(1);
                if response_RT > 4000 % some people responded at 4003ms, so records response despite timeout
                    response_RT = NaN;
                end
            else
                response_n = '';
                response_i = NaN;
                response_RT = NaN;
            end
    
            % get confidence
            confidence_screen_idx = strcmp(trial_raw.ScreenName, 'Confidence');
            if any(confidence_screen_idx)
                confidence_n = trial_raw(confidence_screen_idx, 'Response').Response{1};
                if strcmp(confidence_n, 'High')
                    confidence_i = 1;
                elseif strcmp(confidence_n, 'Low')
                    confidence_i = 0;
                end
        
                % get confidence RT
                confidence_RT = trial_raw(strcmp(trial_raw.ScreenName, 'Confidence'), 'ReactionTime').ReactionTime(1);
            else
                confidence_n = '';
                confidence_i = NaN;
                confidence_RT = NaN;
            end
    
            % Add to table
            sub_table(trialN, 'Trial_N') = {trialN};
            sub_table(trialN, 'Block_N') = {blockN};
            sub_table(trialN, 'Contingency') = {contingency};
            sub_table(trialN, 'Expectedness') = {expectedness};
            sub_table(trialN, 'Cue_idx') = {cue_i};
            sub_table(trialN, 'Cue_name') = {cue_n};
            sub_table(trialN, 'Outcome_idx') = {outcome_i};
            sub_table(trialN, 'Outcome_name') = {outcome_n};
            sub_table(trialN, 'Outcome_p_sad') = {outcome_p_sad};
            sub_table(trialN, 'Response_name') = {response_n};
            sub_table(trialN, 'Response_idx') = {response_i};
            sub_table(trialN, 'Response_RT') = {response_RT};
            sub_table(trialN, 'Confidence_idx') = {confidence_i};
            sub_table(trialN, 'Confidence_name') = {confidence_n};
            sub_table(trialN, 'Confidence_RT') = {confidence_RT};
    
            % save table
            % writetable(sub_table,...
            %     ['..\..\DATA\STE_data2025\', ...
            %     num2str(ID), '_', ...
            %     raw_files(iRaw).group, '_', ...
            %     raw_files(iRaw).condition, '.csv']);
        end
    end
end

