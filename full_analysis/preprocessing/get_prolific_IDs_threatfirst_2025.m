
T = readtable('..\..\DATA\2025\raw\threat_first\data_exp_170385-v5_questionnaire-y4ui.csv');
private_IDs = unique(T.ParticipantPrivateID(~isnan(T.ParticipantPrivateID)));
prolific_IDs = arrayfun(@(x) T(T.ParticipantPrivateID == private_IDs(x) & strcmp(T.ResponseType, 'response'), :).Response{1}, 1:numel(private_IDs), 'UniformOutput', false)';


fname='completed_IDs_threatfirst_2025.txt';
fid=fopen(fname, 'w');
for i=1:numel(prolific_IDs)
    fprintf(fid, '%s\n', prolific_IDs{i});
end
fclose(fid);