function pal_BIC_Fit_comparison(theModels,fileSpec)

%reads in Fit structs saved in .mat files

%input:
% theModels: should be a vertical cell containing strings for model names
    %example: 
        %{'PH'}
        %{'RW'}
        %{'PH_pw'}
    %although a string array should also work, ["PH","RW","PH_pw"]';
    %IMPORTANT: it must be vertical
% fileSpec: input to sprintf for getting .mat file containing the Fit
    
if ~exist('fileSpec','var')
    fileSpec = 'Fitted_%s.mat'; %default setting
end

% Display the resulting cell array
disp(theModels);

Nmodels = length(theModels);

meanBIC = nan(Nmodels,1);
seBIC = nan(Nmodels,1);

for i = 1:Nmodels
    %load the fitted struct
    Fit = load(sprintf(fileSpec,theModels{i}));
    Fit = Fit.Fit;
    
    %
    subNum = Fit.bestFit(:,1); %subject number
    
    %
    allBICs(i,:) = Fit.BIC; %column will be n_subjects; each row a model
        %IMPORTANT: in different row/column to match VBA requirement
    meanBIC(i,:) = mean(Fit.BIC);
    seBIC(i,:) = std(Fit.BIC)/sqrt(length(Fit.BIC));
     
end

tBIC = array2table(allBICs');
tBIC.Properties.VariableNames = theModels;
tBIC.subNum = subNum;

%% doing the aggregate
%bms:
lme = -allBICs/2;
[posterior,out] = VBA_groupBMC(lme);
%[~,~,~,pxp,~,~] = bms(lme); %this is Gershman's old code
pxp = out.pxp;
%size(pxp)
%size(theModels)
%T = array2table(pxp, 'VariableNames', theModels);
T = table(theModels,pxp', 'VariableNames', {'model','PXP'});
disp('Full sample:')
disp(T)

T_modelFreq = table(theModels,out.Ef, 'VariableNames', {'model','estimated model frequencies'});
disp(T_modelFreq)

%visual
figure;
bar(meanBIC)
hold on
errorbar(1:length(theModels), meanBIC, seBIC, '.', 'LineWidth', 1)
hold off
legend('mean BIC across ptp','standard error','Location','Southwest')
xticks(1:length(theModels))
set(gca,'TickLabelInterpreter','none')
xticklabels(theModels)
xtickangle(90)
xlabel('Models')
ylabel('Mean BIC')
title('Comparison based on mean BIC across participants')
ylim([min(meanBIC),max(meanBIC)])

end