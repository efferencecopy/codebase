%% NOTE FOR HOW TO CREATE SPATIAL CONTRAST SENSITIVITY PLOTS IDENTICAL TO THOSE PRODUCED FOR KALI

% PICK AN OBSERVER:
%observer = '/Users/charliehass/LabStuff/NexFiles/Charlie/Kali/text files/questCSFdata.txt' %kali
% observer = 'N:\NexFiles\Charlie\Kali\text files\questCSFnewchamber.txt'; %kali, new axis of action
 observer = '/Users/charliehass/LabStuff/NexFiles/Charlie/Sedna/text files/quest.txt'; %sedna
% observer = '/Users/charliehass/LabStuff/NexFiles/Charlie/Apollo/text files/questTraining.txt'; %Apollo
% observer = 'N:\NexFiles\Charlie\CharliePsychophysics\Text Files\charlieTrainingSet.txt'; %charlie


% STEP ONE: modify the textfile list so that all the data files have the
% same stimulus parameters. Use paramsCheck. Make sure that all the quest
% functions have at least 30 trials per run.
% p = paramsCheck(observer)

% STEP TWO: batch process the data:
nTrials = 20;
perfRange = []; %don't filter on the baisis of performance
[colors, sfs, data] = questBatchProcess(observer, 'mode', nTrials, perfRange);

% STEP THREE: remove spatial frequencies that don't fall on the traditional latice
sfs %find the non-traditional values just to double check
ind = softEq(0.8893, sfs, 4) | softEq(0.9919, sfs, 4)
sfs(ind) = [];
data(:,ind,:) = [];

% STEP FOUR: plot the CSF and make sure that there are enough experiments
% per condition
order = 2;
figure
f = fitCSF(colors, sfs, data, order, 'all');
nExpts = sum(~isnan(data), 3);


