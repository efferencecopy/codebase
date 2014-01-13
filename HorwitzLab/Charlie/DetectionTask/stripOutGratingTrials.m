function [stro gratingStruct] = stripOutGratingTrials(stro)

DTindicies;
gratingTrials = stro.trial(:,trialTypeInd) == 1;
gratingStruct.trial = stro.trial(gratingTrials,:);
gratingStruct.ras = stro.ras(gratingTrials,:);
stro.trial(gratingTrials,:)= [];
stro.ras(gratingTrials,:) = []; 