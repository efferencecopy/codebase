function [colors, sfs, data, tf] = questBatchProcess(txtFile, fitMeth, numTrials, perfRange)
%
% AGGREGATE DATA FROM MULTIPLE QUEST EXPERIMENTS
%
%  EXAMPLE:   [colors, sfs, data] = questBatchProcess(txtFile, fitMeth, [numTrials], [perfRange]);
%
% By default, this code will not exclude experiments on the basis of the
% observer's performance on the the last few trials. If this is desired,
% set numTrials and perfRange to the number of trials over which to compute
% percent correct, and to the tolerance level.
%
% CAH

if ~exist('numTrials', 'var')
    numTrials = 20;
end
if ~exist('perfRange', 'var')
    perfRange = [];
end


fnames = fnamesFromTxt2(txtFile);
nExpts = size(fnames, 1);
colors = [];  
sfs = [];
for a = 1:nExpts
    stro = nex2stro(findfile(fnames{a}{1}));
    disp(fnames{a});
    DTindicies;

    %unpack the data
    switch stro.sum.paradigmID
        case 210 %quest from DT
            [thresh, ~, exptSfs] = DTunpack(stro, lower(fitMeth), numTrials, perfRange);
            exptColors = reshape(stro.sum.exptParams.RF_colors, 3, 3)';
        case 211 %quest from HABIT
            [thresh, ~, exptSfs] = habitUnpack(stro, lower(fitMeth), numTrials, perfRange);
            exptColors = reshape(stro.sum.exptParams.gaborColor, 3, 3)';
    end
    %Since color directions can be represented in two (sign inverted) ways,
    %standardize the representation. The convention used will be to have
    %more + than - cone signs.
    minusLs = sign(exptColors(:,1))+eps;
    exptColors = exptColors.* repmat(minusLs, 1, 3);
    noColor = sum(abs(exptColors), 2) == 0;
    exptColors(noColor,:) = [];
    norms = sqrt(sum(exptColors.^2, 2));
    exptColors = exptColors./repmat(norms, 1, 3);
    
    
    set(gcf,'name',stro.sum.fileName(end-13:end-4)) %figure generated by DTunpack
    close(gcf)
    expt(a).t = thresh;
    expt(a).sfs = exptSfs;
    expt(a).col = exptColors;


    %aggregate data across experimental sessions
    for b = 1:length(exptSfs)
        if find(sfs == exptSfs(b), 1)
            continue
        else
            sfs(end+1) = exptSfs(b);
        end
    end

    for b = 1:size(exptColors, 1);
        sameColor = softEq(exptColors(b,:), colors, [], 'rows');
        if any(sameColor)
            continue
        else
            colors(end+1, :) = exptColors(b,:);
        end
    end
end

% sort the data into a big matrix of dimensions (colors x sfs x nExpts)
data = nan(size(colors, 1), length(sfs), nExpts);
for a = 1:nExpts
    for c = 1:size(expt(a).col, 1)
        rowInd = softEq(expt(a).col(c,:), colors, [], 'rows');
        for s = 1:length(expt(a).sfs)
            colInd = [sfs == expt(a).sfs(s)];
            data(rowInd, colInd, a) = expt(a).t(c, s);
        end
    end
end
