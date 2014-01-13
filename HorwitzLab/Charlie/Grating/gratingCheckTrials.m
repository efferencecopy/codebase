function [newtrlidxs, params] = gratingCheckTrials(dat, tuneFun, oldidxs)

% CHECK GRATING TRIALS
%
% EXAMPLE:  newTrialIndexs = gratingCheckTrials(stro, 'orient', oldTrialIndexes)
%
% The new gratings online gui has the ability to dynamically change the
% spatial frequency, diameter, and color. This means that any particular
% protocol run may be a combination of sf, diams, and colors. For example,
% a protocol 1 run (orientation) may have gratings of a few different
% diameters. This function is designed to identify the sf, diam, and color
% that occur most frequently, and return the idicies to these trials. Thus,
% a tuning function is computed on trials that vary one independent
% variable at a time.
%
% CAH 07.15.10


%spatial frequency
nSfs = length(unique(dat.sfs(oldidxs)));
sfs = unique(dat.sfs(oldidxs));
if nSfs>1
    for a = 1:nSfs
        l_sfs(:,a) = dat.sfs(oldidxs) == sfs(a);
    end
    [~, whichSF] = max(sum(l_sfs,1));
    params.sf = sfs(whichSF);
    l_sfs = l_sfs(:,whichSF);
else
    params.sf = sfs;
    l_sfs = dat.sfs(oldidxs) == sfs;
end

%orientation
nOrients = length(unique(dat.orients(oldidxs)));
orients = unique(dat.orients(oldidxs));
if nOrients>1
    for a = 1:nOrients
        l_orient(:,a) = dat.orients(oldidxs) == orients(a);
    end
    [~, whichOrient] = max(sum(l_orient, 1));
    params.orient = orients(whichOrient);
    l_orient = l_orient(:,whichOrient);
else
    params.orient = orients;
    l_orient = dat.orients(oldidxs) == orients;
end

%diameters
nDiams = length(unique(dat.diams(oldidxs)));
diams = unique(dat.diams(oldidxs));
if nDiams > 1
    for a = 1:nDiams
        l_diam(:,a) = dat.diams(oldidxs) == diams(a);
    end
    [~, whichDiam] = max(sum(l_diam,1));
    params.diam = diams(whichDiam);
    l_diam = l_diam(:,whichDiam);
else
    params.diam = diams;
    l_diam = dat.diams(oldidxs) == diams;
end

%color
nColors = size(unique(dat.colordirections(oldidxs,:), 'rows'),1);
colors = unique(dat.colordirections(oldidxs,:), 'rows');
if nColors > 1
    for a = 1:nColors
        l_color(:,a) = ismember(dat.colordirections(oldidxs,:), colors(a,:), 'rows');
    end
    [~, whichColor] = max(sum(l_color,1));
    params.color = colors(whichColor,:);
    l_color = l_color(:, whichColor);
else
    params.color = colors;
    l_color = ismember(dat.colordirections(oldidxs,:), colors, 'rows');
end

newtrlidxs = oldidxs; %in case there are no conflicts
switch tuneFun
    case 'orient'
        if any([nSfs, nDiams, nColors]>1)
            l_commonTrials = l_sfs & l_diam & l_color;
            newtrlidxs = oldidxs(l_commonTrials);
        end
    case 'sptFreq'
        if any([nOrients, nDiams, nColors]>1)
            l_commonTrials = l_orient & l_diam & l_color;
            newtrlidxs = oldidxs(l_commonTrials);
        end
    case 'diam'
        if  any([nOrients, nSfs, nColors]>1)
            l_commonTrials = l_orient & l_sfs & l_color;
            newtrlidxs = oldidxs(l_commonTrials);
        end
    case 'color'
        if any([nOrients, nSfs, nDiams]>1);
            l_commonTrials = l_orient & l_sfs & l_diam;
            newtrlidxs = oldidxs(l_commonTrials);
        end
end
