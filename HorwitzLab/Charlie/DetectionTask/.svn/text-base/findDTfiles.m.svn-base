function fnames = findDTfiles(topDirectory, fnames)

    %   FIND VALID DT FILES
    %
    % example: fnames = findDTfiles

    %the minimum number of trials (in RF) to consider a valid expt
    MINNUMTRIALS = 8;

    %find a starting directory. When called recursively, this is already
    %defined.
    if ~exist('topDirectory', 'var')
        if isunix
            topDirectory = uigetdir('/Users/charliehass/LabStuff/NexFiles/');
        elseif ispc
            topDirectory = uigetdir('C:\NO BACKUP\NexFiles');
        end
    end


    %initialize a list of files, but only if the command window calls the
    %function. The recursive calls will define 'fnames', and hence NOT
    %initialize the file list
    if ~exist('fnames', 'var')
        fnames = {};
    end


    tree = dir(topDirectory);
    possibleDirs = {};
    for a = 1:length(tree);
        if tree(a).isdir && ~strcmp(tree(a).name, '.') && ~strcmp(tree(a).name, '..')
            possibleDirs{end+1} = [topDirectory, '/', tree(a).name];
        elseif ~tree(a).isdir
            if strcmpi(tree(a).name(end-3:end), '.nex')
                try
                    [topDirectory, '/', tree(a).name] %display the current file name to the screen
                    DT = dtobj(tree(a).name);
                    [~, cell, ~] = DTunpack(DT, 1);
                    nTrialsInRF = min(horzcat(cell.nTrialsIn{:}));
                    close all hidden %all the figures that have been generated (including the pesky waitbars...
                    if nTrialsInRF >= MINNUMTRIALS;
                        fnames{end+1} = tree(a).name
                    end
                catch
                    close all hidden %all the figures that have been generated (including the pesky waitbars...
                    continue %i.e., probably not a DT file
                end
            end
        end
    end

    % the recursive part
    for a = 1:length(possibleDirs)
        fnames = findDTfiles(possibleDirs{a}, fnames);
    end
end








