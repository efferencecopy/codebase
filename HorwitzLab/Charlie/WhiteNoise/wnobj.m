classdef wnobj
    % creates a WN object. Basically the same as a STRO structure except
    % that a few other things are defined here. Indicies to trial events
    % are defined here too.
    properties
        sum = [];
        trial = [];
        ras = [];
        other = [];
        idx = [];
        name = [];
    end
    
    methods
        function obj = wnobj(fileName)
            if ~exist('fileName', 'var') || isempty(fileName)
                currentDir = pwd;
                if isunix
                    cd('/Users/charliehass/LabStuff/NexFiles/Charlie');
                elseif ispc
                    cd('N:\NexFiles\Charlie');
                end
                [diskname,filepath] = uigetfile('*.nex');
                filepath = [filepath,diskname];
                cd(currentDir);
            else
                if strcmp(license, '367516')
                    load 'C:\NO BACKUP\NexFiles\Charlie\Batch Data And Text Files\nexPaths.mat';
                    filepath = findNexPath(nexPaths, fileName);
                elseif ispc %the non-server lab PCs
                    load 'N:\NexFiles\Charlie\Batch Data And Text Files\nexPaths.mat';
                    filepath = findNexPath(nexPaths, fileName);
                elseif strcmpi(license, '380245') %charlie's laptop
                    load '/Users/charliehass/LabStuff/NexFiles/Charlie/Batch Data And Text Files/nexPaths.mat';
                    filepath = findNexPath(nexPaths, fileName);
                elseif strcmpi(license, '633208') % rig 3 mac
                    load '~/Desktop/Charlie Cone Model Stuff/nexPaths.mat';
                    filepath = findNexPath(nexPaths, fileName);
                else
                    filepath = findfile(fileName);
                end
            end
            
            %unpack the nexfile
            stro = nex2stro(filepath);
            if stro.sum.paradigmID ~= 100
                error('The selected .nex file is not a WN expt')
            end
            
            %define the main fields
            obj.sum = stro.sum;
            obj.trial = stro.trial;
            obj.ras = stro.ras;
            obj.other = stro.other;
            slashInd = strfind(stro.sum.fileName, filesep);
            obj.name = stro.sum.fileName(slashInd(end)+1 : end);
            
            %define the indicies
            for a = 1:size(obj.sum.trialFields,2)
                obj.idx.(obj.sum.trialFields{1,a}) = a;
            end
            
            %need to interate over all the available continuous and spike channels.
            if ~isempty(obj.sum.rasterCells)
                obj.idx.spikes = strncmpi('sig', obj.sum.rasterCells(1,:), 3);
                obj.idx.eyepos = ismember(obj.sum.rasterCells(1,:),{'AD11' 'AD12'});
                obj.idx.lfp = strncmpi('AD', obj.sum.rasterCells(1,:),2) & ~obj.idx.eyepos;
                obj.idx.anlgStart = strcmpi('anlgStartTime', obj.sum.rasterCells(1,:));
            end
        end
        
        function oob = checkTrials(obj, lfpidx)
            nADSteps = 2^12;
            anlgChannelName = obj.sum.rasterCells{lfpidx};
            anlgChannelIdx = strcmpi(obj.sum.analog.sigid, anlgChannelName);
            ADtoMV = obj.sum.analog.ADtoMV{anlgChannelIdx};
            maxVal = (ADtoMV * nADSteps) ./ 2; %it's actually +/- this value;
            oob = cellfun(@(x) any(abs(x)>=maxVal), obj.ras(:, lfpidx));
        end
    end
end
