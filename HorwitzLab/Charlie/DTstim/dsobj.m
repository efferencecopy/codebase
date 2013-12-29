classdef dsobj
    % creates a DTstim object. Basically the same as a STRO structure except
    % that a few other things are defined here. Indicies to trial events
    % are defined here too.
    
    properties
        sum = [];
        trial = [];
        ras = [];
        other = [];
        idx = [];
        LMS = [];
        name = [];
    end
    
    methods
        function obj = dsobj(fileName)
            if ~exist('fileName', 'var') || isempty(fileName)
                currentDir = pwd;
                if isunix
                    cd('/Users/charliehass/LabStuff/NexFiles/Charlie');
                elseif ispc
                    cd('N:\NexFiles\Charlie');
                end
                [name,filepath] = uigetfile('*.nex');
                filepath = [filepath,name];
                cd(currentDir);
            else
                if isunix
                    load '/Users/charliehass/LabStuff/NexFiles/Charlie/Batch Data And Text Files/nexPaths.mat';
                elseif ispc
                    load 'N:\NexFiles\Charlie\Batch Data And Text Files\nexPaths.mat';
                end
                filepath = findNexPath(nexPaths, fileName);
            end
            
            %unpack the nexfile
            stro = nex2stro(filepath);
            if stro.sum.paradigmID ~= 212;
                error('The selected .nex file is not a DTstim expt')
            end
            
            %define the main fields
            obj.sum = stro.sum;
            obj.trial = stro.trial;
            obj.ras = stro.ras;
            obj.other = stro.other;
            slashInd = strfind(stro.sum.fileName, filesep);
            obj.name = stro.sum.fileName(slashInd(end)+1 : end);
            
            %define the indicies
            for a = 1:size(obj.sum.trialFields,2);
                obj.idx = setfield(obj.idx, obj.sum.trialFields{1,a}, a);
            end
            
            %need to interate over all the available continuous channels.
            if ~isempty(obj.sum.rasterCells)
                obj.idx.spikes = strncmpi('sig', obj.sum.rasterCells(1,:), 3);
                obj.idx.anlg = strncmpi('AD', obj.sum.rasterCells(1,:),2);
                obj.idx.anlgStart = strcmpi('anlgStartTime', obj.sum.rasterCells(1,:));
            end
            
            %convert the RGB vals to LMS
            obj.LMS = RGB2LMS(obj);
        end
        
        function ccs = RGB2LMS(obj)
            numTrials = size(obj.trial,1);
            bkgndrgb = [obj.sum.exptParams.bkgndr, obj.sum.exptParams.bkgndg, obj.sum.exptParams.bkgndb];
            M = reshape(obj.sum.exptParams.mMtx, 3, 3);
            bkgndlms = M * bkgndrgb(:);
            x = 0:255; %the normal range of the gamma look up table
            xx = linspace(0, 255, 2^16); %the desired quantization of the gammaTable
            g1 = reshape(obj.sum.exptParams.gammaTable, 256, 3);
            gammaTable = [spline(x, g1(:,1), xx)', spline(x, g1(:,2), xx)', spline(x, g1(:,3), xx)'];
            gunVals = [obj.trial(:,obj.idx.Rgun), obj.trial(:,obj.idx.Ggun), obj.trial(:,obj.idx.Bgun)] + 1;
            rgb = [gammaTable(gunVals(:,1), 1), gammaTable(gunVals(:,2), 2), gammaTable(gunVals(:,3), 3)];
            rgb = rgb - repmat(bkgndrgb, numTrials, 1);
            ccs = [(M * rgb') ./ repmat(bkgndlms, 1, numTrials)]';
        end
    end
end
