classdef gtobj
    % creates a GT object. 
    %   Detailed explanation goes here
    
    properties
        sum = [];
        trial = [];
        ras = [];
        other = [];
        idx = [];
        name = [];
    end
    
    methods
        function obj = gtobj(fileName)
            if ~exist('fileName', 'var') || isempty(fileName)
                currentDir = pwd;
                if isunix
                    cd('/Users/charliehass/LabStuff/NexFiles/Charlie');
                elseif ispc
                    cd('N:\NexFiles\Charlie');
                end
                [file_name,file_path] = uigetfile('*.nex');
                file_path = [file_path file_name];
                cd(currentDir);
            else
                if ismac
                    load('/Users/charliehass/LabStuff/NexFiles/Charlie/Batch Data And Text Files/nexPaths.mat');
                elseif strcmpi(license, '367516')
                    load('C:\NO BACKUP\NexFiles\Charlie\Batch Data And Text Files\nexPaths.mat');
                elseif ispc
                    load('N:\NexFiles\Charlie\Batch Data And Text Files\nexPaths.mat');
                end
                file_path = findNexPath(nexPaths, fileName);
            end
            
            %unpack the nexfile
            stro = nex2stro(file_path);
            if stro.sum.paradigmID ~= 150;
                error('The selected .nex file is not a GT expt')
            end

            %define the main fields
            obj.sum = stro.sum;
            obj.trial = stro.trial;
            obj.ras = stro.ras;
            obj.other = stro.other;
            obj.name = fileName;
            
            %define indicies that are useful for analysis, but change on a
            %expt by expt basis.
            obj.idx.protocol = strcmpi('protocol', obj.sum.trialFields(1,:));
            obj.idx.orient = strcmpi('orient', obj.sum.trialFields(1,:));
            obj.idx.lcc = strcmpi('lcont', obj.sum.trialFields(1,:));
            obj.idx.mcc = strcmpi('mcont', obj.sum.trialFields(1,:));
            obj.idx.scc = strcmpi('scont', obj.sum.trialFields(1,:));
            obj.idx.diam = strcmpi('diam', obj.sum.trialFields(1,:));
            obj.idx.fpon = strcmpi('fp_on', obj.sum.trialFields(1,:));
            obj.idx.sf = strcmpi('sf', obj.sum.trialFields(1,:));
            obj.idx.fpacq = strcmpi('fp_acq', obj.sum.trialFields(1,:));
            obj.idx.stimon = strcmpi('stim_on', obj.sum.trialFields(1,:));
            obj.idx.stimoff = strcmpi('stim_off', obj.sum.trialFields(1,:));
            obj.idx.nFrames = strcmpi('nframes', obj.sum.trialFields(1,:));
            obj.idx.lfp = strcmpi('AD01', obj.sum.rasterCells(1,:));
            obj.idx.anlgStart = strcmpi('anlgStartTime', obj.sum.rasterCells(1,:));
            obj.idx.spikes = strcmpi('sig001a', obj.sum.rasterCells(1,:));
            obj.idx.wf = strcmpi('sig001a_wf', obj.sum.rasterCells(1,:));
            
            %standardize the colors contained in the .trial field. The L
            %cones should have a positive sign. S-iso should also be
            %positive.
            LMS = [obj.trial(:,obj.idx.lcc), obj.trial(:,obj.idx.mcc), obj.trial(:,obj.idx.scc)];
            l_negL = LMS(:,1) < 0;
            LMS(l_negL,:) = LMS(l_negL,:) .* -1; %deal with all non-siso colors
            l_negSiso = ismember(sign(LMS), [0 0 -1], 'rows');
            LMS(l_negSiso,:) = LMS(l_negSiso,:) .* -1; %deal with the Siso colors
            obj.trial(:,obj.idx.lcc) = LMS(:,1);
            obj.trial(:,obj.idx.mcc) = LMS(:,2);
            obj.trial(:,obj.idx.scc) = LMS(:,3);
        end
        function list = getTrialList(obj, protocolType, color)
            LMS = [obj.trial(:, obj.idx.lcc), obj.trial(:, obj.idx.mcc), obj.trial(:, obj.idx.scc)];
            l_p = obj.trial(:,obj.idx.protocol) == protocolType;
            if exist('color', 'var')
            l_color = ismember(sign(LMS), color, 'rows');
            else
                l_color = true(size(LMS,1),1);
            end
            list = l_p & l_color;            
        end
        function resp = spikeResp(obj, type)
            t_on = mat2cell(obj.trial(:, obj.idx.stimon), ones(size(obj.trial,1),1));
            t_off = mat2cell(obj.trial(:, obj.idx.stimoff), ones(size(obj.trial,1), 1));
            counts = cellfun(@(on,off,sp) sum((sp>on)&(sp<=off)), t_on, t_off, obj.ras(:, obj.idx.spikes)); %#ok<CPROP>
            if strcmpi(type, 'count')
                resp = counts;
            elseif strcmpi(type, 'rate')
                duration = cellfun(@minus, t_off, t_on);
                resp = counts./duration; %spikes/sec
            else
                error('input <%s> not recognized', type)
            end                
        end
    end
    
end

