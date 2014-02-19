classdef dtobj
    % creates a DT object. Basically the same as a STRO structure except
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
        function obj = dtobj(fileName)
            if ~exist('fileName', 'var') || isempty(fileName)
                currentDir = pwd;
                if strcmpi(license, '380245') %charlie's laptop
                    cd('/Users/charliehass/LabStuff/NexFiles/Charlie');
                elseif ispc %lab PCs
                    cd('N:\NexFiles\Charlie');
                end
                [name,filepath] = uigetfile('*.nex');
                filepath = [filepath,name];
                cd(currentDir);
            else
                [tmppath, ~, tmpext] = fileparts(fileName);
                if ~isempty(tmppath) && strcmpi(tmpext, '.nex'); % the user supplied a valid absolute path
                    filepath = fileName;
                else
                    switch whoami
                        case 'hass_mbp'
                            load '/Users/charliehass/LabStuff/Huskies/NexFiles/Charlie/Batch Data And Text Files/nexPaths.mat';
                            filepath = findNexPath(nexPaths, fileName);
                        case 'nuke'
                            load '\\crash.dhe.duke.edu\charlie\DTcones\Files for server\nexPaths.mat'
                            filepath = findNexPath(nexPaths, fileName);
                        otherwise
                            filepath = findfile(fileName);
                    end
                end
            end
            
            
            %unpack the nexfile
            fprintf('filepath: %s\n', filepath);
            stro = nex2stro(filepath);
            if stro.sum.paradigmID ~= 210;
                error('The selected .nex file is not a DT expt')
            end
            
            %define the main fields
            stro = stripOutGratingTrials(stro);
            obj.sum = stro.sum;
            obj.trial = stro.trial;
            obj.ras = stro.ras;
            obj.other = stro.other;
            slashInd = strfind(stro.sum.fileName, filesep);
            obj.name = stro.sum.fileName(slashInd(end)+1 : end);
            
            %define the indicies
            obj.idx.actflashR  = strmatch('act_flash_R', [stro.sum.trialFields(1,:)]);
            obj.idx.actflashG = strmatch('act_flash_G', [stro.sum.trialFields(1,:)]);
            obj.idx.actflashB = strmatch('act_flash_B', [stro.sum.trialFields(1,:)]);
            obj.idx.choiceTime = strmatch('choice_time', [stro.sum.trialFields(1,:)]);
            obj.idx.colorDir = strmatch('color_dir', [stro.sum.trialFields(1,:)]);
            obj.idx.correct = strmatch('correct', [stro.sum.trialFields(1,:)]);
            obj.idx.cntrstLev = strmatch('cntrst_lev', [stro.sum.trialFields(1,:)]);
            obj.idx.driftRate = strmatch('gabor_speed', [stro.sum.trialFields(1,:)]);
            obj.idx.framesPresent = strmatch('frames_present', [stro.sum.trialFields(1,:)]);
            obj.idx.flashX = strmatch('flash_x', [stro.sum.trialFields(1,:)]);
            obj.idx.flashY = strmatch('flash_y', [stro.sum.trialFields(1,:)]);
            obj.idx.flashR = strmatch('flash_R', [stro.sum.trialFields(1,:)]);
            obj.idx.flashG = strmatch('flash_G', [stro.sum.trialFields(1,:)]);
            obj.idx.flashB = strmatch('flash_B', [stro.sum.trialFields(1,:)]);
            obj.idx.flashOn = strmatch('flash_on', [stro.sum.trialFields(1,:)]);
            obj.idx.repFlashOn = strmatch('rep_flash_on', [stro.sum.trialFields(1,:)]);
            obj.idx.flashOff = strmatch('flash_off', [stro.sum.trialFields(1,:)]);
            obj.idx.fpOn = strmatch('fp_on', [stro.sum.trialFields(1,:)]);
            obj.idx.frameOn = strmatch('frame_on', [stro.sum.trialFields(1,:)]);
            obj.idx.gaborGamma = strmatch('gabor_gamma', [stro.sum.trialFields(1,:)]);
            obj.idx.gaborLambda = strmatch('gabor_lambda', [stro.sum.trialFields(1,:)]);
            obj.idx.gaborSigma = strmatch('gabor_sigma', [stro.sum.trialFields(1,:)]);
            obj.idx.gaborTheta = strmatch('gabor_theta', [stro.sum.trialFields(1,:)]);
            obj.idx.nFrames = strmatch('numframes', [stro.sum.trialFields(1,:)]);
            obj.idx.questMode = strmatch('quest_thresh', [stro.sum.trialFields(1,:)]);
            obj.idx.rewOn = strmatch('rew_time', [stro.sum.trialFields(1,:)]);
            obj.idx.targOn = strmatch('targ_on', [stro.sum.trialFields(1,:)]);
            obj.idx.trialType = strmatch('trial_type', [stro.sum.trialFields(1,:)]);
            
            if ~isempty(obj.sum.rasterCells)
                obj.idx.spikes = strcmpi('sig001a', obj.sum.rasterCells(1,:));
                obj.idx.lfp = strcmpi('AD01', obj.sum.rasterCells(1,:));
                obj.idx.anlgStart = strcmpi('anlgStartTime', obj.sum.rasterCells(1,:));
            end
            
            %create the .LMS field
            obj.LMS = RGB2LMS(obj);
            
        end
        function ccs = RGB2LMS(obj)
            numTrials = size(obj.trial,1);
            bkgndrgb = [obj.sum.exptParams.bkgnd_r, obj.sum.exptParams.bkgnd_g, obj.sum.exptParams.bkgnd_b];
            M = reshape(obj.sum.exptParams.m_mtx, 3, 3);
            bkgndlms = M * bkgndrgb';
            x = 0:255; %the normal range of the gamma look up table
            xx = linspace(0, 255, 2^16); %the desired quantization of the gammaTable
            g1 = reshape(obj.sum.exptParams.gamma_table, 256, 3);
            gammaTable = [spline(x, g1(:,1), xx)', spline(x, g1(:,2), xx)', spline(x, g1(:,3), xx)'];
            gunVals = [obj.trial(:,obj.idx.flashR), obj.trial(:,obj.idx.flashG), obj.trial(:,obj.idx.flashB)] + 1;
            rgb = [gammaTable(gunVals(:,1), 1), gammaTable(gunVals(:,2), 2), gammaTable(gunVals(:,3), 3)];
            rgb = rgb - repmat(bkgndrgb, numTrials, 1);
            ccs = [(M * rgb') ./ repmat(bkgndlms, 1, numTrials)]';
            ccs(ccs(:,1)<0,:) = -ccs(ccs(:,1)<0,:); % flip the sign on the L-cones if it's negative
            l_siso = softEq([0,0,1], abs(ccs), 4, 'rows');
            if any(sum(ccs(l_siso,:),2) < 0);
                ccs(l_siso,:) = -ccs(l_siso,:);
            end
        end
        
        
        function resp = spikeResp(obj, type)
            t_on = mat2cell(obj.trial(:, obj.idx.repFlashOn), ones(size(obj.trial,1),1));
            t_off = mat2cell(obj.trial(:, obj.idx.flashOff), ones(size(obj.trial,1), 1));
            counts = cellfun(@(on,off,sp) sum((sp>on)&(sp<=off)), t_on, t_off, obj.ras(:, obj.idx.spikes));
            if strcmpi(type, 'count')
                resp = counts;
            elseif strcmpi(type, 'rate')
                duration = cellfun(@minus, t_off, t_on);
                resp = counts./duration; %spikes/sec
            else
                error('input <%s> not recognized', type)
            end
        end
    end %methods
    
end

