function [trace, info] = fiberVolleyAnalysis(exptList, exptWorkbook, PLOTFIGURES, RMLINENOISE, EXPTTYPE)


if ~exist('PLOTFIGURES', 'var')
    PLOTFIGURES = true;
end



% figure out the file names, etc.
mouseNames = exptWorkbook(exptList, strcmpi(exptWorkbook(1,:), 'Mouse Name'));
siteNumber = exptWorkbook(exptList, strcmpi(exptWorkbook(1,:), 'site'));
fnames = exptWorkbook(exptList, strcmpi(exptWorkbook(1,:), 'file name'));
exptConds = exptWorkbook(exptList, strcmpi(exptWorkbook(1,:), 'drugs'));
channels = exptWorkbook(exptList, cellfun(@(x) ~isempty(x), regexpi(exptWorkbook(1,:), 'CH\d')));
rmsweeps = exptWorkbook(exptList, strcmpi(exptWorkbook(1,:), 'rmSweeps'));



% ERROR CHECKING: check to make sure that the exptConds are correct
validConds = {'none',...
              'cd2',...
              'nbqx_apv',...
              'nbqx_apv_cd2',...
              'nbqx_apv_cd2_ttx',...
              'nbqx_apv_ttx',...
              'ttx',...
              'cd2_ttx'};
          
exptConds = cellfun(@(x) regexprep(x, '\+', '_'), exptConds, 'uniformoutput', false);
idx = cellfun(@(x,y) regexpi(x,y), exptConds, repmat({validConds}, size(exptConds)), 'uniformoutput', false);
idx = cellfun(@(x) any([x{:}]), idx);
assert(all(idx), 'ERROR: at least one experimental condition is not recognized');

% ERROR CHECKING: make sure that the tmp.head.validChans is consistent
% in the chanels used for each input file
channelConfigs = cell2mat(channels);
channelConfigs = unique(channelConfigs, 'rows');
channelList = find(channelConfigs);
assert(size(channelConfigs,1)==1, 'ERROR: The recorded channels changes from file to file')



% read in all the data
clc
fprintf('reading in the data for fiber volley analysis\n')
for i_fid = 1:numel(fnames)
    fprintf('   reading <%s>\n', fnames{i_fid})
    tmp = abfobj(fnames{i_fid});
    tmp.head.validChans = cat(1, channels{i_fid,:});
    
    % remove some sweeps if need be
    if ~isnan(rmsweeps{i_fid})
        goodSweeps = true(size(tmp.dat,3),1);
        badSweeps = eval(['[',rmsweeps{i_fid},']']);
        goodSweeps(badSweeps) = false;
        tmp.dat = tmp.dat(:,:,goodSweeps);
        tmp.wf = tmp.wf(:,:,goodSweeps);
    end
    
    % store the data, grouped by unique conditions. First, I need to figure
    % out the appropriate field name for each condition, to do this, I need
    % to figure out which channel has the optical stimulation, and then run
    % 'outerleave.m'
    [ledidx, laseridx] = deal([]);
    if isfield(tmp.idx, 'LED_470')
        ledidx = tmp.idx.LED_470;
    end
    if isfield(tmp.idx, 'Laser')
        laseridx = tmp.idx.Laser;
    end
    both_idx_defined = any(ledidx) && any(laseridx);
    assert(~both_idx_defined, 'ERROR: LED and Laser were used simultaneously');
    
    if any(ledidx)
        optostimidx = ledidx;
    elseif any(laseridx)
        optostimidx = laseridx;
    end
    assert(numel(optostimidx) == 1, 'ERROR: too many optostim channels defined')
    
    tmpWF = tmp.dat(:, optostimidx, :);
    sampRate = tmp.head.sampRate;
    tdict = outerleave(tmpWF, sampRate);
    nSweepTypes = size(tdict.conds,1);
    
    
    % what types of params are around?
    fileAmps = round(tdict.conds(:,1));
    Namps = numel(unique(fileAmps));
    fileWidths = round(tdict.conds(:,2));
    Nwidths = numel(unique(fileWidths));
    fileTFs = round(tdict.conds(:,3));
    fileRecovs = round(tdict.conds(:,4)); % in ms
    Nrecovs = numel(unique(fileRecovs));
    
    % display what things are interleaved
    fprintf('     Interleaved tfs: %s\n', num2str(unique(fileTFs)'));
    if Namps > 1
        fprintf('     Interleaved pAmps: %s\n', num2str(unique(fileAmps)'));
    end
    if Nwidths > 1
        fprintf('     Interleaved pWidths: %s\n', num2str(unique(fileWidths)'));
    end
    if Nrecovs > 1
        fprintf('     Interleaved recovery times: %s\n', num2str(unique(fileRecovs)'));
    end
    
    
    for i_sweepType = 1:nSweepTypes
        
        % new field names for the 'sweepType'
        if all([Namps, Nwidths, Nrecovs] == 1); % nothing but TF
            swpType = ['tf', num2str(fileTFs(i_sweepType))];
            
        elseif Namps>1 && all([Nwidths, Nrecovs] == 1); % multiple amps
            swpType = ['tf', num2str(fileTFs(i_sweepType)), '_amp', num2str(fileAmps(i_sweepType))];
            
        elseif Nwidths>1 && all([Namps, Nrecovs] == 1); % multiple pWidths
            swpType = ['tf', num2str(fileTFs(i_sweepType)), '_width', num2str(fileWidths(i_sweepType))];
            
        elseif Nrecovs>1 && all([Namps, Nwidths] == 1); % multiplerecovery times
            swpType = ['tf', num2str(fileTFs(i_sweepType)), '_recov', num2str(fileRecovs(i_sweepType))];
        end
        
        drugType = exptConds{i_fid};
        
        % make a franken-abf struct
        ax.(swpType).(drugType).head = tmp.head;
        ax.(swpType).(drugType).dat = tmp.dat(:,:,tdict.trlList == i_sweepType);
        ax.(swpType).(drugType).idx = tmp.idx;
        ax.(swpType).(drugType).pAmp = tdict.conds(i_sweepType,1);
        ax.(swpType).(drugType).pWidth = tdict.conds(i_sweepType,2);
        ax.(swpType).(drugType).pTF = fileTFs(i_sweepType);
        ax.(swpType).(drugType).tRecov = fileRecovs(i_sweepType);
        ax.(swpType).(drugType).realTrialNum = find(tdict.trlList == i_sweepType);
    end

end


% pull out the raw data, filter, and organize all the raw traces.
fprintf('Filtering and organizing raw data\n');
trace = [];
sweepTypeFields = fieldnames(ax);
for i_sweepType = 1:numel(sweepTypeFields)
    
    swpType = sweepTypeFields{i_sweepType};
    drugFields = fieldnames(ax.(swpType));
    
    for i_drug = 1:numel(drugFields)
        
        % generate some field names for extraction and saving
        drugType = drugFields{i_drug};
        
        sampFreq = ax.(swpType).(drugType).head.sampRate;
        %tt = (0:size(ax.(tType).(drugType).dat, 1)-1) ./ sampFreq .* 1000;
        
        
        % grab the raw data (channel by channel)
        validChans = ax.(swpType).(drugType).head.validChans;
        validChans = find(validChans == 1);
        for i_ch = 1:numel(validChans);
            
            % figure out the appropriate index to the data. It should have
            % the correct "HSx_" prefix, and have the correct units
            whichChan = validChans(i_ch);
            switch EXPTTYPE
                case 'Intracellular'
                    correctUnits = strcmpi('pA', ax.(swpType).(drugType).head.recChUnits);
                otherwise
                    correctUnits = strcmpi('mV', ax.(swpType).(drugType).head.recChUnits);
            end
            correctHS = strncmpi(sprintf('HS%d_', whichChan), ax.(swpType).(drugType).head.recChNames, 3);
            chIdx = correctUnits & correctHS;
            assert(sum(chIdx)==1, 'ERROR: incorrect channel selection')
            
            tmp = ax.(swpType).(drugType).dat(:, chIdx,:);
            tmp = squeeze(tmp);
            
            
            % baseline subtract, determine when the pulses came on...
            thresh = 0.05;
            [ledidx, laseridx] = deal([]);
            if isfield(ax.(swpType).(drugType).idx, 'LED_470')
                ledidx = ax.(swpType).(drugType).idx.LED_470;
            end
            if isfield(ax.(swpType).(drugType).idx, 'Laser')
                laseridx = ax.(swpType).(drugType).idx.Laser;
            end
            both_idx_defined = any(ledidx) && any(laseridx);
            assert(~both_idx_defined, 'ERROR: LED and Laser were used simultaneously');
            
            if any(ledidx)
                optostimidx = ledidx;
            elseif any(laseridx)
                optostimidx = laseridx;
            end
            assert(numel(optostimidx) == 1, 'ERROR: too many optostim channels defined')

            template = ax.(swpType).(drugType).dat(:,optostimidx,1);
            template = template > thresh;
            template = [0; diff(template)];
            storedCrossings_on{i_sweepType} = template == 1;
            tmp_off = template == -1;
            storedCrossings_off{i_sweepType} = [tmp_off(2:end); false]; % i think there is an OBO error otherwise...
            pulseOnset = find(template==1, 1, 'first');
            
            bkgndTime = min([0.100, pulseOnset./sampFreq]); % guards against there not being enough data in front of the frist pulse...  
            bkgndSamps = ceil(bkgndTime .* sampFreq);
            tmp = bsxfun(@minus, tmp, mean(tmp(pulseOnset-bkgndSamps:pulseOnset-1, :),1));
            
            
            % filter out the high frequency stuff. Filtering shouldn't go
            % below 1500 b/c you'll start to carve out the fiber volley
            % (which is only a few ms wide)
            lp_freq = 1500;
            filtered = butterfilt(tmp, lp_freq, sampFreq, 'low', 1);
            filtered = bsxfun(@minus, filtered, mean(filtered(pulseOnset-bkgndSamps:pulseOnset-1, :),1));
            
            % take the mean
            average = mean(filtered,2);
            trace.(swpType).(drugType)(:,i_ch) = average;
            
        end % i_ch
        
        % store the pulse onset times for the population analysis
        info.(swpType).(drugType).realTrialNum = ax.(swpType).(drugType).realTrialNum;
        info.(swpType).(drugType).pulseOn_idx = storedCrossings_on{i_sweepType};
        info.(swpType).(drugType).pulseOff_idx = storedCrossings_off{i_sweepType};
        info.(swpType).(drugType).sampRate = ax.(swpType).(drugType).head.sampRate;
        info.(swpType).(drugType).pWidth = ax.(swpType).(drugType).pWidth;
        info.(swpType).(drugType).pAmp = ax.(swpType).(drugType).pAmp;
        info.(swpType).(drugType).pTF = ax.(swpType).(drugType).pTF;
        
    end % i_drug
end % i_swpType



fprintf('Calculating fiber volley\n')
sweepTypeFields = fieldnames(ax);
for i_sweepType = 1:numel(sweepTypeFields)
    
    % rules:
    %
    % nbqx_apv - nbqx_apv_ttx           =>   fiber volley with Na and Ca2 and presynpatic mGluR
    % nbqx_apv - nbqx_apv_cd_ttx        =>   fiber volley with Na and Ca2 and presynpatic mGluR
    % nbqx_apv_cd2 - nbqx_apv_cd2_ttx   =>   fiber volley with just Na
    % none - nbqx_apv                   =>   synapticTransmission
    % ttx - cd2_ttx                     =>   direct release from terminals?
    % none                              =>   control
    
    
    % generate some field names
    swpType = sweepTypeFields{i_sweepType};
    
    % is there control data and nbqx_apv?
    control_Present = isfield(trace.(swpType), 'none');
    nbqx_apv_Present = isfield(trace.(swpType), 'nbqx_apv');
    if control_Present && nbqx_apv_Present
        
        pWidthMatch = info.(swpType).none.pWidth == info.(swpType).nbqx_apv.pWidth;
        pAmpMatch = info.(swpType).none.pAmp == info.(swpType).nbqx_apv.pAmp;
        assert(pWidthMatch && pAmpMatch, 'ERROR: pulse amp or width are not consistent');
        
        trace.(swpType).synapticTransmission = trace.(swpType).none - trace.(swpType).nbqx_apv;
        
        info.(swpType).synapticTransmission.pulseOn_idx = info.(swpType).none.pulseOn_idx;
        info.(swpType).synapticTransmission.pulseOff_idx = info.(swpType).none.pulseOff_idx;
        info.(swpType).synapticTransmission.sampRate = info.(swpType).none.sampRate;
        info.(swpType).synapticTransmission.pWidth = info.(swpType).none.pWidth;
        info.(swpType).synapticTransmission.pAmp = info.(swpType).none.pAmp;
        info.(swpType).synapticTransmission.pTF = info.(swpType).none.pTF;
        info.(swpType).synapticTransmission.realTrialNum = info.(swpType).none.realTrialNum;
    end
    
    % is there a ttx condition that can be used to define fiber volley with
    % sodium and calcium currents?
    nbqx_apv_ttx_Present = isfield(trace.(swpType), 'nbqx_apv_ttx');
    nbqx_apv_cd2_ttx_Present = isfield(trace.(swpType), 'nbqx_apv_cd2_ttx');
    assert(~all([nbqx_apv_ttx_Present, nbqx_apv_cd2_ttx_Present]), 'ERROR: multiple TTX files defined');
    if nbqx_apv_Present && nbqx_apv_ttx_Present
        
        pWidthMatch = info.(swpType).nbqx_apv.pWidth == info.(swpType).nbqx_apv_ttx.pWidth;
        pAmpMatch = info.(swpType).nbqx_apv.pAmp == info.(swpType).nbqx_apv_ttx.pAmp;
        assert(pWidthMatch && pAmpMatch, 'ERROR: pulse amp or width are not consistent');
        
        trace.(swpType).FV_Na_Ca2_mGluR = trace.(swpType).nbqx_apv - trace.(swpType).nbqx_apv_ttx;
        
        info.(swpType).FV_Na_Ca2_mGluR.pulseOn_idx = info.(swpType).nbqx_apv.pulseOn_idx;
        info.(swpType).FV_Na_Ca2_mGluR.pulseOff_idx = info.(swpType).nbqx_apv.pulseOff_idx;
        info.(swpType).FV_Na_Ca2_mGluR.sampRate = info.(swpType).nbqx_apv.sampRate;
        info.(swpType).FV_Na_Ca2_mGluR.pWidth = info.(swpType).nbqx_apv.pWidth;
        info.(swpType).FV_Na_Ca2_mGluR.pAmp = info.(swpType).nbqx_apv.pAmp;
        info.(swpType).FV_Na_Ca2_mGluR.pTF = info.(swpType).nbqx_apv.pTF;
        info.(swpType).FV_Na_Ca2_mGluR.realTrialNum = info.(swpType).nbqx_apv.realTrialNum;
    
    elseif nbqx_apv_Present && nbqx_apv_cd2_ttx_Present
        
        pWidthMatch = info.(swpType).nbqx_apv.pWidth == info.(swpType).nbqx_apv_cd2_ttx.pWidth;
        pAmpMatch = info.(swpType).nbqx_apv.pAmp == info.(swpType).nbqx_apv_cd2_ttx.pAmp;
        if ~(pWidthMatch && pAmpMatch)
            fprintf('>>>   Omiting FV_Na_Ca2_mGluR due to mismatch in pWidth or pAmp\n')
        else
            
            trace.(swpType).FV_Na_Ca2_mGluR = trace.(swpType).nbqx_apv - trace.(swpType).nbqx_apv_cd2_ttx;
            
            info.(swpType).FV_Na_Ca2_mGluR.pulseOn_idx = info.(swpType).nbqx_apv.pulseOn_idx;
            info.(swpType).FV_Na_Ca2_mGluR.pulseOff_idx = info.(swpType).nbqx_apv.pulseOff_idx;
            info.(swpType).FV_Na_Ca2_mGluR.sampRate = info.(swpType).nbqx_apv.sampRate;
            info.(swpType).FV_Na_Ca2_mGluR.pWidth = info.(swpType).nbqx_apv.pWidth;
            info.(swpType).FV_Na_Ca2_mGluR.pAmp = info.(swpType).nbqx_apv.pAmp;
            info.(swpType).FV_Na_Ca2_mGluR.pTF = info.(swpType).nbqx_apv.pTF;
            info.(swpType).FV_Na_Ca2_mGluR.realTrialNum = info.(swpType).nbqx_apv.realTrialNum;
        end
    end
    
    % is there a ttx+cd condition that can be used to define the fiber
    % volley with just Na+
    nbqx_apv_cd2_Present = isfield(trace.(swpType), 'nbqx_apv_cd2');
    if nbqx_apv_cd2_Present && nbqx_apv_cd2_ttx_Present
        
        pWidthMatch = info.(swpType).nbqx_apv_cd2.pWidth == info.(swpType).nbqx_apv_cd2_ttx.pWidth;
        pAmpMatch = info.(swpType).nbqx_apv_cd2.pAmp == info.(swpType).nbqx_apv_cd2_ttx.pAmp;
        assert(pWidthMatch && pAmpMatch, 'ERROR: pulse amp or width are not consistent');
        
        trace.(swpType).FV_Na = trace.(swpType).nbqx_apv_cd2 - trace.(swpType).nbqx_apv_cd2_ttx;
        
        info.(swpType).FV_Na.pulseOn_idx = info.(swpType).nbqx_apv_cd2.pulseOn_idx;
        info.(swpType).FV_Na.pulseOff_idx = info.(swpType).nbqx_apv_cd2.pulseOff_idx;
        info.(swpType).FV_Na.sampRate = info.(swpType).nbqx_apv_cd2.sampRate;
        info.(swpType).FV_Na.pWidth = info.(swpType).nbqx_apv_cd2.pWidth;
        info.(swpType).FV_Na.pAmp = info.(swpType).nbqx_apv_cd2.pAmp;
        info.(swpType).FV_Na.pTF = info.(swpType).nbqx_apv_cd2.pTF;
        info.(swpType).FV_Na.realTrialNum = info.(swpType).nbqx_apv_cd2.realTrialNum;
    end
    
    
    % are there conditions that can be used to assess directRelease due to
    % depolarization of the terminals?
    ttx_only_present = isfield(trace.(swpType), 'ttx');
    cd2_ttx_present = isfield(trace.(swpType), 'cd2_ttx');
    if ttx_only_present && cd2_ttx_present
        
        pWidthMatch = info.(swpType).ttx.pWidth == info.(swpType).cd2_ttx.pWidth;
        pAmpMatch = info.(swpType).ttx.pAmp == info.(swpType).cd2_ttx.pAmp;
        assert(pWidthMatch && pAmpMatch, 'ERROR: pulse amp or width are not consistent');
        
        trace.(swpType).directRelease = trace.(swpType).ttx - trace.(swpType).cd2_ttx;
        
        info.(swpType).directRelease.pulseOn_idx = info.(swpType).ttx.pulseOn_idx;
        info.(swpType).directRelease.pulseOff_idx = info.(swpType).ttx.pulseOff_idx;
        info.(swpType).directRelease.sampRate = info.(swpType).ttx.sampRate;
        info.(swpType).directRelease.pWidth = info.(swpType).ttx.pWidth;
        info.(swpType).directRelease.pAmp = info.(swpType).ttx.pAmp;
        info.(swpType).directRelease.pTF = info.(swpType).ttx.pTF;
        info.(swpType).directRelease.realTrialNum = info.(swpType).ttx.realTrialNum;
    end
    
    
    
    if RMLINENOISE
        % try to reduce line noise from a few of the traces
        conds = {'FV_Na_Ca2_mGluR', 'FV_Na', 'ttx', 'cd2_ttx', 'directRelease', 'synapticTransmission'};
        for i_cond = 1:numel(conds)
            
            if isfield(trace.(swpType), conds{i_cond})
                
                tmp_trace = trace.(swpType).(conds{i_cond});
                
                lines = 60;
                winEnd_idx = size(tmp_trace,1);
                if info.(swpType).(conds{i_cond}).pTF >= 40;
                    % just look at the data following the last pulse.
                    lastpulse = find(info.(swpType).(conds{i_cond}).pulseOff_idx==1, 1, 'last');
                    sampRate = info.(swpType).(conds{i_cond}).sampRate;
                    winStart_idx = lastpulse + ceil(0.15 .* sampRate);
                else
                    % this takes all the data (even the pulses) but is better
                    % then just taking the data after the last pulse b/c long
                    % trains have essentially no data after them.
                    winStart_idx = 1;
                end
                
                for i_ch = 1:size(tmp_trace,2)
                    tmp_trace(:,i_ch) = rmhum(tmp_trace(:,i_ch), sampFreq, winStart_idx, winEnd_idx, lines);
                end
                
                trace.(swpType).(conds{i_cond}) = tmp_trace;
                
            end
        end
    end
    
    % filter some of the traces more agressively
    conds = {'synapticTransmission'};
    for i_cond = 1:numel(conds)
        if isfield(trace.(swpType), conds{i_cond})
            tmp_trace = trace.(swpType).(conds{i_cond});
            tmp_trace = butterfilt(tmp_trace, 1000, sampFreq, 'low', 1);
            trace.(swpType).(conds{i_cond}) = tmp_trace;
        end
    end
    
end





%
% plot the results (one summary figure for each TF)
%
if PLOTFIGURES
    sweepTypeFields = fieldnames(ax);
    for i_sweepType = 1:numel(sweepTypeFields)
        
        % generate the tab labels
        swpType = sweepTypeFields{i_sweepType};
        tabLabels = fieldnames(trace.(swpType));
        
        % one figure for each channel
        for i_ch = 1:sum(channelConfigs);
            
            % create tabbed GUI
            hFig = figure;
            set(gcf, 'position', [40 48 972 711]);
            set(gcf, 'name', sprintf('%s, site %.1f, %s, chan: %d', mouseNames{1}, siteNumber{1}, sweepTypeFields{i_sweepType}, channelList(i_ch)))
            s = warning('off', 'MATLAB:uitabgroup:OldVersion');
            hTabGroup = uitabgroup('Parent',hFig);
            
            for i_cond = 1:numel(tabLabels)
                hTabs(i_cond) = uitab('Parent', hTabGroup, 'Title', tabLabels{i_cond});
                hAx(i_cond) = axes('Parent', hTabs(i_cond));
                hold on,
                sampFreq = info.(swpType).(tabLabels{i_cond}).sampRate;
                tt = (0:size(trace.(swpType).(tabLabels{i_cond}), 1)-1) ./ sampFreq .* 1000;
                firstPulseIdx = find(info.(swpType).(tabLabels{i_cond}).pulseOn_idx == 1, 1, 'first');
                tt = tt-tt(firstPulseIdx);
                plot(tt, trace.(swpType).(tabLabels{i_cond})(:,i_ch), 'k', 'linewidth', 3)
                crossings_on = storedCrossings_on{i_sweepType};
                crossings_off = storedCrossings_off{i_sweepType};
                plot(tt(crossings_on), zeros(1,sum(crossings_on)), 'ro', 'markerfacecolor', 'r')
                plot(tt(crossings_off), zeros(1,sum(crossings_off)), 'mo', 'markerfacecolor', 'r')
                xlabel('time (ms)')
                ylabel('LFP amplitude')
                axis tight
                %xlim([-75, tt(find(crossings_off==1, 1,'last'))+75])
                hold off
                
            end
        end
        
    end
end

