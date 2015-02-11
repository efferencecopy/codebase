function [trace, info] = fiberVolleyAnalysis(exptList, exptWorkbook, PLOTFIGURES)


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
tfs = exptWorkbook(exptList, strcmpi(exptWorkbook(1,:), 'TF'));



% ERROR CHECKING: check to make sure that the exptConds are correct
validConds = {'none',...
              'cd2',...
              'nbqx_apv',...
              'nbqx_apv_cd2',...
              'nbqx_apv_cd2_ttx',...
              'nbqx_apv_ttx'};
exptConds = cellfun(@(x) regexprep(x, '\+', '_'), exptConds, 'uniformoutput', false);
idx = cellfun(@(x,y) regexpi(x,y), exptConds, repmat({validConds}, size(exptConds)), 'uniformoutput', false);
idx = cellfun(@(x) any([x{:}]), idx);
assert(all(idx), 'ERROR: at least one experimental condition is not recognized');

% ERROR CHECKING: make sure that there tmp.head.validChans is consistency
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
    
    % store the data, grouped by TF. First, I need to figure out the
    % appropriate field name for each TF instance
    
    tdict = outerleave(tmp, tmp.idx.LED_470);
    
    fileTFs = round(tdict.conds(:,3));
    fileAmps = round(tdict.conds(:,1));
    fprintf('     Interleaved tfs: %s\n     Interleaved amps: %s\n', num2str(fileTFs'), num2str(fileAmps'));
    
    for i_tf = 1:numel(fileTFs)
        
        % new field names
        tType = ['tf', num2str(fileTFs(i_tf)), '_amp', num2str(fileAmps(i_tf))];
        field_expt = exptConds{i_fid};
        
        % make a franken-abf struct
        ax.(tType).(field_expt).head = tmp.head;
        ax.(tType).(field_expt).dat = tmp.dat(:,:,tdict.trlList == i_tf);
        ax.(tType).(field_expt).idx = tmp.idx;
        ax.(tType).(field_expt).pAmp = tdict.conds(i_tf,1);
        ax.(tType).(field_expt).pWidth = tdict.conds(i_tf,2);
        ax.(tType).(field_expt).pTF = fileTFs(i_tf);
    end

end


% pull out the raw data, filter, and organize all the raw traces.
fprintf('Filtering and organizing raw data\n');
trace = [];
structFields = fieldnames(ax);
for i_tf = 1:numel(structFields)
    
    conditionFields = fieldnames(ax.(structFields{i_tf}));
    
    for i_cond = 1:numel(conditionFields)
        
        % generate some field names for extraction and saving
        tType = structFields{i_tf};
        field_expt = conditionFields{i_cond};
        
        sampFreq = ax.(tType).(field_expt).head.sampRate;
        %tt = (0:size(ax.(tType).(field_expt).dat, 1)-1) ./ sampFreq .* 1000;
        
        
        % grab the raw data (channel by channel)
        validChans = ax.(tType).(field_expt).head.validChans;
        validChans = find(validChans == 1);
        for i_ch = 1:numel(validChans);
            
            % figure out the appropriate index to the data. It should have
            % the correct "HSx_" prefix, and have the correct units
            whichChan = validChans(i_ch);
            correctUnits = strcmpi('mV', ax.(tType).(field_expt).head.recChUnits);
            correctHS = strncmpi(sprintf('HS%d_', whichChan), ax.(tType).(field_expt).head.recChNames, 3);
            chIdx = correctUnits & correctHS;
            assert(sum(chIdx)==1, 'ERROR: incorrect channel selection')
            
            tmp = ax.(tType).(field_expt).dat(:, chIdx,:);
            tmp = squeeze(tmp);
            
            
            % baseline subtract, determine when the pulses came on...
            thresh = 0.05;
            ledIdx = ax.(tType).(field_expt).idx.LED_470;
            template = ax.(tType).(field_expt).dat(:,ledIdx,1);
            template = template > thresh;
            template = [0; diff(template)];
            storedCrossings_on{i_tf} = template == 1;
            tmp_off = template == -1;
            storedCrossings_off{i_tf} = [tmp_off(2:end); false]; % i think there is an OBO error otherwise...
            pulseOnset = find(template==1, 1, 'first');
            
            bkgndTime = 0.150; 
            bkgndSamps = ceil(bkgndTime .* sampFreq);
            tmp = bsxfun(@minus, tmp, mean(tmp(pulseOnset-bkgndSamps:pulseOnset-1, :),1));
            
            
            % filter out the high frequency stuff. Filtering shouldn't go
            % below 2500 b/c you'll start to carve out the fiber volley
            % (which is only a few ms wide)
            lp_freq = 2500;
            filtered = butterfilt(tmp, lp_freq, sampFreq, 'low', 1);
            
            % take the mean
            average = mean(filtered,2);
            trace.(tType).(field_expt)(:,i_ch) = average;
            
        end % i_ch
        
        % store the pulse onset times for the population analysis
        info.(tType).(field_expt).pulseOn_idx = storedCrossings_on{i_tf};
        info.(tType).(field_expt).pulseOff_idx = storedCrossings_off{i_tf};
        info.(tType).(field_expt).sampRate = ax.(tType).(field_expt).head.sampRate;
        info.(tType).(field_expt).pWidth = ax.(tType).(field_expt).pWidth;
        info.(tType).(field_expt).pAmp = ax.(tType).(field_expt).pAmp;
        info.(tType).(field_expt).pTF = ax.(tType).(field_expt).pTF;
        
    end % i_cond
end % i_tf



fprintf('Calculating fiber volley\n')
structFields = fieldnames(ax);
for i_tf = 1:numel(structFields)
    
    % rules:
    %
    % nbqx_apv - nbqx_apv_ttx           =>   fiber volley with Na and Ca2 and presynpatic mGluR
    % nbqx_apv - nbqx_apv_cd_ttx        =>   fiber volley with Na and Ca2 and presynpatic mGluR
    % nbqx_apv_cd2 - nbqx_apv_cd2_ttx   =>   fiber volley with just Na
    % none - nbqx_apv                   =>   synapticTransmission
    % none                              =>   control
    
    
    % generate some field names
    tType = structFields{i_tf};
    
    % is there control data and nbqx_apv?
    control_Present = isfield(trace.(tType), 'none');
    nbqx_apv_Present = isfield(trace.(tType), 'nbqx_apv');
    if control_Present && nbqx_apv_Present
        
        pWidthMatch = info.(tType).none.pWidth == info.(tType).nbqx_apv.pWidth;
        pAmpMatch = info.(tType).none.pAmp == info.(tType).nbqx_apv.pAmp;
        assert(pWidthMatch && pAmpMatch, 'ERROR: pulse amp or width are not consistent');
        
        trace.(tType).synapticTransmission = trace.(tType).none - trace.(tType).nbqx_apv;
        
        info.(tType).synapticTransmission.pulseOn_idx = info.(tType).none.pulseOn_idx;
        info.(tType).synapticTransmission.pulseOff_idx = info.(tType).none.pulseOff_idx;
        info.(tType).synapticTransmission.sampRate = info.(tType).none.sampRate;
        info.(tType).synapticTransmission.pWidth = info.(tType).none.pWidth;
        info.(tType).synapticTransmission.pAmp = info.(tType).none.pAmp;
        info.(tType).synapticTransmission.pTF = info.(tType).none.pTF;
    end
    
    % is there a ttx condition that can be used to define fiber volley with
    % sodium and calcium currents?
    nbqx_apv_ttx_Present = isfield(trace.(tType), 'nbqx_apv_ttx');
    nbqx_apv_cd2_ttx_Present = isfield(trace.(tType), 'nbqx_apv_cd2_ttx');
    assert(~all([nbqx_apv_ttx_Present, nbqx_apv_cd2_ttx_Present]), 'ERROR: multiple TTX files defined');
    if nbqx_apv_Present && nbqx_apv_ttx_Present
        
        pWidthMatch = info.(tType).nbqx_apv.pWidth == info.(tType).nbqx_apv_ttx.pWidth;
        pAmpMatch = info.(tType).nbqx_apv.pAmp == info.(tType).nbqx_apv_ttx.pAmp;
        assert(pWidthMatch && pAmpMatch, 'ERROR: pulse amp or width are not consistent');
        
        trace.(tType).FV_Na_Ca2_mGluR = trace.(tType).nbqx_apv - trace.(tType).nbqx_apv_ttx;
        
        info.(tType).FV_Na_Ca2_mGluR.pulseOn_idx = info.(tType).nbqx_apv.pulseOn_idx;
        info.(tType).FV_Na_Ca2_mGluR.pulseOff_idx = info.(tType).nbqx_apv.pulseOff_idx;
        info.(tType).FV_Na_Ca2_mGluR.sampRate = info.(tType).nbqx_apv.sampRate;
        info.(tType).FV_Na_Ca2_mGluR.pWidth = info.(tType).nbqx_apv.pWidth;
        info.(tType).FV_Na_Ca2_mGluR.pAmp = info.(tType).nbqx_apv.pAmp;
        info.(tType).FV_Na_Ca2_mGluR.pTF = info.(tType).nbqx_apv.pTF;
    
    elseif nbqx_apv_Present && nbqx_apv_cd2_ttx_Present
        
        pWidthMatch = info.(tType).nbqx_apv.pWidth == info.(tType).nbqx_apv_cd2_ttx.pWidth;
        pAmpMatch = info.(tType).nbqx_apv.pAmp == info.(tType).nbqx_apv_cd2_ttx.pAmp;
        assert(pWidthMatch && pAmpMatch, 'ERROR: pulse amp or width are not consistent');
        
        trace.(tType).FV_Na_Ca2_mGluR = trace.(tType).nbqx_apv - trace.(tType).nbqx_apv_cd2_ttx;
        
        info.(tType).FV_Na_Ca2_mGluR.pulseOn_idx = info.(tType).nbqx_apv.pulseOn_idx;
        info.(tType).FV_Na_Ca2_mGluR.pulseOff_idx = info.(tType).nbqx_apv.pulseOff_idx;
        info.(tType).FV_Na_Ca2_mGluR.sampRate = info.(tType).nbqx_apv.sampRate;
        info.(tType).FV_Na_Ca2_mGluR.pWidth = info.(tType).nbqx_apv.pWidth;
        info.(tType).FV_Na_Ca2_mGluR.pAmp = info.(tType).nbqx_apv.pAmp;
        info.(tType).FV_Na_Ca2_mGluR.pTF = info.(tType).nbqx_apv.pTF;
    end
    
    % is there a ttx+cd condition that can be used to define the fiber
    % volley with just Na+
    nbqx_apv_cd2_Present = isfield(trace.(tType), 'nbqx_apv_cd2');
    if nbqx_apv_cd2_Present && nbqx_apv_cd2_ttx_Present
        
        pWidthMatch = info.(tType).nbqx_apv_cd2.pWidth == info.(tType).nbqx_apv_cd2_ttx.pWidth;
        pAmpMatch = info.(tType).nbqx_apv_cd2.pAmp == info.(tType).nbqx_apv_cd2_ttx.pAmp;
        assert(pWidthMatch && pAmpMatch, 'ERROR: pulse amp or width are not consistent');
        
        trace.(tType).FV_Na = trace.(tType).nbqx_apv_cd2 - trace.(tType).nbqx_apv_cd2_ttx;
        
        info.(tType).FV_Na.pulseOn_idx = info.(tType).nbqx_apv_cd2.pulseOn_idx;
        info.(tType).FV_Na.pulseOff_idx = info.(tType).nbqx_apv_cd2.pulseOff_idx;
        info.(tType).FV_Na.sampRate = info.(tType).nbqx_apv_cd2.sampRate;
        info.(tType).FV_Na.pWidth = info.(tType).nbqx_apv_cd2.pWidth;
        info.(tType).FV_Na.pAmp = info.(tType).nbqx_apv_cd2.pAmp;
        info.(tType).FV_Na.pTF = info.(tType).nbqx_apv_cd2.pTF;
    end
    
    
    % try to reduce line noise from a few of the traces
    conds = {'FV_Na_Ca2_mGluR', 'FV_Na'};
    for i_cond = 1:numel(conds)
        
        if isfield(trace.(tType), conds{i_cond})
            lines = [5.5, 60, 120, 180];
            winStart_idx = 1; % this takes all the data (even the pulses) but is better then just taking the data after the last pulse b/c long trains have essentially no data after them. 
            winEnd_idx = numel(storedCrossings_off{i_tf});
            tmp_trace = trace.(tType).(conds{i_cond});
            
            for i_ch = 1:size(tmp_trace,2)
                warning('RMHUM not in the loop')
                %tmp_trace(:,i_ch) = rmhum(tmp_trace(:,i_ch), sampFreq, winStart_idx, winEnd_idx, lines);
            end
            
            trace.(tType).(conds{i_cond}) = tmp_trace;
            
        end
    end
        
    
end





%
% plot the results (one summary figure for each TF)
%
if PLOTFIGURES
    for i_tf = 1:numel(structFields)
        
        % generate the tab labels
        tType = structFields{i_tf};
        tabLabels = fieldnames(trace.(tType));
        
        % one figure for each channel
        for i_ch = 1:sum(channelConfigs);
            
            % create tabbed GUI
            hFig = figure;
            set(gcf, 'position', [40 48 972 711]);
            set(gcf, 'name', sprintf('%s, site %d, %s, chan: %d', mouseNames{1}, siteNumber{1}, structFields{i_tf}, channelList(i_ch)))
            s = warning('off', 'MATLAB:uitabgroup:OldVersion');
            hTabGroup = uitabgroup('Parent',hFig);
            
            for i_cond = 1:numel(tabLabels)
                
                hTabs(i_cond) = uitab('Parent', hTabGroup, 'Title', tabLabels{i_cond});
                hAx(i_cond) = axes('Parent', hTabs(i_cond));
                hold on,
                sampFreq = info.(tType).(tabLabels{i_cond}).sampRate;
                tt = (0:size(ax.(tType).(field_expt).dat, 1)-1) ./ sampFreq .* 1000;
                firstPulseIdx = find(info.(tType).(tabLabels{i_cond}).pulseOn_idx == 1, 1, 'first');
                tt = tt-tt(firstPulseIdx);
                plot(tt, trace.(tType).(tabLabels{i_cond})(:,i_ch), 'k', 'linewidth', 3)
                crossings_on = storedCrossings_on{i_tf};
                crossings_off = storedCrossings_off{i_tf};
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

