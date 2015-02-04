function [trace, info] = fiberVolleyAnalysis(exptList, exptWorkbook, PLOTFIGURES)


if ~exist('PLOTFIGURES', 'var')
    PLOTFIGURES = true;
end



% figure out the file names, etc.
mouseNames = exptWorkbook(exptList,1);
siteNumber = exptWorkbook(exptList,2);
fnames = exptWorkbook(exptList,3);
exptConds = exptWorkbook(exptList,5);
channels = exptWorkbook(exptList, 6:7);
rmsweeps = exptWorkbook(exptList,9);
tfs = exptWorkbook(exptList, 4);



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
    fprintf('     Interleaved tfs: %s\n', num2str(fileTFs'));
    
    for i_tf = 1:numel(fileTFs)
        
        % new field names
        field_tf = ['tf_', num2str(fileTFs(i_tf))];
        field_expt = exptConds{i_fid};
        
        % make a franken-abf struct
        ax.(field_tf).(field_expt).head = tmp.head;
        ax.(field_tf).(field_expt).dat = tmp.dat(:,:,tdict.trlList == i_tf);
        ax.(field_tf).(field_expt).idx = tmp.idx;
        ax.(field_tf).(field_expt).pAmp = tdict.conds(i_tf,1);
        ax.(field_tf).(field_expt).pWidth = tdict.conds(i_tf,2);
    end

end


% pull out the raw data, filter, and organize all the raw traces.
fprintf('Filtering and organizing raw data\n');
trace = [];
TFfields = fieldnames(ax);
for i_tf = 1:numel(TFfields)
    
    conditionFields = fieldnames(ax.(TFfields{i_tf}));
    
    for i_cond = 1:numel(conditionFields)
        
        % generate some field names for extraction and saving
        field_tf = TFfields{i_tf};
        field_expt = conditionFields{i_cond};
        
        sampFreq = ax.(field_tf).(field_expt).head.sampRate;
        %tt = (0:size(ax.(field_tf).(field_expt).dat, 1)-1) ./ sampFreq .* 1000;
        
        
        % grab the raw data (channel by channel)
        validChans = ax.(field_tf).(field_expt).head.validChans;
        validChans = find(validChans == 1);
        for i_ch = 1:numel(validChans);
            
            % figure out the appropriate index to the data. It should have
            % the correct "HSx_" prefix, and have the correct units
            whichChan = validChans(i_ch);
            correctUnits = strcmpi('mV', ax.(field_tf).(field_expt).head.recChUnits);
            correctHS = strncmpi(sprintf('HS%d_', whichChan), ax.(field_tf).(field_expt).head.recChNames, 3);
            chIdx = correctUnits & correctHS;
            assert(sum(chIdx)==1, 'ERROR: incorrect channel selection')
            
            tmp = ax.(field_tf).(field_expt).dat(:, chIdx,:);
            tmp = squeeze(tmp);
            
            
            % baseline subtract, determine when the pulses came on...
            thresh = 0.05;
            ledIdx = ax.(field_tf).(field_expt).idx.LED_470;
            template = ax.(field_tf).(field_expt).dat(:,ledIdx,1);
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
            trace.(field_tf).(field_expt)(:,i_ch) = average;
            
        end % i_ch
        
        % store the pulse onset times for the population analysis
        info.(field_tf).(field_expt).pulseOn_idx = storedCrossings_on{i_tf};
        info.(field_tf).(field_expt).pulseOff_idx = storedCrossings_off{i_tf};
        info.(field_tf).(field_expt).sampRate = ax.(field_tf).(field_expt).head.sampRate;
        info.(field_tf).(field_expt).pWidth = ax.(field_tf).(field_expt).pWidth;
        info.(field_tf).(field_expt).pAmp = ax.(field_tf).(field_expt).pAmp;
        
    end % i_cond
end % i_tf



fprintf('Calculating fiber volley\n')
TFfields = fieldnames(ax);
for i_tf = 1:numel(TFfields)
    
    % rules:
    %
    % nbqx_apv - nbqx_apv_ttx           =>   fiber volley with Na and Ca2 and presynpatic mGluR
    % nbqx_apv - nbqx_apv_cd_ttx        =>   fiber volley with Na and Ca2 and presynpatic mGluR
    % nbqx_apv_cd2 - nbqx_apv_cd2_ttx   =>   fiber volley with just Na
    % none - nbqx_apv                   =>   synapticTransmission
    % none                              =>   control
    
    
    % generate some field names
    field_tf = TFfields{i_tf};
    
    % is there control data and nbqx_apv?
    control_Present = isfield(trace.(field_tf), 'none');
    nbqx_apv_Present = isfield(trace.(field_tf), 'nbqx_apv');
    if control_Present && nbqx_apv_Present
        
        pWidthMatch = info.(field_tf).none.pWidth == info.(field_tf).nbqx_apv.pWidth;
        pAmpMatch = info.(field_tf).none.pAmp == info.(field_tf).nbqx_apv.pAmp;
        assert(pWidthMatch && pAmpMatch, 'ERROR: pulse amp or width are not consistent');
        
        trace.(field_tf).synapticTransmission = trace.(field_tf).none - trace.(field_tf).nbqx_apv;
        
        info.(field_tf).synapticTransmission.pulseOn_idx = info.(field_tf).none.pulseOn_idx;
        info.(field_tf).synapticTransmission.pulseOff_idx = info.(field_tf).none.pulseOff_idx;
        info.(field_tf).synapticTransmission.sampRate = info.(field_tf).none.sampRate;
        info.(field_tf).synapticTransmission.pWidth = info.(field_tf).none.pWidth;
        info.(field_tf).synapticTransmission.pAmp = info.(field_tf).none.pAmp;
    end
    
    % is there a ttx condition that can be used to define fiber volley with
    % sodium and calcium currents?
    nbqx_apv_ttx_Present = isfield(trace.(field_tf), 'nbqx_apv_ttx');
    nbqx_apv_cd2_ttx_Present = isfield(trace.(field_tf), 'nbqx_apv_cd2_ttx');
    assert(~all([nbqx_apv_ttx_Present, nbqx_apv_cd2_ttx_Present]), 'ERROR: multiple TTX files defined');
    if nbqx_apv_Present && nbqx_apv_ttx_Present
        
        pWidthMatch = info.(field_tf).nbqx_apv.pWidth == info.(field_tf).nbqx_apv_ttx.pWidth;
        pAmpMatch = info.(field_tf).nbqx_apv.pAmp == info.(field_tf).nbqx_apv_ttx.pAmp;
        assert(pWidthMatch && pAmpMatch, 'ERROR: pulse amp or width are not consistent');
        
        trace.(field_tf).FV_Na_Ca2_mGluR = trace.(field_tf).nbqx_apv - trace.(field_tf).nbqx_apv_ttx;
        
        info.(field_tf).FV_Na_Ca2_mGluR.pulseOn_idx = info.(field_tf).nbqx_apv.pulseOn_idx;
        info.(field_tf).FV_Na_Ca2_mGluR.pulseOff_idx = info.(field_tf).nbqx_apv.pulseOff_idx;
        info.(field_tf).FV_Na_Ca2_mGluR.sampRate = info.(field_tf).nbqx_apv.sampRate;
        info.(field_tf).FV_Na_Ca2_mGluR.pWidth = info.(field_tf).nbqx_apv.pWidth;
        info.(field_tf).FV_Na_Ca2_mGluR.pAmp = info.(field_tf).nbqx_apv.pAmp;
    
    elseif nbqx_apv_Present && nbqx_apv_cd2_ttx_Present
        
        pWidthMatch = info.(field_tf).nbqx_apv.pWidth == info.(field_tf).nbqx_apv_cd2_ttx.pWidth;
        pAmpMatch = info.(field_tf).nbqx_apv.pAmp == info.(field_tf).nbqx_apv_cd2_ttx.pAmp;
        assert(pWidthMatch && pAmpMatch, 'ERROR: pulse amp or width are not consistent');
        
        trace.(field_tf).FV_Na_Ca2_mGluR = trace.(field_tf).nbqx_apv - trace.(field_tf).nbqx_apv_cd2_ttx;
        
        info.(field_tf).FV_Na_Ca2_mGluR.pulseOn_idx = info.(field_tf).nbqx_apv.pulseOn_idx;
        info.(field_tf).FV_Na_Ca2_mGluR.pulseOff_idx = info.(field_tf).nbqx_apv.pulseOff_idx;
        info.(field_tf).FV_Na_Ca2_mGluR.sampRate = info.(field_tf).nbqx_apv.sampRate;
        info.(field_tf).FV_Na_Ca2_mGluR.pWidth = info.(field_tf).nbqx_apv.pWidth;
        info.(field_tf).FV_Na_Ca2_mGluR.pAmp = info.(field_tf).nbqx_apv.pAmp;
    end
    
    % is there a ttx+cd condition that can be used to define the fiber
    % volley with just Na+
    nbqx_apv_cd2_Present = isfield(trace.(field_tf), 'nbqx_apv_cd2');
    if nbqx_apv_cd2_Present && nbqx_apv_cd2_ttx_Present
        
        pWidthMatch = info.(field_tf).nbqx_apv_cd2.pWidth == info.(field_tf).nbqx_apv_cd2_ttx.pWidth;
        pAmpMatch = info.(field_tf).nbqx_apv_cd2.pAmp == info.(field_tf).nbqx_apv_cd2_ttx.pAmp;
        assert(pWidthMatch && pAmpMatch, 'ERROR: pulse amp or width are not consistent');
        
        trace.(field_tf).FV_Na = trace.(field_tf).nbqx_apv_cd2 - trace.(field_tf).nbqx_apv_cd2_ttx;
        
        info.(field_tf).FV_Na.pulseOn_idx = info.(field_tf).nbqx_apv_cd2.pulseOn_idx;
        info.(field_tf).FV_Na.pulseOff_idx = info.(field_tf).nbqx_apv_cd2.pulseOff_idx;
        info.(field_tf).FV_Na.sampRate = info.(field_tf).nbqx_apv_cd2.sampRate;
        info.(field_tf).FV_Na.pWidth = info.(field_tf).nbqx_apv_cd2.pWidth;
        info.(field_tf).FV_Na.pAmp = info.(field_tf).nbqx_apv_cd2.pAmp;
    end
    
    
    % try to reduce line noise from a few of the traces
    conds = {'FV_Na_Ca2_mGluR', 'FV_Na'};
    for i_cond = 1:numel(conds)
        
        if isfield(trace.(field_tf), conds{i_cond})
            lines = [5.5, 60, 120, 180];
            winStart_idx = 1; % this takes all the data (even the pulses) but is better then just taking the data after the last pulse b/c long trains have essentially no data after them. 
            winEnd_idx = numel(storedCrossings_off{i_tf});
            tmp_trace = trace.(field_tf).(conds{i_cond});
            
            for i_ch = 1:size(tmp_trace,2)
                tmp_trace(:,i_ch) = rmhum(tmp_trace(:,i_ch), sampFreq, winStart_idx, winEnd_idx, lines);
            end
            
            trace.(field_tf).(conds{i_cond}) = tmp_trace;
            
        end
    end
        
    
end





%
% plot the results (one summary figure for each TF)
%
if PLOTFIGURES
    for i_tf = 1:numel(TFfields)
        
        % generate the tab labels
        field_tf = TFfields{i_tf};
        tabLabels = fieldnames(trace.(field_tf));
        
        % one figure for each channel
        for i_ch = 1:sum(channelConfigs);
            
            % create tabbed GUI
            hFig = figure;
            set(gcf, 'position', [40 48 972 711]);
            set(gcf, 'name', sprintf('%s, site %d, %s, chan: %d', mouseNames{1}, siteNumber{1}, TFfields{i_tf}, channelList(i_ch)))
            s = warning('off', 'MATLAB:uitabgroup:OldVersion');
            hTabGroup = uitabgroup('Parent',hFig);
            
            for i_cond = 1:numel(tabLabels)
                
                hTabs(i_cond) = uitab('Parent', hTabGroup, 'Title', tabLabels{i_cond});
                hAx(i_cond) = axes('Parent', hTabs(i_cond));
                hold on,
                sampFreq = info.(field_tf).(tabLabels{i_cond}).sampRate;
                tt = (0:size(ax.(field_tf).(field_expt).dat, 1)-1) ./ sampFreq .* 1000;
                firstPulseIdx = find(info.(field_tf).(tabLabels{i_cond}).pulseOn_idx == 1, 1, 'first');
                tt = tt-tt(firstPulseIdx);
                plot(tt, trace.(field_tf).(tabLabels{i_cond})(:,i_ch), 'k', 'linewidth', 3)
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

