classdef abfobj
    % creates an abf-object.
    
    
    properties
        head = [];
        dat = [];
        wf = [];
        tt = [];
        idx = [];
        name = [];
    end
    
    methods
        function obj = abfobj(fileName, mdb)
            %
            % find and import the data. mdb is an optional argument. If
            % it's not supplied, than the mdb will be generated de novo.
            %
            global GL_DATPATH
            
            if ~exist('fileName', 'var') || isempty(fileName) % no file name supplied
                currentDir = pwd;
                cd(GL_DATPATH)
                [fileName,fpath] = uigetfile('*.abf');
                fpath = [fpath,fileName];
                cd(currentDir);
                
            elseif any(regexpi(fileName, filesep))
                
                fpath = fileName; % the user supplied a fully qualified path
                
            else
                % narrow down the search for findfile.m
                if ~exist('mdb', 'var')
                    suppressVerbose = true;
                    mdb = initMouseDB('update', suppressVerbose);
                end
                mouseNames = mdb.search(fileName(1:10)); %ignore the exact file name and just look at the date
                assert(~isempty(mouseNames), 'ABFOBJ ERROR: could not find data file <%s>', fileName);
                
                % use find file to locate the .abf file
                for a = 1:numel(mouseNames)
                    fpath{a} = findfile(fileName, [GL_DATPATH, mouseNames{a}], '.abf');
                end
                validPaths = cellfun(@(x) ~isempty(x), fpath);
                assert(sum(validPaths) > 0, 'ABFOBJ ERROR: could not find data file <%s>', fileName)
                assert(sum(validPaths) <= 1, 'ABFOBJ ERROR: too many matches for data file <%s>', fileName)
                fpath = fpath{validPaths};
                
            end
            
            % convert the .abf file to a matlab structure. Define a few
            % other things.
            [obj.dat, obj.head, obj.wf] = my_abfload(fpath);
            obj.name = fileName;
            obj.tt = [0:size(obj.dat,1)-1]./obj.head.sampRate;
            
            
            % define the indices into the columns of the aquired signals (obj.dat)
            for a = 1:numel(obj.head.recChNames)
                obj.idx.(obj.head.recChNames{a}) = a;
            end
            
            % define the indices into the columns of the waveforms (obj.wf)
            if isfield(obj.head, 'DACchNames')
                for a = 1:numel(obj.head.DACchNames)
                    fieldname = deblank(obj.head.DACchNames{a}); % sometimes the ouput ch is improperly named....
                    fieldname = fieldname(~isspace(fieldname));
                    if cell2mat(regexpi(fieldname(1), {'\d', '_'})); %fixes leading underscores and numbers
                        fieldname = fieldname(2:end);
                    end
                    obj.idx.(fieldname) = a;
                end
            end
            
        end
        
        function obj = removeSweeps(obj, idx)
            
            %error checking
            assert(size(obj.dat,3)>1, 'ABFOBJ ERROR: only one sweep present, nothing to delete')
            
            % make a logical vector
            newidx = true(1,size(obj.dat,3));
            newidx(idx) = false;
            
            obj.dat = obj.dat(:,:,newidx);
            obj.wf = obj.wf(:,:,newidx);
        end
        
        function out = getvals(obj, ch, sweep, timeStart, timeEnd)
            
            idx = (obj.tt >= timeStart) & (obj.tt < timeEnd);
            out = obj.dat(idx,ch,sweep);
            
        end
        
        function [idx, time] = threshold(obj, thresh, index, direction)
            % index should contain the 2nd, 3rd, 4th, etc dimensions of the
            % sweep to take. The first dimension will be time, and I'll
            % take all time points. For example, if the index is:
            %  wf(:, index(1), index(2))
            % than I'm taking all time points from the "index(1) channel"
            % and the "index(2) sweep"
            
            indexString = '(:';
            for a = 1:numel(index)
                indexString = [indexString, ',', num2str(index(a))];
            end
            indexString = [indexString, ')'];
            
            tmpWF = eval(['obj.wf', indexString]);
            
            above = tmpWF > thresh;
            change = [0; diff(above)];
            
            switch direction
                case 'u'
                    idx = change == 1;
                    time = obj.tt(idx);
                case 'd'
                    idx = change == -1;
                    time = obj.tt(idx);
                otherwise
                    error('ABFOBJ: Threshold crossing direction not specified')
            end
        end
        
        function out = getRa(obj, method)

            % figure out which fitting method to use.
            if ~exist('method', 'var')
                    method = 'quick';
            end
            
            % figure out which channels were turned on, and set to Vclamp
            Vclamp_idx = strcmpi(obj.head.recChUnits, 'pA'); % units are telegraphed from Multiclamp and more accurate than "names"
            secCh = cellfun(@any, regexpi(obj.head.recChNames, '_sec'));
            Vclamp_idx = Vclamp_idx & ~secCh;
            Vclamp_names = obj.head.recChNames(Vclamp_idx);
            Vclamp_names_cmd = obj.head.DACchNames;
            
            % loop over channels and sweeps. Calculate the Ra and Vhold.
            % The return arguments will have the same dimensionality
            % convention as obj.dat (time x channels x sweeps). So the
            % outputs will be [1 x 2 x Nsweeps]
            Nsweeps = size(obj.dat, 3);
            out.Ra = nan(1, 2, Nsweeps);
            out.Verr = nan(1, 2, Nsweeps);
            out.Rinput = nan(1, 2, Nsweeps);
            for i_ch = 1:2
                
                % make sure data are present for this recording channel and
                % initialize the outputs
                HSname = sprintf('HS%d_', i_ch);
                HSpresent = strncmp(Vclamp_names, HSname, numel(HSname));
                cmdname_idx = strncmp(Vclamp_names_cmd, HSname,  numel(HSname));
                assert(sum(HSpresent)<=1, 'ERROR: too many matches');
                if ~any(HSpresent)
                    out.chNames{i_ch} = '';
                    continue
                else
                    HSname = Vclamp_names{HSpresent};  % need to modify in cases where Clampex thinks Iclamp but multiclamp set to Vclamp
                    CMDname = Vclamp_names_cmd{cmdname_idx};
                    out.chNames{i_ch} = HSname;
                end
                
                
                for i_swp = 1:Nsweeps;
                    
                    % find the relevant time points with respect to the
                    % command wf
                    cmdwf = obj.wf(:, obj.idx.(CMDname), i_swp);
                    cmdwf = squeeze(cmdwf);
                    
                    minval = min(cmdwf);
                    assert(abs(minval+5)<1e-6, 'ERROR: may not have found test pulse')
                    thresh = minval .* 0.9;
                    crossing_on = [false; diff(cmdwf<thresh)==1];
                    assert(sum(crossing_on)==1, 'ERROR: too many crossings');
                    idxOnset = find(crossing_on);
                    crossing_off = [false; diff(cmdwf<thresh)==-1];
                    assert(sum(crossing_off)==1, 'ERROR: too many crossings');
                    idxOffset = find(crossing_off);
                    
                    idx_pulse = false(size(cmdwf));
                    idx_pulse(idxOnset:idxOffset) = true;
                    
                    % the Im doesn't start until after the onset of the
                    % Vcmd pulse. Find the most negative point following the
                    % Vcmd pulse
                    minVal = min(obj.dat(idx_pulse, obj.idx.(HSname), i_swp));
                    eq2minval = squeeze(obj.dat(:, obj.idx.(HSname), i_swp)) == minVal;
                    respStart = find([idx_pulse & eq2minval], 1, 'first');
                    t_idx = respStart:idxOffset;
                    Im_pulse = obj.dat(t_idx, obj.idx.(HSname), i_swp);
                    tt_pulse = obj.tt(t_idx);
                    
                    %calculate the baseline
                    N_10ms = round(0.010 .* obj.head.sampRate);
                    idx_baseline = (idxOnset-N_10ms) : (idxOnset-10);
                    Im_baseline = mean(obj.dat(idx_baseline, obj.idx.(HSname), i_swp));
                    
                    %caluclate the magnitude of the pulse
                    Vm_baseline = mean(obj.wf(idx_baseline, obj.idx.(CMDname), i_swp));
                    Vm_pulse = mean(obj.wf(idx_pulse, obj.idx.(CMDname), i_swp));
                    pulse_mv = Vm_pulse - Vm_baseline;
                    
                    switch method
                        case 'quick'
                            delta_pa = minVal - Im_baseline;
                            delta_na = delta_pa ./ 1000;
                            out.Ra(1, i_ch, i_swp) = pulse_mv ./ delta_na;
                            
                        case 'linear'
                            pred = [tt_pulse(1:20)', ones(20,1)];
                            resp = Im_pulse(1:20);
                            betas = pred\resp;
                            tt_On = obj.tt(idxOnset);
                            Im_atOnset = betas(1) .* tt_On + betas(2);
                            delta_pa = Im_atOnset - Im_baseline;
                            delta_na = delta_pa ./ 1000;
                            out.Ra(1, i_ch, i_swp) = pulse_mv ./ delta_na;
                        case 'exp'
                            % currently doesn't do anything
                    end
                    
                    %
                    % calculate the Vclamp error due to holding current as
                    % a percentage of the holding potential
                    %%%%%%%%%%%%%%%%%%%%%%%
                    Verr_volts = (out.Ra(1, i_ch, i_swp) .* 10^6) .* (Im_baseline .* 10^-12);
                    Verr_mv = Verr_volts .* 1000;
                    out.Verr(1, i_ch, i_swp) = abs(Verr_mv);
                    
                    
                    %
                    % calculate the input resistance. Assume a perfect
                    % seal, and understand that the Rin will depend on the
                    % holding potential (which can recruit voltage
                    % dependent currents)
                    %
                    %%%%%%%%%%%%%%%%%%%%%%%
                    Im_pulse = obj.dat(idx_pulse, obj.idx.(HSname), i_swp);
                    Im_pulse(end-5:end) = []; % in case a few samples of the offset made it this far
                    Im_steadyState = mean(Im_pulse(end-N_10ms : end));
                    delta_pa = Im_steadyState - Im_baseline;
                    delta_na = delta_pa ./ 1000;
                    out.Rinput(1, i_ch, i_swp) = pulse_mv ./ delta_na;
                        
                    
                    
                end
            end
        end
        
        function quickPlot(obj)
                        
            % find all secondary channels
            l_secCh = regexp(obj.head.recChNames, '_sec');
            l_secCh = ~cellfun(@isempty, l_secCh);
            
            % find the primary recording from channel 1
            l_ch1 = regexp(obj.head.recChNames, 'HS1_');
            l_ch1 = ~cellfun(@isempty, l_ch1);
            l_ch1 = l_ch1 & ~l_secCh;
            assert(sum(l_ch1)<=1, 'Error: too many primary channels')
            if any(l_ch1)
                raw_ch1 = obj.dat(:,l_ch1,:);
                raw_ch1 = permute(raw_ch1, [1,3,2]);
                
                % grab the WF data
                l_ch1_wf = regexp(obj.head.DACchNames, 'HS1_');
                l_ch1_wf = ~cellfun(@isempty, l_ch1_wf);
                wf_ch1 = obj.wf(:, l_ch1_wf, :);
                wf_ch1 = permute(wf_ch1, [1,3,2]);
            end
            
            
            % find the primary recording from channel 2
            l_ch2 = regexp(obj.head.recChNames, 'HS2_');
            l_ch2 = ~cellfun(@isempty, l_ch2);
            l_ch2 = l_ch2 & ~l_secCh;
            assert(sum(l_ch2)<=1, 'Error: too many primary channels')
            if any(l_ch2)
                raw_ch2 = obj.dat(:,l_ch2,:);
                raw_ch2 = permute(raw_ch2, [1,3,2]);
                
                % grab the WF data
                l_ch2_wf = regexp(obj.head.DACchNames, 'HS2_');
                l_ch2_wf = ~cellfun(@isempty, l_ch2_wf);
                wf_ch2 = obj.wf(:, l_ch2_wf, :);
                wf_ch2 = permute(wf_ch2, [1,3,2]);
            end
            
            % figure out how many plots to make
            f = figure;
            set(f, 'color', [1 1 1])
            hzm = zoom(f);
            set(hzm,'ActionPostCallback', @quickPlot_adjustAxis);
            if any(l_ch1) && any(l_ch2)
                
                set(f, 'position', [239     5   983   801])
                
                subplot(2,2,1)
                quickPlot_rawData(raw_ch1, 1, l_ch1)
                
                subplot(2,2,3)
                quickPlot_commandWF(wf_ch1, 1, l_ch1_wf)
                
                subplot(2,2,2)
                quickPlot_rawData(raw_ch2, 2, l_ch2)
                
                subplot(2,2,4)
                quickPlot_commandWF(wf_ch2, 2, l_ch2_wf)
                
            elseif any(l_ch1)
                
                subplot(2,1,1)
                quickPlot_rawData(raw_ch1, 1, l_ch1)
                
                subplot(2,1,2)
                quickPlot_commandWF(wf_ch1, 1, l_ch1_wf)
            
            elseif any(l_ch2)
                
                subplot(2,1,1)
                quickPlot_rawData(raw_ch2, 2, l_ch2)
                
                subplot(2,1,2)
                quickPlot_commandWF(wf_ch2, 2, l_ch2_wf)
                
            end
            
            
            function quickPlot_rawData(datToPlot, channel, chIdx)
                plot(obj.tt, datToPlot)
                set(gca, 'fontsize', 16)
                xlabel('time');
                ylabel(sprintf('Channel %d (%s)', channel, obj.head.recChUnits{chIdx}))
                xlim([obj.tt(1) obj.tt(end)])
                %ylim([min(datToPlot(:)).*.95 max(datToPlot(:)).*1.05])
                box off
                set(gca, 'userdata', channel);
            end
            
            function quickPlot_commandWF(wfToPlot, channel, chIdx)
                plot(obj.tt, wfToPlot)
                set(gca, 'fontsize', 16)
                xlabel('time')
                ylabel(sprintf('Channel %d (%s)', channel, obj.head.DACchUnits{chIdx}))
                xlim([obj.tt(1) obj.tt(end)])
                ylim([min(wfToPlot(:)).*.95 max(wfToPlot(:)).*1.05])
                box off
                set(gca, 'userdata', channel);
            end
            
            function quickPlot_adjustAxis(~, h)
                
                %grab the new xlims
                newXlims = get(h.Axes, 'xlim');
                
                % figure out which channel was manipulated
                channel = get(h.Axes, 'userdata');
                

                % figure out which axis need to be changed.
                Naxes = numel(get(gcf, 'children'));
                Nchannels = sum(double(l_ch1 | l_ch2));
                Nrows = Naxes ./ Nchannels;
                
                if (Nchannels == 1) && (Nrows == 2);
                    l_axes = [1,2];
                else
                    l_axes = channel:Nchannels:Nchannels+Nrows;
                end
                
                % enforce the new xlims on the WF plot
                for a = 1:numel(l_axes)
                    subplot(Nrows,Nchannels,l_axes(a))
                    set(gca, 'xlim', newXlims)
                end
                
            end
            
            
        end
        
        
    end %methods
    
end