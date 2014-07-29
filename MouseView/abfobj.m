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
            else
                % narrow down the search for findfile.m
                if ~exist('mdb', 'var')
                    mdb = initMouseDB(false, true);
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
                    obj.idx.(obj.head.DACchNames{a}) = a;
                end
            end
            
        end
        
        function obj = removeSweeps(obj, idx)
            
            assert(size(obj.dat,3)>1, 'ABFOBJ ERROR: only one sweep present, nothing to delete')
            obj.dat(:,:,idx) = [];
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
        
        function Ra = getRa(obj, method)

            
            [idx_Im, idx_Vclamp] = deal([]);
            for a = 1:numel(obj.head.recChNames);
                
                % only consider the primary channel (for example, the
                % secondary channel will record pA when in Current Clamp)
                secCh = regexp(obj.head.recChNames{a}, '_sec', 'match');
                if ~isempty(secCh); continue; end
                
                
                % bail if the channel is non-neural
                validChName = regexp(obj.head.recChNames{a}, 'HS\d_', 'match');                
                if isempty(validChName); continue; end
                
                % bail if the channel is not Vclamp
                vclamp = regexpi(obj.head.recChUnits{a}, 'pA');
                if isempty(vclamp); continue; end
                
                % the index to the recorded membrane current
                tmp = eval(['obj.idx.', validChName{1}, 'Im']);
                idx_Im = [idx_Im, tmp];
                
                % the index to the comand waveform
                tmp = eval(['obj.idx.', validChName{1}, 'Vclamp']);
                idx_Vclamp = [idx_Vclamp, tmp];              
                
            end
            
            % figure out which fitting method to use. Once exponential fits
            % are up and running, than I could institute exp fits for cases
            % where the sampling rate was slow, or there was a filter that
            % slows the kinetics. Until then, just make the most negative
            % value the estimate. This likely overestimates Ra, but that's
            % not such a bad thing...
            if ~exist('method', 'var')
                    method = 'quick';
            end
            
            % loop over the valid channels and sweeps. Calculate the Ra as
            % you go. The return argument will have the same dimensionality
            % convention as obj.dat (time x channels x sweeps);
            Nsweeps = size(obj.dat, 3);
            Nchannels = numel(idx_Im);
            Ra = nan(1, Nchannels, Nsweeps);
            for ch = 1:numel(idx_Im)
                for swp = 1:Nsweeps;
                    
                    % find the relevant time points with respect to the
                    % command wf
                    keyboard
                    idxOnset = obj.threshold(-0.1, idx_Vclamp(ch), swp, 'd');
                    idxOffset = obj.threshold(-0.1, idx_Vclamp(ch), swp, 'u');
                    idx_pulse = find(idxOnset) : find(idxOffset);
                    
                    % the Im doesn't start until after the onset of the
                    % Vcmd pulse. Find the most negative point folloing the
                    % Vcmd pulse
                    minVal = min(obj.dat(idx_pulse, idx_Im(ch), swp));
                    respStart = find(obj.dat(:, idx_Im(ch), swp) == minVal, 1, 'first');
                    t_idx = respStart:find(idxOffset);
                    Im_pulse = obj.dat(t_idx, idx_Im(ch), swp);
                    tt_pulse = obj.tt(t_idx);
                    
                    %calculate the baseline
                    idx_baseline = find(idxOnset)-110 : find(idxOnset)-10;
                    Im_baseline = mean(obj.dat(idx_baseline, idx_Im(ch), swp));
                    
                    %caluclate the magnitude of the pulse
                    Vm_baseline = mean(obj.wf(idx_baseline, idx_Vclamp(ch), swp));
                    Vm_pulse = mean(obj.wf(idx_pulse, idx_Vclamp(ch), swp));
                    pulse_mv = Vm_pulse - Vm_baseline;
                    
                    switch method
                        case 'quick'
                            delta_pa = minVal - Im_baseline;
                            delta_na = delta_pa ./ 1000;
                            Ra(1,idx_Im(ch), swp) = pulse_mv ./ delta_na;
                            
                        case 'linear'
                            pred = [tt_pulse(1:20)', ones(20,1)];
                            resp = Im_pulse(1:20);
                            betas = pred\resp;
                            tt_On = obj.tt(idxOnset);
                            Im_atOnset = betas(1) .* tt_On + betas(2);
                            delta_pa = Im_atOnset - Im_baseline;
                            delta_na = delta_pa ./ 1000;
                            Ra(1,idx_Im(ch), swp) = pulse_mv ./ delta_na;
                            
                            % plot everything
%                             figure,
%                             subplot(2,1,1), hold on,
%                             plot(obj.tt, obj.dat(:, idx_Im(ch), swp), 'k')
%                             plot(tt_pulse, Im_pulse, 'ro')
%                             plot(obj.tt(idx_baseline), obj.dat(idx_baseline, idx_Im(ch), swp), 'co')
%                             plot([pred(:,1);tt_On], betas(1).*[pred(:,1);tt_On]+betas(2), 'b', 'linewidth', 2)
%                             xlim([0.139 0.1394])
%                             subplot(2,1,2)
%                             plot(obj.tt, obj.wf(:, idx_Vclamp(ch), swp), 'k')
%                             xlim([0.139 0.1394])
                            
                        case 'exp'
                            % currently doesn't do anything
                    end
                    
                    
                    
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
                box off
                set(gca, 'userdata', channel);
            end
            
            function quickPlot_commandWF(wfToPlot, channel, chIdx)
                plot(obj.tt, wfToPlot)
                set(gca, 'fontsize', 16)
                xlabel('time')
                ylabel(sprintf('Channel %d (%s)', channel, obj.head.DACchUnits{chIdx}))
                xlim([obj.tt(1) obj.tt(end)])
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