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
                [fileName,fpath] = uigetfile('*.nex');
                fpath = [fpath,fileName];
                cd(currentDir);
            else
                % narrow down the search for findfile.m
                if ~exist('mdb', 'var')
                    mdb = initMouseDB('update', 'notext');
                end
                mouseName = mdb.search(fileName);
                assert(~isempty(mouseName), 'ABFOBJ ERROR: could not find data file <%s>', fileName);
                assert(numel(mouseName)==1, 'ABFOBJ ERROR: too many data files match the input <%s>', fileName);
                
                % use find file to locate the .abf file
                fpath = findfile(fileName, [GL_DATPATH, mouseName{1}], '.abf');
                assert(~isempty(fpath), 'ABFOBJ ERROR: could not locate <%s> in directory <%s>', fileName, mouseName{1})
                
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
        
        function out = fitexp(raw)
            % do this as a cell fun type of thing?
            % Operate separately for cell inputs (cellfun) vs simple
            % inputs?
        end
        
        function out = getvals(obj, ch, sweep, timeStart, timeEnd)
            
            idx = (obj.tt >= timeStart) & (obj.tt < timeEnd);
            out = obj.dat(idx,ch,sweep);
            
        end
        
        function [idx, time] = threshold(obj, thresh, ch, sweep, direction)
            
            above = obj.wf(:,ch,sweep) > thresh;
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
        
        function [Ra, Rin] = getRaRin(obj, method)

            
            [idx_Im, idx_Vclamp] = deal([]);
            for a = 1:numel(obj.head.recChNames);
                
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
            
            % figure out which fitting method to use
            if ~exist('method', 'var')
                lowSampRate = obj.head.sampRate < 40e3;
                slowFilter = any(obj.head.fSignalLowpassFilter(idx_Im)<10e3);
                if lowSampRate || slowFilter
                    method = 'linear';
                else
                    method = 'quick';
                end
            end
            
            % loop over the valid channels and sweeps. Calculate the Ra as
            % you go. The return argument will have the same dimensionality
            % convention as obj.dat (time x channels x sweeps);
            Nsweeps = size(obj.dat, 3);
            Nchannels = numel(idx_Im);
            [Ra, Rin] = deal(nan(1, Nchannels, Nsweeps));
            for ch = 1:numel(idx_Im)
                for swp = 1:Nsweeps;
                    
                    % find the relevant time points with respect to the
                    % command wf
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
                            
%                             % plot everything
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
            end
            
            
            % find the primary recording from channel 2
            l_ch2 = regexp(obj.head.recChNames, 'HS2_');
            l_ch2 = ~cellfun(@isempty, l_ch2);
            l_ch2 = l_ch2 & ~l_secCh;
            assert(sum(l_ch2)<=1, 'Error: too many primary channels')
            if any(l_ch2)
                raw_ch2 = obj.dat(:,l_ch2,:);
                raw_ch2 = permute(raw_ch2, [1,3,2]);
            end
            
            % figure out how many plots to make
            figure            
            if any(l_ch1) && any(l_ch2)
                
                set(gcf, 'position', [239 386 1064 420]);
                
                subplot(1,2,1)
                plot(obj.tt, raw_ch1)
                set(gca, 'fontsize', 16)
                xlabel('time')
                ylabel(sprintf('Channel 1 (%s)', obj.head.recChUnits{l_ch1}))
                xlim([obj.tt(1) obj.tt(end)])
                
                subplot(1,2,2)
                plot(obj.tt, raw_ch2)
                set(gca, 'fontsize', 16)
                xlabel('time')
                ylabel(sprintf('Channel 2 (%s)', obj.head.recChUnits{l_ch2}))
                xlim([obj.tt(1) obj.tt(end)])
                
            elseif any(l_ch1)
                
                plot(obj.tt, raw_ch1)
                set(gca, 'fontsize', 16)
                xlabel('time')
                ylabel(sprintf('Channel 1 (%s)', obj.head.recChUnits{l_ch1}))
                xlim([obj.tt(1) obj.tt(end)])
            
            elseif any(l_ch2)
            
                plot(obj.tt, raw_ch2)
                set(gca, 'fontsize', 16)
                xlabel('time')
                ylabel(sprintf('Channel 2 (%s)', obj.head.recChUnits{l_ch2}))
                xlim([obj.tt(1) obj.tt(end)])
            end
            
        end
        
        
    end %methods
    
end