function params = anlyMod_avgOuterleave(params)

    %
    % CALCULATE THE AVERAGE RESPONSE
    %
    % Iterate over each experimental group, and then over each data file
    % within group.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % construct a library of the available data files
    abfLibrary = {};
    for a = 1:numel(params.ax)
        abfLibrary{a} = params.ax{a}.name;
    end

    
    nExptGroups = size(params.groups,1);
    for i_group = 1:nExptGroups;
        
        groupName = params.groups{i_group,1};
        
        % determine the files that contribute to this group, and their
        % index into params.ax
        prefix = params.files{1};
        list = params.groups{i_group,2};
        groupFiles = {};
        for i = 1:numel(list);
            suffix = num2str(list(i));
            nZerosNeeded = 4-numel(suffix);
            suffix = [repmat('0',1,nZerosNeeded), suffix];
            groupFiles{i} = [prefix, suffix];
        end
        
        
        for i_ax = 1:numel(groupFiles)
            
            axIdx = cellfun(@(x) ~isempty(x), regexpi(groupFiles{i_ax}, abfLibrary));
            
            ledIdx = params.ax{axIdx}.idx.LED_470;
            tmpWF = params.ax{axIdx}.dat(:, ledIdx,:);
            sampFreq = params.ax{axIdx}.head.sampRate;
            tdict = outerleave(tmpWF, sampFreq);
            params.(groupName).tdict{i_ax} = tdict;
            params.(groupName).head{i_ax} = params.ax{axIdx}.head;
            
            
            % calculate the pulse times, which will be used later
            tmpWF = params.ax{axIdx}.dat(:,ledIdx,:);
            tmpWF = permute(tmpWF, [1,3,2]);
            [Nsamps, Nsweeps] = size(tmpWF);
            thresh = 0.03; % should greatly exceed the noise in the sampled LED output
            above = tmpWF > thresh;
            change = [zeros(1,Nsweeps); diff(above, 1, 1)];
            threshCrossing = change == 1;
            
            
            % figure out how many channels there were
            l_secCh = cellfun(@(x) ~isempty(x), regexpi(params.ax{axIdx}.head.recChNames, 'sec'));
            l_hs1 = cellfun(@(x) ~isempty(x), regexpi(params.ax{axIdx}.head.recChNames, 'hs1'));
            l_hs2 = cellfun(@(x) ~isempty(x), regexpi(params.ax{axIdx}.head.recChNames, 'hs2'));
            primaryChIdx = [find(l_hs1 & ~l_secCh), find(l_hs2 & ~l_secCh)];
            secondaryChIdx = [find(l_hs1 & l_secCh), find(l_hs2 & l_secCh)];
            assert(numel(primaryChIdx)<=2, 'ERROR: too many channels')
            
            for i_ch = 1:numel(primaryChIdx);
                
                chIdx = primaryChIdx(i_ch);
                
                % generate a trial list for this specific channel and .abf file
                % based off the "exclusion list" defined in physiology_notes.m
                l_goodSweeps = findGoodSweeps(params, axIdx, i_ch);
                if strcmpi(l_goodSweeps, 'exclude_file')
                    continue
                end
                
                % loop over the conditions and calculate the mean current
                for i_cond = 1:size(tdict.conds,1);
                    
                    % trial list based off tdict.trlList, then combine the two
                    % lists into a single list
                    l_cond = tdict.trlList == i_cond;
                    l_valid = l_goodSweeps & l_cond;
                    
                    if sum(l_valid) == 0
                        params.(groupName).avg.trace_pA{i_ax}{i_cond, i_ch} = nan(Nsamps, 1);
                        continue
                    end
                    
                    % grab the raw data, mean subtract, and then average across sweeps
                    tmp_raw = params.ax{axIdx}.dat(:,chIdx,l_valid);
                    tmp_raw = permute(tmp_raw, [1,3,2]);
                    dim = 1;
                    tmp_raw = butterfilt(tmp_raw, params.filter, params.ax{axIdx}.head.sampRate, 'low', dim);
                    pOn_idx = sum(threshCrossing(:,l_valid),2); % error checking to make sure the pulses for each sweep are at the same time
                    [pks, pOn_idx]= findpeaks(double(pOn_idx));
                    assert(~isempty(pOn_idx), 'ERROR: could not locate pulses')
                    assert(all(pks == sum(l_valid)), 'ERROR: Pulse times may not be syncd or there may be different numbers of pulses each sweep')
                    
                    
                    % store the pulse on times for down stream analysis
                    params.(groupName).tdict{i_ax}.pOnIdx{i_cond, i_ch} = pOn_idx;
                    
                    
                    % now continue with the mean subtraction and averaging...
                    baselineSamples = ceil(0.150 .* params.ax{axIdx}.head.sampRate);
                    baseline_idx = false(Nsamps, 1);
                    baseline_idx(pOn_idx(1)-baselineSamples : pOn_idx(1)) = true;
                    
                    baseline_raw = mean(tmp_raw(baseline_idx,:),1); % the actual baseline
                    tmp_raw = bsxfun(@minus, tmp_raw, baseline_raw);
                    
                    params.(groupName).avg.trace_pA{i_ax}{i_cond, i_ch} = mean(tmp_raw,2); % average across sweeps.
                    
                    % store the holding potential (if Vclamp)
                    if ~isempty(secondaryChIdx)
                        idx = secondaryChIdx(i_ch);
                        vclamp = strcmpi(params.ax{axIdx}.head.recChUnits{idx}, 'mv');
                        if vclamp
                            tmp_raw = params.ax{axIdx}.dat(:,idx,l_valid);
                            tmp_raw = permute(tmp_raw, [1,3,2]);
                            vhold = mean(tmp_raw(baseline_idx,:),1); % vhold for each sweep
                            assert(range(vhold)<=1, 'ERROR: difference in Vhold exceeds tolerance');
                            
                            params.(groupName).avg.vhold{i_ax}{i_cond, i_ch} = mean(vhold);
                        end
                    else
                        params.(groupName).avg.vhold{i_ax}{i_cond, i_ch} = nan;
                    end
                    
                end
                
            end
            
        end % files
    end % groups
    
    
    %
    % PLOT THE AVERAGE RESPONSE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i_grp = 1:size(params.groups,1)
        
        groupName = params.groups{i_grp,1};
        
        figure
        set(gcf, 'name',groupName, 'position', [298 27 394 757]);
        
        % Figure out which files are in each group, and
        % the index to params.avg.trace_pA
        nFiles = numel(params.(groupName).avg.trace_pA);
        nTraces = 0;
        for i_fid = 1:nFiles
            nTraces = nTraces + size(params.(groupName).tdict{i_fid}.conds, 1);
        end
        clrs = colormap('lines'); 
        cidx = round(linspace(1,size(clrs,1), nTraces));
        clrs = clrs(cidx,:);
        cidx = 1;
        xmax = [0 0]; % one for each channel
        for i_fid = 1:nFiles;
            nConds = size(params.(groupName).avg.trace_pA{i_fid},1);
            
            for i_cond = 1:nConds
                
                nCh = size(params.(groupName).avg.trace_pA{i_fid},2);
                for i_ch = 1:nCh;
                    subplot(nCh,1,i_ch), hold on,
                    
                    if ~all(isnan(params.(groupName).avg.trace_pA{i_fid}{i_cond, i_ch})) % excludes conditions for which there were no sweeps
                        N = numel(params.(groupName).avg.trace_pA{i_fid}{i_cond, i_ch});
                        sampFreq = params.(groupName).head{i_fid}.sampRate;
                        tt = (0:N-1) ./ sampFreq;
                        t_on_idx = params.(groupName).tdict{i_fid}.pOnIdx{i_cond, i_ch}(1);
                        t_on_ms = t_on_idx./sampFreq;
                        t_off_idx = params.(groupName).tdict{i_fid}.pOnIdx{i_cond, i_ch}(end);
                        t_off_ms = t_off_idx ./ sampFreq;
                        xmax(i_ch) = max([t_off_ms, xmax(i_ch)]); % determine the new t_off
                        plot(tt, params.(groupName).avg.trace_pA{i_fid}{i_cond, i_ch}, '-', 'color', clrs(cidx,:))
                        xlim([t_on_ms-0.025, xmax(i_ch)+0.075])
                    end
                    
                    xlabel('time (sec)')
                    ylabel('current (pA)')
                    
                end
                
                % change the color on a condition by condition basis
                cidx = cidx+1;
            end
        end
    end
    
    
    
end % function


function out = findGoodSweeps(params, axIdx, i_ch)
    
    
    % determine the data files for this particular condition
    prefix = params.files{1,1};
    list = params.files{1,2};
    completeNames = {};
    for i_fid = 1:numel(list);
        suffix = num2str(list(i_fid));
        nZerosNeeded = 4-numel(suffix);
        suffix = [repmat('0',1,nZerosNeeded), suffix];
        completeNames{i_fid} = [prefix, suffix];
    end
    
    
    % ignore data files for specific channels if need be
    switch i_ch
        case 1
            possibleExclusions = params.excludeHS1;
        case 2
            possibleExclusions = params.excludeHS2;
    end

    exclude_file = 0;
    ex = 1;
    while ex <= numel(possibleExclusions)
        if iscell(possibleExclusions{ex})
            exclude_file = ~isempty(regexpi(completeNames{axIdx}(end-4:end), possibleExclusions{ex}{1}));
            if exclude_file
                exclude_sweeps = possibleExclusions{ex}{2};
            end
        else
            exclude_file = ~isempty(regexpi(completeNames{axIdx}(end-4:end), possibleExclusions{ex}));
            if exclude_file
                exclude_sweeps = 1:size(params.ax{axIdx}.dat,3);
            end
        end



        % break if there was a match
        if exclude_file
            break
        end

        % increment the counter
        ex = ex+1;
    end

    totalSweeps = size(params.ax{axIdx}.dat, 3);
    if exclude_file && (totalSweeps == numel(exclude_sweeps))
        fprintf('excluding ch %d from file %s \n', i_ch, completeNames{axIdx});
        out = 'exclude_file';
        
    elseif  exclude_file && (totalSweeps > numel(exclude_sweeps))
        fprintf('excluding %d sweeps from ch %d from file %s \n', numel(exclude_sweeps), i_ch, completeNames{axIdx});
        out = true(totalSweeps, 1);
        out(exclude_sweeps) = false;
        
    else
        out = true(totalSweeps, 1);
    end


end
