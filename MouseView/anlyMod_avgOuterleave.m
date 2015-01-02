function params = anlyMod_avgOuterleave(params)

    %
    % CALCULATE THE AVERAGE RESPONSE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i_ax = 1:numel(params.ax)

        ledIdx = params.ax{i_ax}.idx.LED_470;
        tdict = outerleave(params.ax{i_ax}, ledIdx);
        params.tdict{i_ax} = tdict;       
        
        
        % calculate the pulse times, which will be used later
        tmpWF = params.ax{i_ax}.dat(:,ledIdx,:);
        tmpWF = permute(tmpWF, [1,3,2]);
        [Nsamps, Nsweeps] = size(tmpWF);
        thresh = 0.025;
        above = tmpWF > thresh;
        change = [zeros(1,Nsweeps); diff(above, 1, 1)];
        threshCrossing = change == 1;


        % figure out how many channels there were
        l_secCh = cellfun(@(x) ~isempty(x), regexpi(params.ax{i_ax}.head.recChNames, 'sec'));
        l_hs1 = cellfun(@(x) ~isempty(x), regexpi(params.ax{i_ax}.head.recChNames, 'hs1'));
        l_hs2 = cellfun(@(x) ~isempty(x), regexpi(params.ax{i_ax}.head.recChNames, 'hs2'));
        primaryChIdx = [find(l_hs1 & ~l_secCh), find(l_hs2 & ~l_secCh)];
        secondaryChIdx = [find(l_hs1 & l_secCh), find(l_hs2 & l_secCh)];
        assert(numel(primaryChIdx)<=2, 'ERROR: too many channels')

        for i_ch = 1:numel(primaryChIdx);

            chIdx = primaryChIdx(i_ch);

            % generate a trial list for this specific channel and .abf file
            l_goodSweeps = findGoodSweeps(params, i_ax, i_ch);
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
                    continue
                end
                
                % grab the raw data, mean subtract, and then average across sweeps
                tmp_raw = params.ax{i_ax}.dat(:,chIdx,l_valid);
                tmp_raw = permute(tmp_raw, [1,3,2]);
                pOn_idx = sum(threshCrossing(:,l_valid),2);
                [pks, pOn_idx]= findpeaks(double(pOn_idx));
                assert(~isempty(pOn_idx), 'ERROR: could not locate pulses')
                assert(all(pks == sum(l_valid)), 'ERROR: Pulse times may not be syncd or there may be different numbers of pulses each sweep')
                
                
                % store the pulse on times for down stream analysis
                params.tdict{i_ax}.pOnIdx{i_cond, i_ch} = pOn_idx;
                
                
                % now continue with the mean subtraction and averaging...
                baselineSamples = ceil(0.150 .* params.ax{i_ax}.head.sampRate);
                baseline_idx = false(Nsamps, 1);
                baseline_idx(pOn_idx(1)-baselineSamples : pOn_idx(1)) = true;
                
                baseline_raw = mean(tmp_raw(baseline_idx,:),1); % the actual baseline
                tmp_raw = bsxfun(@minus, tmp_raw, baseline_raw);
                
                params.avg.trace_pA{i_ax}{i_cond, i_ch} = mean(tmp_raw,2); % average across sweeps.
                
                % store the holding potential (if Vclamp)
                idx = secondaryChIdx(i_ch);
                vclamp = strcmpi(params.ax{i_ax}.head.recChUnits{idx}, 'mv');
                if vclamp
                    tmp_raw = params.ax{i_ax}.dat(:,idx,l_valid);
                    tmp_raw = permute(tmp_raw, [1,3,2]);
                    vhold = mean(tmp_raw(baseline_idx,:),1); % vhold for each sweep
                    assert(range(vhold)<=1, 'ERROR: difference in Vhold exceeds tolerance');
                    
                    params.avg.vhold{i_ax}{i_cond, i_ch} = mean(vhold);
                end

            end

        end

    end
    
    
    %
    % PLOT THE AVERAGE RESPONSE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i_grp = 1:size(params.groups,1)
        
        figure
        set(gcf, 'name', params.groups{i_grp,1}, 'position', [298 27 394 757]);
        
        % stuff in the params.avg.trace_pA array is ordered in the same way
        % as params.files. Figure out which files are in each group, and
        % the index to params.avg.trace_pA
        nFiles = numel(params.groups{i_grp,2});
        clrs = pmkmp(nFiles+1,'IsoL'); % pmkmp bonks when asked for a single color. adding one to avoid the bonking...
        clrs = clrs(randperm(nFiles), :);
        for i_fid = 1:nFiles;
            
            fid_idx = params.files{2} == params.groups{i_grp,2}(i_fid);
            nConds = size(params.avg.trace_pA{fid_idx},1);
            
            for i_cond = 1:nConds
                
                nCh = size(params.avg.trace_pA{fid_idx},2);
                for i_ch = 1:nCh;
                    subplot(nCh,1,i_ch), hold on,
                    try
                    plot(params.ax{fid_idx}.tt, params.avg.trace_pA{fid_idx}{i_cond, i_ch}, '-', 'color', clrs(i_fid,:))
                    catch
                        keyboard
                    end
                    t_on = params.ax{fid_idx}.tt(params.tdict{fid_idx}.pOnIdx{i_cond, i_ch}(1));
                    t_off = params.ax{fid_idx}.tt(params.tdict{fid_idx}.pOnIdx{i_cond, i_ch}(end));
                    xlim([t_on, t_off+0.100])
                    xlabel('time (sec)')
                    ylabel('current (pA)')
                end
            end
        end
    end
    
    
    
end % function


function out = findGoodSweeps(params, i_ax, i_ch)
    
    
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
            exclude_file = ~isempty(regexpi(completeNames{i_ch}(end-4:end), possibleExclusions{ex}{1}));
            if exclude_file
                exclude_sweeps = possibleExclusions{ex}{2};
            end
        else
            exclude_file = ~isempty(regexpi(completeNames{i_ch}(end-4:end), possibleExclusions{ex}));
            if exclude_file
                exclude_sweeps = 1:size(ax.dat,3);
            end
        end



        % break if there was a match
        if exclude_file
            break
        end

        % increment the counter
        ex = ex+1;
    end

    totalSweeps = size(params.ax{i_ax}.dat, 3);
    if exclude_file && (totalSweeps == numel(exclude_sweeps))
        fprintf('excluding ch %d from file %s \n', i_ch, completeNames{i_ax});
        out = 'exclude_file';
        
    elseif  exclude_file && (totalSweeps > numel(exclude_sweeps))
        fprintf('excluding %d sweeps from ch %d from file %s \n', numel(exclude_sweeps), i_ch, completeNames{i_ax});
        out = true(totalSweeps, 1);
        out(exclude_sweeps) = false;
        
    else
        out = true(totalSweeps, 1);
    end


end
