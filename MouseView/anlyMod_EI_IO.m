function params = anlyMod_EI_IO(params)


    % Game plan:
    %
    % 1) Several files could be in the queue for analysis. They should already
    % be abf2obj'd, but they will likely have an unknown number of interleaved
    % stimulus types. I need some way of account for all the different stimulus
    % types.
    %
    %



    % iterate over the files.
    %
    % outerleave to retrieve the number of conditions
    %
    % go through each condition CH by CH. the list of valid channels should be
    % given by the outerleave tdict AND by the params.exclude field.
    %
    % calculate the average of each group, store it in the
    % params.isolatedCurrents field along with the holding potential.
    %
    % plot the raw values (baseline subtracted) and the average.
    %
    % calculate conductance on the fly (also send this to the isolatedData
    % fields)
    %
    % calculate peaks on a peak by peak basis


    for i_ax = 1:numel(params.ax)

        ledIdx = params.ax{i_ax}.idx.LED_470;
        tdict = outerleave(params.ax{i_ax}, ledIdx);
        
        
        
        
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
        assert(numel(primaryChIdx)<=2, 'ERROR: too many channels')
        secondaryChIdx = [find(l_hs1 & l_secCh), find(l_hs2 & l_secCh)];

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
                
                % grab the raw data, mean subtract, and then average across
                % sweeps
                tmp_raw = params.ax{i_ax}.dat(:,chIdx,l_valid);
                tmp_raw = permute(tmp_raw, [1,3,2]);
                pOn_idx = sum(threshCrossing(:,l_valid),2);
                [pks, pOn_idx]= findpeaks(double(pOn_idx));
                assert(~isempty(pOn_idx), 'ERROR: could not locate pulses')
                assert(all(pks == sum(l_valid)), 'ERROR: Pulse times may not be syncd or there may be different numbers of pulses each sweep')
                
                baselineSamples = ceil(0.150 .* params.ax{i_ax}.head.sampRate);
                baseline_idx = false(Nsamps, 1);
                baseline_idx(pOn_idx-baselineSamples : pOn_idx) = true;
                
                baseline_raw = mean(tmp_raw(baseline_idx,:),1); % the actual baseline
                tmp_raw = bsxfun(@minus, tmp_raw, baseline_raw);
                
                avg = mean(tmp_raw,2); % average across sweeps.
                
                

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
        fprintf('excluding ch %d from file %s \n', ch, completeNames{i_ax});
        out = 'exclude_file';
        
    elseif  exclude_file && (totalSweeps > numel(exclude_sweeps))
        fprintf('excluding %d sweeps from ch %d from file %s \n', numel(exclude_sweeps), ch, completeNames{i_ax});
        out = true(totalSweeps, 1);
        out(exclude_sweeps) = false;
        
    else
        out = true(totalSweeps, 1);
    end


end
