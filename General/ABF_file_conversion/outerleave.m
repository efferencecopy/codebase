function tDict = outerleave(stimWF, sampRate, BLACKROCKCORRECTION)
    %
    %  EXAMPLE    tDict = outerleave(stimWF, sampRate, [BLACKROCKCORRECTION])
    % 
    %  stimWF should be a [Ntime x 1 x Ntrials] matrix (clampex), or
    %         several cell arrays (blackrock).
    %
    %  BLACKROCKCORRECTION is optional, but should be set to TRUE when the
    %                      data come from blackrock
    %
    %
    %  tDict.vars     => {'pAmp', 'pWidth', 'tFreq', 'tRecov', 'tRIT'}  a reference to colums in tDict.conds
    %  tDict.conds    => matrix of values. each ROW corresponds to a unique type
    %  tDict.trlList  => matrix of scalars that map each sweep onto a 'condition', one column for each channel recorded
    %
    %  pAmp      is specified for all sweep types, in Volts
    %  pWidth    is specified for all sweep types, in Seconds
    %  tFreq     valid for normal and recovery trains. Equals zero for RITs
    %  tRecov    for recov trains, equals the recovery period. Else, equals zero
    %  tRIT      a unique number for each version of a RIT, else equals zero
    %
    %  ** Be advised that this function will round pulse widths to the nearest
    %     100 usec, becuse there is jitter in blacrock's timing of events
    %     relative to how clampex aquires signals.
    %
    % C.Hass 2015
    %

    
    if ~exist('BLACKROCKCORRECTION', 'var') || isempty(BLACKROCKCORRECTION)
        BLACKROCKCORRECTION = false;
    end

    % condition the inputs if the data come from a .abf file
    if ~iscell(stimWF)
        stimWF = permute(stimWF, [1,3,2]);
        stimWF = mat2cell(stimWF, size(stimWF, 1), ones(1,size(stimWF,2)));
    end


    Nsweeps = numel(stimWF);
    si = 1./sampRate; % the sample interval...

    % iterate over the sweeps determining the pulse amplitude, width, freq
    [pWidth, pAmp, tFreq] = deal(nan(Nsweeps, 1));
    tRecov = zeros(Nsweeps, 1); % using a numeric value as the default b/c nans will make 'unique' bonk later in the function
    tRIT = zeros(Nsweeps, 1); % zero for non-RIT trains, else, a unique number for each RIT version
    RIT_versions = {}; % to keep track of IPIs for each RIT version.
    for swp = 1:Nsweeps

        tmpWF = stimWF{swp};
        tmpWF = tmpWF(:);

        % amplitude first (which sets the threshold)
        pAmp(swp) = max(tmpWF);
        pAmp(swp) = round(pAmp(swp).*100) ./ 100; % round to the one hundreths place

        % now find pulse width and IPI
        thresh = pAmp(swp) .* 0.8;
        above = tmpWF > thresh;
        change = [0; diff(above)];
        xUp = change == 1;
        xDown = change == -1;

        tt = [0:numel(xUp)-1] ./ sampRate;
        pOnTimes = tt(xUp);
        pOffTimes = tt(xDown);



        tmpWidth = pOffTimes-pOnTimes;
        if range(tmpWidth) < (1e-6) && ~BLACKROCKCORRECTION% pulse by pulse diffs are < 1us, typically only for clampex data so exclude blackrock things
            tmpWidth = mean(tmpWidth);
            tmpWidth = round(tmpWidth,6); % round to keep things roughly consistent
        elseif all(tmpWidth > 100e-6) && BLACKROCKCORRECTION
            %blackrock timing is inconsistent pulse to pulse...but this fix
            %should only be used for pulses that are wide.
            tmpWidth = mean(tmpWidth);
            tmpWidth = tmpWidth .* 1e4; % in hundreds of usec.
            tmpWidth = floor(tmpWidth) ./ 1e4;
            assert(si<=50e-6, 'ERROR: pulse widths may be unreliable')
        else
            error('Found too many pWidths, but can not fix it')
        end

        pWidth(swp) = tmpWidth;

        % figure out if the stimulus was simple trains or random pulses. For
        % simple trains, there should only be 1 (or 2) different IPIs (2 if
        % there is a recovery train). For RITs, there should be several
        % different IPIs.
        is_train = numel(pOnTimes)>1;
        if is_train
            IPIs = diff(pOnTimes);
            unique_ipis = unique(round(IPIs, 6), 'stable'); % force unique to preserve order and not sort
            Nipi = numel(unique_ipis);
            
            % all trains have at least one unique IPI. Assume this
            % corresponds to the train freq
            train_tf =  1./unique_ipis(1);
            tFreq(swp) = round(train_tf); % round to the nearest whole number
            
            % determine the recovery time if the stimulus was Recovery
            % or Zucker
            is_zucker = Nipi == 3;
            is_recov = Nipi == 2;
            if is_zucker || is_recov
                recovIPI = round(unique_ipis(2).*1000); % in ms
                if recovIPI > ((1./tFreq(swp))+0.010)*1000 % needs to be 10 ms longer than the typical IPI
                    tRecov(swp) = recovIPI;
                else
                    error('ERROR: Nominal recovery pulse was <10ms from previous induction pulse')
                end
            end
            
            % if the stimulus was a RIT, then define the RIT_version
            % number
            is_rit = Nipi > 3;
            if is_rit
                [RIT_num, RIT_versions] = define_RIT_version(pOnTimes, RIT_versions);
                tRIT(swp) = RIT_num;
                tFreq(swp) = 0; % needs to be defined so that unique doesn't bonk below
            end
            
        else
            tFreq(swp) = 0; % a single pulse
        end
    end

    % now determine the number of unique trial types
    tDict.vars = {'pAmp', 'pWidth', 'tFreq', 'tRecov', 'tRIT'};
    tDict.conds = unique([pAmp, pWidth, tFreq, tRecov, tRIT], 'rows');

    Nconds = size(tDict.conds, 1);

    tDict.trlList = nan(Nsweeps, 1);
    for a = 1:Nconds
        tmp = tDict.conds(a,:);
        l_cond = ismember([pAmp, pWidth, tFreq, tRecov, tRIT], tmp, 'rows');
        tDict.trlList(l_cond) = a;
    end
    assert(~any(isnan(tDict.trlList)), 'ERROR: some of the trials were not defined');
end


function [RIT_num, RIT_versions] = define_RIT_version(pOnTimes, RIT_versions)

    % does this RIT match any of the previous versions? If so, figure out
    % which one it is. Otherwise, make a new RIT version number and store
    % the IPIs
    
    if isempty(RIT_versions)
        RIT_num = 1;
        RIT_versions{RIT_num} = pOnTimes;        
    else
        
        match = cellfun(@(x) (numel(x)==numel(pOnTimes)) && all(x == pOnTimes), RIT_versions);
        assert(sum(match)<=1, 'ERROR: too many matches found')
        
        if ~any(match)
            RIT_num = numel(RIT_versions)+1;
            RIT_versions{end+1} = pOnTimes;
        else
            RIT_num = find(match);
        end
    end
    
end


