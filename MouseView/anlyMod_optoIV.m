function params = anlyMod_optoIV(params)

% a structure array for the IV curve data:
ivdat = []; 

% construct a library of the available data files
abfLibrary = {};
for a = 1:numel(params.ax)
    abfLibrary{a} = params.ax{a}.name;
end


% iterate over groups pulling out the raw data, and making an average trace
% centered on the LED pulse
Ngroups = size(params.groups,1);
for a = 1:Ngroups
   
    % determine the data files for this particular condition
    prefix = params.files{1};
    list = params.groups{a,2};
    groupFiles = {};
    for i = 1:numel(list);
        suffix = num2str(list(i));
        nZerosNeeded = 4-numel(suffix);
        suffix = [repmat('0',1,nZerosNeeded), suffix];
        groupFiles{i} = [prefix, suffix];
    end
    
    % iterate over the data files and extract stuff for HS1 and HS2
    for i = 1:numel(groupFiles);
        fIdx = cellfun(@(x) ~isempty(x), regexpi(groupFiles{i}, abfLibrary));
        ax = params.ax{fIdx};
        
        % iterate over recorded channels. figure out how many channels
        % there were.
        l_secCh = cellfun(@(x) ~isempty(x), regexpi(ax.head.recChNames, 'sec'));
        l_hs1 = cellfun(@(x) ~isempty(x), regexpi(ax.head.recChNames, 'hs1'));
        l_hs2 = cellfun(@(x) ~isempty(x), regexpi(ax.head.recChNames, 'hs2'));
        recChIdx = [find(l_hs1 & ~l_secCh), find(l_hs2 & ~l_secCh)];
        assert(numel(recChIdx)<=2, 'ERROR: too many channels')
        
        % the time of the LED pulse is the same for both channels, so just
        % grab it now.
        ledChIdx = cellfun(@(x) ~isempty(x), regexpi(ax.head.DACchNames, 'LED'));
        [~, pulseTime] = ax.threshold(0.1, [find(ledChIdx), 1], 'u');
        
        % extract the data from each available channel. store the baseline
        % subtracted average
        preTime = 0.100;
        baselinePoints = preTime .* ax.head.sampRate;
        postTime = 0.300;
        for ch = 1:numel(recChIdx);
            sweeps = ax.getvals(recChIdx(ch), 1:size(ax.dat,3), pulseTime-preTime, pulseTime+postTime);
            sweeps = permute(sweeps, [3,1,2]);
            baseline = mean(sweeps(:, 1:baselinePoints), 2);
            sweeps = bsxfun(@minus, sweeps, baseline);
            ivdat.(params.groups{a,1}).raw{i, ch} = mean(sweeps, 1);
            
            % grab the holding potential here, store it in ivdat
            
        end
    end
end


% pause here and plot the raw data for each channel. Make a new
% plot, and shift them vertically according to the holding
% potential.

% now grab the PCS magnitude. I don't know how to deal with the
% mixed E and I currents yet...
% determine the appropriate time region to analyze, specify this in
% the physiology notes, or use the ginput field for each channel.
% only take the

    
    
    

