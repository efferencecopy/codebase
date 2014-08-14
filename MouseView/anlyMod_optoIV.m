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
        primaryChIdx = [find(l_hs1 & ~l_secCh), find(l_hs2 & ~l_secCh)];
        assert(numel(primaryChIdx)<=2, 'ERROR: too many channels')
        secondaryChIdx = [find(l_hs1 & l_secCh), find(l_hs2 & l_secCh)];
        
        % the time of the LED pulse is the same for both channels, so just
        % grab it now.
        ledChIdx = cellfun(@(x) ~isempty(x), regexpi(ax.head.DACchNames, 'LED'));
        [~, pulseTime] = ax.threshold(0.1, [find(ledChIdx), 1], 'u');
        
        % extract the data from each available channel. store the baseline
        % subtracted average
        preTime = 0.100;
        baselinePoints = preTime .* ax.head.sampRate;
        postTime = 0.300;
        
        for ch = 1:numel(primaryChIdx);
            
            % ignore data files for specific channels if need be
            switch ch
                case 1
                    possibleExclusions = params.excludeHS1;
                case 2
                    possibleExclusions = params.excludeHS2;
            end
            
            exclude_file = 0;
            ex = 1;
            while ex < numel(possibleExclusions)
                if iscell(possibleExclusions{ex})
                    exclude_file = ~isempty(regexpi(groupFiles{i}(end-4:end), possibleExclusions{ex}{1}));
                    if exclude_file
                        exclude_sweeps = possibleExclusions{ex}{2};
                    end
                else
                    exclude_file = ~isempty(regexpi(groupFiles{i}(end-4:end), possibleExclusions{ex}));
                    if exclude_file
                        exclude_sweeps = 1:size(ax.dat,3);
                    end
                end
                
                
                
                % break if there was a match
                if exclude_file
                    disp('found match')
                    break
                end
                
                % increment the counter
                ex = ex+1;
            end
            
            if exclude_file
                fprintf('excluding ch %d from file %s \n', ch, groupFiles{i});
                continue
            end
            
            % extract the data (assuming the file or sweeps didn't get
            % excluded.
            sweeps = ax.getvals(primaryChIdx(ch), 1:size(ax.dat,3), pulseTime-preTime, pulseTime+postTime);
            sweeps = permute(sweeps, [3,1,2]);
            baseline = mean(sweeps(:, 1:baselinePoints), 2);
            sweeps = bsxfun(@minus, sweeps, baseline);
            ivdat.(params.groups{a,1}).raw{ch}{i} = mean(sweeps, 1);
            
            % grab the holding potential here, store it in ivdat
            vdat = ax.getvals(secondaryChIdx(ch), 1:size(ax.dat,3), pulseTime-preTime, pulseTime+postTime);
            vdat = mean(vdat, 1); % mean across time;
            vdat = mean(vdat, 3); % mean across sweeps
            ivdat.(params.groups{a,1}).vhold{ch}{i} = vdat;
        end
    end
end


% pause here and plot the raw data for each channel. Make a new
% plot, and shift them vertically according to the holding
% potential.

groups = fieldnames(ivdat);
for a = 1:numel(groups)
    
    figure
    set(gcf, 'name', groups{a}, 'position', [24    10   560   779]);
    Nchannels = size(ivdat.(groups{a}).raw, 2);
    
    for ch = 1:Nchannels
        
        subplot(Nchannels, 1, ch), hold on,
        xlabel('time (sec)')
        ylabel('current (pA)')
        Nvholds = numel(ivdat.(groups{a}).raw{ch});
        vholds = cat(1, ivdat.(groups{a}).vhold{ch}{:});
        [~,vholdOrder] = sort(vholds);
        clrs = pmkmp(Nvholds,'IsoL');
        clrs = clrs(randperm(Nvholds), :);
        legtext = {};
        
        for i = 1:Nvholds;
            idx = vholdOrder(i);
            pA = ivdat.(groups{a}).raw{ch}{idx};
            mV = ivdat.(groups{a}).vhold{ch}{idx};
            tvec = linspace(-preTime, postTime, numel(pA));
            plot(tvec, pA+mV, '-', 'linewidth', 2, 'color', clrs(i,:))
            legtext{i} = sprintf('%.2f mV', mV);
        end
        legend(legtext)
        legend boxoff
    end
end

% store some useful things in the params structure, and return to the
% calling function
ivdat.tvec = tvec;
params.ivdat = ivdat;


    
    
    

