%% designate an experiment to analyze, import the data

error('need to review this code to make sure it is doing what I want')

fin

NSX = 'ns5'; % the version with the LFP continuous data
STIMTYPE = 'train'; % train or sinusoid
NOISEMETHOD = 'none'; % 'subtract' or 'filter'


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Load in the experimental metadata
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~, ~, wb_info] = xlsread('invivoOptostimMetadata.xlsx', 1);
wbidx.mouse_name = strcmpi(wb_info(1,:), 'mouse name');
wbidx.site = strcmpi(wb_info(1,:), 'site');
wbidx.analyze = strcmpi(wb_info(1,:), 'include in analysis');
wbidx.dat_fname = strcmpi(wb_info(1,:), 'data file');
wbidx.shank_pos = strcmpi(wb_info(1,:), 'distance to pia for first contact per shank (um)');
wbidx.shank_area = strcmpi(wb_info(1,:), 'brain area for each shank');

% fix the shank related idx so that all four shanks get analyzed
wbidx.shank_pos(find(wbidx.shank_pos):find(wbidx.shank_pos)+3) = true;
wbidx.shank_area(find(wbidx.shank_area):find(wbidx.shank_area)+3) = true;

% delete the headers
wb_info(1:2,:) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Iterate over all the mice and pull in the raw data for each recording
% location. Concatenate as need be. Save analysis params, and the
% down-sampled data to disk (with a date stamp)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dat = {};
Nexpts = size(wb_info, 1);
for i_ex = 4:Nexpts;
    
    fprintf('Loading in %s, site: %d\n', wb_info{i_ex,wbidx.mouse_name}, wb_info{i_ex, wbidx.site})
    
    tmp_names = wb_info(i_ex, wbidx.dat_fname);
    fnames = {};
    if tmp_names{1}(1) == '{'
        eval(sprintf('tmp_names = %s;', tmp_names{1}));
        for i_fid = 1:numel(tmp_names{2})
            fnames{i_fid} = [tmp_names{1}, num2str(tmp_names{2}(i_fid))];
        end
    else
        fnames{1} = tmp_names{1};
    end
    
    clear tmp % start fresh each time
    for i_fid = 1:numel(fnames)
        tmp{i_fid} = blkobj(fnames{i_fid}, blkcodes_pulseTrains);
    end
    
    blk = tmp{1};
    
    % remove the SU and MU channels
    anlg_ch_idx = cellfun(@(x,y) ~isempty(regexpi(x, y)), blk.sum.rasterFields, repmat({NSX}, size(blk.sum.rasterFields)));
    blk.ras(:, ~anlg_ch_idx) = [];
    blk.sum.rasterFields(~anlg_ch_idx) = [];
    
    for i_fid = 2:numel(tmp)
        
        % remove non-lfp channels
        anlg_ch_idx = cellfun(@(x,y) ~isempty(regexpi(x, y)), tmp{i_fid}.sum.rasterFields, repmat({NSX}, size(tmp{i_fid}.sum.rasterFields)));
        tmp{i_fid}.ras(:, ~anlg_ch_idx) = [];
        tmp{i_fid}.sum.rasterFields(:, ~anlg_ch_idx) = [];
        
        match = cellfun(@(x,y) strcmp(x,y), blk.sum.trialFields, tmp{i_fid}.sum.trialFields);
        assert(all(match), 'ERROR: trial fields do not match');
        match = cellfun(@(x,y) strcmp(x,y), blk.sum.rasterFields, tmp{i_fid}.sum.rasterFields);
        assert(all(match), 'ERROR: raster fields do not match');
        
        blk.ras = cat(1, blk.ras, tmp{i_fid}.ras);
        blk.trial = cat(1, blk.trial, tmp{i_fid}.trial);
    end
    
    clear tmp; % clear the big array
    
    % remake the stro.sum.idx fields
    blk.sum.idx = [];
    N = numel(blk.sum.rasterFields);
    for a = 1:N
        tmp = blk.sum.rasterFields{a};
        blk.sum.idx.(tmp) = a;
    end
    N = numel(blk.sum.trialFields);
    for a = 1:N
        tmp = blk.sum.trialFields{a};
        blk.sum.idx.(tmp) = a;
    end
    
    
    % re-define these logical indicies b/c they may have changed
    anlg_ch_idx = cellfun(@(x,y) ~isempty(regexpi(x, y)), blk.sum.rasterFields, repmat({NSX}, size(blk.sum.rasterFields)));
    neural_ch_idx = cellfun(@(x,y) ~isempty(regexpi(x, y)), blk.sum.rasterFields, repmat({'chan'}, size(blk.sum.rasterFields)));
    lfp_ch_idx = neural_ch_idx & anlg_ch_idx;
    sampFreq_nsx = blk.sum.(NSX).MetaTags.SamplingFreq;
    
    
    % add the header info to the data structure output
    dat{i_ex}.sum = blk.sum;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % PRE-PROCESSING FOR LFP ANALYSIS
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Define a new sampling rate for the lfp analysis. Pick one that will allow
    % for the freq components of interest to be preserved and one that's evenly
    % divided into sampFreq_nsx. The LFP lowpass filter will be constrained
    % to be <= sampFreq_lfp / 5;
    sampFreq_lfp = 15000;
    
    pool = gcp('nocreate');
    if isempty(pool)
        pool = parpool(32);
    end
    
    dat{i_ex} = extractLFPandMU(blk, lfp_ch_idx, sampFreq_lfp, NOISEMETHOD, STIMTYPE, NSX);
    
end

%% PLOT: AVERAGE LFP SIGNAL

PLOTTYPE = 'mean'; % 'errbar', 'all', 'mean'

% make a tmp variable (duplicate of the preprocessed LFP data)
trialSnips = dat{5}.lfpsnips;

% pull in the map of electrode positions
trodeMap = etrodemap('a32-4x8');
trodePltIdx = reshape(trodeMap', [], 1);


ptypes = fieldnames(trialSnips);
i_bad = strcmpi(ptypes, 'preTime') | strcmpi(ptypes, 'postTime') | strcmpi(ptypes, 'pulseOnT_sec');
ptypes(i_bad) = [];
Nptypes = numel(ptypes);

for i_type = 1:Nptypes
    
    fldname_ptype = ptypes{i_type};
    
    f = figure;
    f.Position = [45 -3 1090 789];
    for i_ch = 1:numel(trialSnips.(fldname_ptype));
        
        tmp = trialSnips.(fldname_ptype){i_ch};
        xbar = nanmean(tmp, 1);
        sem = stderr(tmp, 1);
        
        pltidx = find(trodePltIdx == i_ch);
        subplot(8, 4, pltidx, 'align')
        hold on,
        tt = ([0:numel(xbar)-1] ./ sampFreq_lfp) - trialSnips.preTime;
        
        switch PLOTTYPE
            case 'all'
                plot(tt, tmp')
            case 'mean'
                plot(tt(:), xbar(:), 'k');
            case 'errbar'
%                 plot(tt(:), xbar(:), 'k', 'linewidth' , 1);
%                 plot(tt(:), xbar(:)+sem(:), ':k');
%                 plot(tt(:), xbar(:)-sem(:), ':k');
                shadedErrorBar(tt(:), xbar(:), sem(:), {}, 1);
        end
        
        axis tight
        drawnow
        
    end
end



%% CURRENT SOURCE DENSITY ANALYSIS
close all; clc


METHOD = 'kcsd'; % could be 'normal', or 'kcsd'

EX_NUM = 6;

% make a tmp variable (duplicate of the preprocessed LFP data)
csdSnips = dat{EX_NUM}.lfpsnips;



%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load the metadata
%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~, ~, wb_info] = xlsread('invivoOptostimMetadata.xlsx', 1);
wbidx.mouse_name = strcmpi(wb_info(1,:), 'mouse name');
wbidx.site = strcmpi(wb_info(1,:), 'site');
wbidx.analyze = strcmpi(wb_info(1,:), 'include in analysis');
wbidx.dat_fname = strcmpi(wb_info(1,:), 'data file');
wbidx.shank_pos = strcmpi(wb_info(1,:), 'distance to pia for top contact per shank (um)');
wbidx.shank_area = strcmpi(wb_info(1,:), 'brain area for each shank');

% fix the shank related idx so that all four shanks get analyzed
wbidx.shank_pos(find(wbidx.shank_pos):find(wbidx.shank_pos)+3) = true;
wbidx.shank_area(find(wbidx.shank_area):find(wbidx.shank_area)+3) = true;

% delete the headers
wb_info(1:2,:) = [];

% define the etrode depths
depths = wb_info(EX_NUM,wbidx.shank_pos);
if any(cellfun(@isempty, depths))
    elPos_y = repmat([0.030:0.100:.730]', 1, 4);
    warning('depths not defined')
else
    depths_mm = -1.*[depths{:}]./1000;
    depths_mm = bsxfun(@plus, repmat([0:0.100:.700]', 1, 4), depths_mm);
end




%
% Now do the data analysis
%
ptypes = fieldnames(csdSnips);
i_bad = strcmpi(ptypes, 'preTime') | strcmpi(ptypes, 'postTime') | strcmpi(ptypes, 'pulseOnT_sec');
ptypes(i_bad) = [];

% rearrange the data according to the location of each channel on the
% electrode array {shank}[position, Ntime, Ntrials]
trodeMap = etrodemap('a32-4x8');  % pull in the map of electrode positions
csd_shank = [];
csd_final = [];
for i_cond = 1:numel(ptypes)
    
    % initialize the output
    Ntime = unique(cellfun(@(x) size(x,2), csdSnips.(ptypes{i_cond})));
    Ntrials = unique(cellfun(@(x) size(x,1), csdSnips.(ptypes{i_cond})));
    assert(numel(Ntime) == 1, 'ERROR: did not identify the number of time points')
    assert(numel(Ntrials) == 1, 'ERROR: did not identify the number of trials')
    for i_shank = 1:4;
        csd_shank.(ptypes{i_cond}){i_shank} = nan(8,Ntime,Ntrials);
    end
    
    % rearrange the data
    for i_ch = 1:numel(csdSnips.(ptypes{i_cond}))
        
        % find the new position
        [pos, shank] = find(trodeMap== i_ch);
        
        % rearrange the data
        tmp_data = csdSnips.(ptypes{i_cond}){i_ch};
        tmp_data = permute(tmp_data, [3,2,1]);
        
        % put the new data into the csd_shank array
        csd_shank.(ptypes{i_cond}){shank}(pos,1:Ntime,1:Ntrials) = tmp_data;
        
    end
    
    % compute the 2nd spatial derivitave, then average
    for i_shank = 1:4
        
        % compute the mean across trials
        tmp_data = csd_shank.(ptypes{i_cond}){i_shank};
        tmp_data = mean(tmp_data,3); % mean across trials
        
        % re-subtract the baseline level from the trial avgerage
        preSamps = round(csdSnips.preTime .* sampFreq_lfp);
        bkgnd = mean(tmp_data(:, 1:preSamps), 2);
        tmp_data = bsxfun(@minus, tmp_data, bkgnd);
        
        % store the final processed version of the LFP
        lfp_final.(ptypes{i_cond}){i_shank} = tmp_data;
        
        
        switch METHOD
            case 'normal'
                % compute the 2nd spatial derivitive
                csd_final.(ptypes{i_cond}){i_shank} = -diff(tmp_data, 2, 1); % difference in the spatial dimension (down columns) on a 3D array
                
                % upsample in time and space
                csd_final.(ptypes{i_cond}){i_shank} = interp2(csd_final.(ptypes{i_cond}){i_shank}, 'linear');
                
                % remove oversampling in the time dimension
                newN = size(csd_final.(ptypes{i_cond}){i_shank}, 2);
                idx = 1:2:newN;
                csd_final.(ptypes{i_cond}){i_shank} = csd_final.(ptypes{i_cond}){i_shank}(:,idx);
                
            case 'kcsd'
                
                elPos = depths_mm(:, i_shank)';
                X = 0 : 0.010 : max(elPos)+0.050;
                pots = tmp_data;
                
                % delete recording sites that are outside the brain
                l_oob = elPos < 0;
                pots(l_oob,:) = [];
                elPos(l_oob) = [];
                
                k = kCSD1d(elPos, pots, 'X', X);
                k.estimate;
                csd_final.(ptypes{i_cond}){i_shank}.csd = -1 .* k.csdEst; % invert the polarity of source sink so that warm colors are sinks
                csd_final.(ptypes{i_cond}){i_shank}.elPos = k.elPos;
                csd_final.(ptypes{i_cond}){i_shank}.X = k.X;
        end
        
    end
end

% Try to compute the 2D CSD
if strcmpi(METHOD, 'kcsd');
    
    for i_cond = 1:numel(ptypes)
        
        elPos_y = depths_mm;
        elPos_y = elPos_y(:);
        elPos_x = repmat([0:0.400:1.20], 8, 1);
        elPos_x = elPos_x(:);
        elPos = [elPos_x, elPos_y];
        
        pots = [];
        for i_shank = 1:4
            pots = cat(1, pots, lfp_final.(ptypes{i_cond}){i_shank});
        end
        
        % delete any recording sites that are outside the brain
        l_oob = elPos_y < 0;
        elPos(l_oob, :) = [];
        pots(l_oob,:) = [];
        
        [X, Y] = meshgrid([-0.100:0.005:1.30], [0:0.005:max(elPos_y)+0.025]);
        
        k = kcsd2d(elPos, pots, 'X', X, 'Y', Y);
        k.plot_CSD;
        
        % fix the plot labels
        f = gcf;
        f.Name = sprintf('Expt %d, ptype: %s', EX_NUM, ptypes{i_cond});
    end
end


%
% plot the resulting CSD
%%%%%%%%%%%%
for i_cond = 1:numel(ptypes)
    
    f = figure;
    f.Position = [233 18 679 768];
    f.Name = ptypes{i_cond};
    
    for i_shank = 1:4
        
        s=subplot(4,1,i_shank);
        
        % add a few lines of data (Duplicates) b/c pcolor destroys a few
        % lines
        plotDat = csd_final.(ptypes{i_cond}){i_shank}.csd;
        plotDat = [plotDat(1,:); plotDat];
        plotDat = [plotDat, plotDat(:,1)];
        
        newN = size(plotDat, 2);
        tt_sec = ([0:newN-1] ./ sampFreq_nsx) - csdSnips.preTime;
        pos = linspace(0,8,size(plotDat, 1));
        switch METHOD
            case 'normal'
                pcolor(tt_sec, pos, flipud(plotDat));
                shading interp;
                s.YTickLabel = flipud(s.YTickLabel);
            case 'kcsd'
                imagesc(plotDat);
        end
        
    end
end


%
% plot the average CSD to the first pulse (across TFs)
%%%%%%%%%%%%
f = figure;
f.Position = [361 17 242 758];
for i_shank = 1:4
    
    t_zero_idx = ceil(csdSnips.preTime .* sampFreq_lfp);
    preSamps = ceil(0.005 .* sampFreq_lfp);
    postSamps = ceil(0.024 .* sampFreq_lfp);
    idx = [t_zero_idx - preSamps : t_zero_idx + postSamps];
    
    % concatenate across TF conditions
    plotDat = [];
    for i_cond = 1:numel(ptypes)
        tmptf = str2double(ptypes{i_cond}(end-1:end));
        if tmptf <= 40
            plotDat = cat(3, plotDat, csd_final.(ptypes{i_cond}){i_shank}.csd(:,idx));
        end
    end
    
    % average across TFs
    plotDat = mean(plotDat, 3);
    newN = size(plotDat,2);
    tt_sec = ((0:newN-1)-preSamps) ./ sampFreq_lfp;
    s=subplot(4,1,i_shank);
    
    switch METHOD
        case 'normal'
            % add row/col that will be deleted by pcolor
            plotDat = [plotDat(1,:); plotDat];
            plotDat = [plotDat, plotDat(:,1)];
            
            
            % do the plotting
            pos = linspace(0,size(plotDat, 1),size(plotDat, 1));
            p=pcolor(tt_sec, pos, flipud(plotDat));
            shading interp;
            s.YTickLabel = flipud(s.YTickLabel);
        case 'kcsd'
            imagesc(plotDat);
            pos = csd_final.(ptypes{i_cond}){i_shank}.X;
            pos = [min(pos) : 0.100 : max(pos)]; 
    end
    colorbar
end



%
% plot the finalized LFP traces (as a sanity check)
%%%%%%%%%%%%
for i_cond = 1:numel(ptypes)
    
    f = figure;
    f.Position = [233 18 679 768];
    f.Name = ptypes{i_cond};
    
    for i_shank = 1:4
        
        s=subplot(4,1,i_shank);
        
        % add a few lines of data (Duplicates) b/c pcolor destroys a few
        % lines
        plotDat = lfp_final.(ptypes{i_cond}){i_shank};
        
        newN = size(plotDat, 2);
        tt_sec = ([0:newN-1] ./ sampFreq_lfp) - csdSnips.preTime;
        plot(tt_sec, plotDat')
        
    end
end


%% CROSS-CORRELATION ANALYSIS

% the motivation is to assess how well different layers track one-another,
% if L2/3 receives the input from V1, then it's possible that the 



%% PLOT PSTH OF MULTI-UNIT SPIKE TIMES

if ~exist('dat', 'var')
    fin
    load 'invivo_dat_160205.mat'
end


BINWIDTH = .250e-3;
DELETEARTIFACTS = false; % deletes spike times in a specified window immediately after the optostim
TOOSOONTIME = 0.5e-3;
STANDARDIZE_YLIMS = true;
DISPLAYDEPTH = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load the metadata
%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~, ~, wb_info] = xlsread('invivoOptostimMetadata.xlsx', 1);
wbidx.mouse_name = strcmpi(wb_info(1,:), 'mouse name');
wbidx.site = strcmpi(wb_info(1,:), 'site');
wbidx.analyze = strcmpi(wb_info(1,:), 'include in analysis');
wbidx.dat_fname = strcmpi(wb_info(1,:), 'data file');
wbidx.shank_pos = strcmpi(wb_info(1,:), 'distance to pia for first contact per shank (um)');
wbidx.shank_area = strcmpi(wb_info(1,:), 'brain area for each shank');

% fix the shank related idx so that all four shanks get analyzed
wbidx.shank_pos(find(wbidx.shank_pos):find(wbidx.shank_pos)+3) = true;
wbidx.shank_area(find(wbidx.shank_area):find(wbidx.shank_area)+3) = true;

% delete the headers
wb_info(1:2,:) = [];

uniqueMice = unique(wb_info(:, wbidx.mouse_name));

% rearrange the data according to position on the electrode array
trodeMap = etrodemap('a32-4x8');  % pull in the map of electrode positions
[Nsites, Nshanks] = size(trodeMap);


for i_mouse = 2;%1:numel(uniqueMice)
    
    % determine which data files correspond to this mouse
    l_mouse = strcmpi(uniqueMice{i_mouse}, wb_info(:, wbidx.mouse_name));
    l_valid = logical(cell2mat(wb_info(:, wbidx.analyze)));
    l_mouse = find(l_mouse & l_valid);
    Nrecpos = numel(l_mouse);
    
    
    mu_shank = [];
    mu_shank_info = {};
    for i_pos = 1:Nrecpos;
        
        idx = l_mouse(i_pos);
        tmp_spikeTimes = dat{idx}.spikeTimes;
        ptypes = fieldnames(tmp_spikeTimes);
        i_bad = strcmpi(ptypes, 'preTime') | strcmpi(ptypes, 'postTime') | strcmpi(ptypes, 'pulseOnT_sec');
        ptypes(i_bad) = [];
        tfidx = cellfun(@(x,y) regexpi(x, '_tf'), ptypes, 'uniformoutput', false);
        ptypes_tfonly = cellfun(@(x,y) x(y+1:end), ptypes, tfidx, 'uniformoutput', false);
        Nptypes = numel(ptypes);
        
        % grab the meta-data for this recording location
        mu_shank_info{i_pos}.firstSiteDepth = wb_info(idx, wbidx.shank_pos);
        mu_shank_info{i_pos}.brainAreas = wb_info(idx, wbidx.shank_area);
        
        for i_type = 1:Nptypes
            
            for i_ch = 1:numel(tmp_spikeTimes.(ptypes{i_type}).tt)
                
                % find the new position
                [site, shank] = find(trodeMap == i_ch);
                
                % grab the data
                type_spiketimes = tmp_spikeTimes.(ptypes{i_type}).tt{i_ch};
                
                if DELETEARTIFACTS
                    pulseOn_sec = tmp_spikeTimes.pulseOnT_sec.(ptypes{i_type});
                    for i_pulse = 1:numel(pulseOn_sec)
                        l_artifacts = (type_spiketimes>pulseOn_sec(i_pulse)) & (type_spiketimes < pulseOn_sec(i_pulse) + TOOSOONTIME);
                        type_spiketimes(l_artifacts) = [];
                    end
                end
                
                % put the data into the appropriate position
                mu_shank.(ptypes_tfonly{i_type}){i_pos}.spiketimes{site, shank} = type_spiketimes;
                mu_shank.(ptypes_tfonly{i_type}){i_pos}.Ntrials{site, shank} = tmp_spikeTimes.(ptypes{i_type}).Ntrials{i_ch};
                mu_shank.(ptypes_tfonly{i_type}){i_pos}.trialIDs = tmp_spikeTimes.(ptypes{i_type}).trialIDs{i_ch};
                mu_shank.(ptypes_tfonly{i_type}){i_pos}.pulseOnT_sec = tmp_spikeTimes.pulseOnT_sec.(ptypes{i_type});
                
            end
            
        end
        
    end
    
    
    %
    % Plotting: PSTHs
    %
    startTime = dat{idx}.spikeTimes.preTime;
    hFig = figure;
    hFig.Units = 'normalized';
    hFig.Position = [0 0 1 1];
    hTabGroup = uitabgroup(hFig);
    for i_type = 1:numel(ptypes_tfonly)
        
        maxTime = max(mu_shank.(ptypes_tfonly{i_type}){i_pos}.pulseOnT_sec);
        endTime = maxTime + startTime;
        bins = -startTime : BINWIDTH : endTime;
        
        % make a tab
        htab(i_type) = uitab(hTabGroup,'Title',ptypes_tfonly{i_type});
        hax(i_type) = axes('Parent', htab(i_type));
        Nrows = 8;
        Ncols = 4*Nrecpos;
        hBar = [];
        for i_pos = 1:Nrecpos
            for i_shank = 1:Nshanks;
                for i_site = 1:Nsites;
                    
                    pltDat = mu_shank.(ptypes_tfonly{i_type}){i_pos}.spiketimes{i_site, i_shank};
                    Ntrials = mu_shank.(ptypes_tfonly{i_type}){i_pos}.Ntrials{i_site, i_shank};
                    
                    if isempty(pltDat);
                        pltDat = inf; % a trick to make fake data that won't plot
                    end
                    
                    counts= histc(pltDat, bins);
                    counts = (counts./Ntrials) ./ BINWIDTH;
                    
                    plotCol = ((i_pos - 1) .* 4) + i_shank;
                    plotRow = i_site;
                    pltidx = sub2ind([Ncols, Nrows], plotCol, plotRow); % transpose row and col to get the correct linear index for plots
                    hBar(end+1) = subplot(Nrows, Ncols, pltidx, 'align');
                    bar(bins*1000, counts, 0.8, 'histc');                    
                    ylim([0, max(counts).* 1.1+eps])
                    axis tight
                    if i_site == 1
                        title(sprintf('area: %s', mu_shank_info{i_pos}.brainAreas{i_shank}))
                    end
                    
                    if DISPLAYDEPTH
                        % figure out the depth of the etrode
                        firstSiteDepth = mu_shank_info{i_pos}.firstSiteDepth{i_shank};
                        siteDepth = firstSiteDepth - ((i_site-1)*100);
                        ylabel(sprintf('%d um', siteDepth));
                        
                        set(gca, 'yticklabel', {}, 'xticklabel', {}, 'tickdir', 'out')
                    end
                end
            end
        end
        
        % standardize the ylims
        if STANDARDIZE_YLIMS
            for i_col = 1:Ncols
                ymax = 0;
                for i_row = 1:Nrows
                    pltidx = sub2ind([Ncols, Nrows], i_col, i_row);
                    subplot(Nrows, Ncols, pltidx);
                    yvals = get(gca, 'ylim');
                    ymax = max([ymax, yvals(2)]);
                end
                for i_row = 1:Nrows
                    pltidx = sub2ind([Ncols, Nrows], i_col, i_row);
                    subplot(Nrows, Ncols, pltidx);
                    set(gca, 'ylim', [0 ymax])
                end
            end
        end
    end
    
    
    %
    % Plotting: PSTHs of the first pulse (avg across TF)
    %
    bins = -0.001 : BINWIDTH : 0.014;
    
    % make a tab
    htab(i_type) = uitab(hTabGroup,'Title', 'first pulse');
    hax(i_type) = axes('Parent', htab(i_type));
    Nrows = 8;
    Ncols = 4*Nrecpos;
    hBar = [];
    for i_pos = 1:Nrecpos
        for i_shank = 1:Nshanks;
            for i_site = 1:Nsites;
                
                % average data across TFs
                counts = [];
                Ns = [];
                for i_type = 1:numel(ptypes_tfonly)
                    
                     % check to make sure I don't include too high of TFs
                    if (1000/str2double(ptypes_tfonly{i_type}(end-1:end))) < BINWIDTH
                        continue
                    end
                    
                    pltDat = mu_shank.(ptypes_tfonly{i_type}){i_pos}.spiketimes{i_site, i_shank};
                    Ntrials = mu_shank.(ptypes_tfonly{i_type}){i_pos}.Ntrials{i_site, i_shank};
                    
                    if isempty(pltDat);
                        pltDat = inf;
                    end
                    
                    tmp_counts= histc(pltDat, bins);
                    counts = cat(1, counts, tmp_counts);
                    Ns = cat(1, Ns, Ntrials);
                end
                
                counts = sum(counts, 1) ./ sum(Ns); % in counts/bin
                counts = counts ./ BINWIDTH;
                
                plotCol = ((i_pos - 1) .* 4) + i_shank;
                plotRow = i_site;
                pltidx = sub2ind([Ncols, Nrows], plotCol, plotRow); % transpose row and col to get the correct linear index for plots
                hBar(end+1) = subplot(Nrows, Ncols, pltidx, 'align');
                bar(bins*1000, counts, 0.8, 'histc');
                axis tight
                if i_site == 1
                    title(sprintf('area: %s', mu_shank_info{i_pos}.brainAreas{i_shank}))
                end
                set(gca, 'yticklabel', {})
                if DISPLAYDEPTH
                    % figure out the depth of the etrode
                    firstSiteDepth = mu_shank_info{i_pos}.firstSiteDepth{i_shank};
                    siteDepth = firstSiteDepth - ((i_site-1)*100);
                    ylabel(sprintf('%d um', siteDepth));
                    
                    set(gca, 'yticklabel', {}, 'xticklabel', {}, 'tickdir', 'out')
                end
            end
        end
    end
    
    % standardize the ylims
    if STANDARDIZE_YLIMS
        for i_col = 1:Ncols
            ymax = 0;
            for i_row = 1:Nrows
                pltidx = sub2ind([Ncols, Nrows], i_col, i_row);
                subplot(Nrows, Ncols, pltidx);
                yvals = get(gca, 'ylim');
                ymax = max([ymax, yvals(2)]);
            end
            for i_row = 1:Nrows
                pltidx = sub2ind([Ncols, Nrows], i_col, i_row);
                subplot(Nrows, Ncols, pltidx);
                set(gca, 'ylim', [0 ymax])
            end
        end
    end
    
end % loop over mice



%% PSTH FOR EACH PULSE

if ~exist('dat', 'var')
    fin
    load 'invivo_dat_160205.mat'
end


BINWIDTH = 0.250e-3;
DELETEARTIFACTS = false; % deletes spike times in a specified window immediately after the optostim
TOOSOONTIME = 1e-3;
FILTERPSTH = false;
STANDARDIZE_YLIMS = true;
DISPLAYDEPTH = false;


% design a filter
A_smooth = 1;
B_smooth =  normpdf(-3:3, 0, 1);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load the metadata
%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~, ~, wb_info] = xlsread('invivoOptostimMetadata.xlsx', 1);
wbidx.mouse_name = strcmpi(wb_info(1,:), 'mouse name');
wbidx.site = strcmpi(wb_info(1,:), 'site');
wbidx.analyze = strcmpi(wb_info(1,:), 'include in analysis');
wbidx.dat_fname = strcmpi(wb_info(1,:), 'data file');
wbidx.shank_pos = strcmpi(wb_info(1,:), 'distance to pia for first contact per shank (um)');
wbidx.shank_area = strcmpi(wb_info(1,:), 'brain area for each shank');

% fix the shank related idx so that all four shanks get analyzed
wbidx.shank_pos(find(wbidx.shank_pos):find(wbidx.shank_pos)+3) = true;
wbidx.shank_area(find(wbidx.shank_area):find(wbidx.shank_area)+3) = true;

% delete the headers
wb_info(1:2,:) = [];

uniqueMice = unique(wb_info(:, wbidx.mouse_name));

% rearrange the data according to position on the electrode array
trodeMap = etrodemap('a32-4x8');  % pull in the map of electrode positions
[Nsites, Nshanks] = size(trodeMap);


for i_mouse = 2;%1:numel(uniqueMice)
    
    % determine which data files correspond to this mouse
    l_mouse = strcmpi(uniqueMice{i_mouse}, wb_info(:, wbidx.mouse_name));
    l_valid = logical(cell2mat(wb_info(:, wbidx.analyze)));
    l_mouse = find(l_mouse & l_valid);
    Nrecpos = numel(l_mouse);
    
    
    mu_shank = [];
    mu_shank_info = {};
    for i_pos = 1:Nrecpos;
        
        idx = l_mouse(i_pos);
        tmp_spikeTimes = dat{idx}.spikeTimes;
        ptypes = fieldnames(tmp_spikeTimes);
        i_bad = strcmpi(ptypes, 'preTime') | strcmpi(ptypes, 'postTime') | strcmpi(ptypes, 'pulseOnT_sec');
        ptypes(i_bad) = [];
        tfidx = cellfun(@(x,y) regexpi(x, '_tf'), ptypes, 'uniformoutput', false);
        ptypes_tfonly = cellfun(@(x,y) x(y+1:end), ptypes, tfidx, 'uniformoutput', false);
        Nptypes = numel(ptypes);
        
        % grab the meta-data for this recording location
        mu_shank_info{i_pos}.firstSiteDepth = wb_info(idx, wbidx.shank_pos);
        mu_shank_info{i_pos}.brainAreas = wb_info(idx, wbidx.shank_area);
        
        for i_type = 1:Nptypes
            
            for i_ch = 1:numel(tmp_spikeTimes.(ptypes{i_type}).tt)
                
                % grab the spike times for all trials
                type_spiketimes = tmp_spikeTimes.(ptypes{i_type}).tt{i_ch};
                
                % grab the Number of trials for this condition
                tmp_N = tmp_spikeTimes.(ptypes{i_type}).Ntrials{i_ch};
                
                % iterate over pulses, and make a PSTH for each pulse
                pulseOn_sec = tmp_spikeTimes.pulseOnT_sec.(ptypes{i_type});
                Npulses = numel(pulseOn_sec);
                preTime_pulse = 0.002;
                postTime_pulse = 0.014;
                bins = -preTime_pulse : BINWIDTH : postTime_pulse;
                pulse_psth = nan(Npulses, numel(bins));
                for i_pulse = 1:Npulses
                    
                    % adjust spike times to center them around the Nth
                    % pulse (i_pulse). Then compute the counts per bin.
                    adj_spiketimes = type_spiketimes-pulseOn_sec(i_pulse);
                    if DELETEARTIFACTS
                        l_toosoon = (adj_spiketimes <= TOOSOONTIME) & (adj_spiketimes > 0);
                        adj_spiketimes(l_toosoon) = [];
                    end
                    tmp_counts = histc(adj_spiketimes, bins);
                    tmp_Hz = (tmp_counts./tmp_N) ./ BINWIDTH;
                    
                    % store the data if they exist, otherwise, leave the
                    % pre-allocated nans.
                    if ~isempty(tmp_Hz)
                        if FILTERPSTH
                            tmp_Hz = filtfilt(B_smooth, A_smooth, tmp_Hz);
                        end
                        pulse_psth(i_pulse, 1:numel(bins)) = tmp_Hz;
                    end
                    
                end
                
                % put the data into the appropriate position
                [site, shank] = find(trodeMap == i_ch);
                mu_shank.(ptypes_tfonly{i_type}){i_pos}.psth_by_pulse{site, shank} = pulse_psth;
                mu_shank.(ptypes_tfonly{i_type}){i_pos}.bins_by_pulse{site, shank} = bins;
                
            end
            
        end
        
    end
    
    
     
    %
    % Plotting: PSTHs
    %
    hFig = figure;
    hFig.Units = 'normalized';
    hFig.Position = [0 0 1 1];
    hTabGroup = uitabgroup(hFig);
    for i_type = 1:numel(ptypes_tfonly)
        

        % make a tab
        htab(i_type) = uitab(hTabGroup,'Title',ptypes_tfonly{i_type});
        hax(i_type) = axes('Parent', htab(i_type));
        Nrows = 8;
        Ncols = 4*Nrecpos;
        for i_pos = 1:Nrecpos
            for i_shank = 1:Nshanks;
                for i_site = 1:Nsites;
                    
                    pltDat = mu_shank.(ptypes_tfonly{i_type}){i_pos}.psth_by_pulse{i_site, i_shank};
                    tt = mu_shank.(ptypes_tfonly{i_type}){i_pos}.bins_by_pulse{i_site, i_shank};
                    
                    Npulses = size(pltDat, 1);
                    cmap = copper(Npulses);
                    
                    plotCol = ((i_pos - 1) .* 4) + i_shank;
                    plotRow = i_site;
                    pltidx = sub2ind([Ncols, Nrows], plotCol, plotRow); % transpose row and col to get the correct linear index for plots
                    subplot(Nrows, Ncols, pltidx);
                    set(gca, 'colororder', cmap, 'NextPlot', 'replacechildren')
                    plot(bins*1000, pltDat');
                    axis tight
                    if i_site == 1
                        title(sprintf('area: %s', mu_shank_info{i_pos}.brainAreas{i_shank}))
                    end
                    
                    if DISPLAYDEPTH
                        % figure out the depth of the etrode
                        firstSiteDepth = mu_shank_info{i_pos}.firstSiteDepth{i_shank};
                        siteDepth = firstSiteDepth - ((i_site-1)*100);
                        ylabel(sprintf('%d um', siteDepth));
                        
                        set(gca, 'yticklabel', {}, 'xticklabel', {}, 'tickdir', 'out')
                    end
                end
            end
        end
        
        % standardize the ylims
        if STANDARDIZE_YLIMS
            for i_col = 1:Ncols
                ymax = 0;
                for i_row = 1:Nrows
                    pltidx = sub2ind([Ncols, Nrows], i_col, i_row);
                    subplot(Nrows, Ncols, pltidx);
                    yvals = get(gca, 'ylim');
                    ymax = max([ymax, yvals(2)]);
                end
                for i_row = 1:Nrows
                    pltidx = sub2ind([Ncols, Nrows], i_col, i_row);
                    subplot(Nrows, Ncols, pltidx);
                    set(gca, 'ylim', [0 ymax])
                end
            end
        end
        
    end
    drawnow
    
end % loop over mice


%% NOTES FOR BATCH PROCESSING

% //
% // need to know which pixels to analyze in the kCSD plots.
% 
%  * store the [x,y] location in pix or um units?
%


% current plan:
%
% Save downsampled versions of:
%   * lfpSnips (currently this gets turned into other variables for CSD and such)
%   * mu_snips



