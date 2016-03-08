%% designate an experiment to analyze, import the data
fin

profile on

NSX = 'ns5'; % the version with the LFP continuous data
STIMTYPE = 'train'; % train or sinusoid
NOISEMETHOD = 'filter'; % 'subtract' or 'filter'


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
for i_ex = 1:Nexpts;
    
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
    % divided into sampFreq_nsx. The lp_filter will be sampFreq_lfp/5
    sampFreq_lfp = 1250;
    
    pool = gcp('nocreate');
    if isempty(pool)
        pool = parpool(32);
    end
  
    dat{i_ex} = extractLFPandMU(blk, lfp_ch_idx, sampFreq_lfp, NOISEMETHOD, STIMTYPE, NSX);
    
end


%% PLOT: AVERAGE LFP SIGNAL

PLOTTYPE = 'errbar'; % 'errbar', 'all', 'mean'

% make a tmp variable (duplicate of the preprocessed LFP data)
trialSnips = lfpSnips;

% pull in the map of electrode positions
trodeMap = etrodemap('a32-4x8');
trodePltIdx = reshape(trodeMap', [], 1);


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

METHOD = 'kcsd'; % could be 'normal', or 'kcsd'

% make a tmp variable (duplicate of the preprocessed LFP data)
csdSnips = lfpSnips;



ptypes = fieldnames(csdSnips);
i_bad = strcmpi(ptypes, 'preTime') | strcmpi(ptypes, 'postTime');
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
                
                elPos = 0.050:0.100:0.750;
                X = 0:0.01:0.900;
                pots = tmp_data;
                
                k = kCSD1d(elPos, pots, 'X', X);
                k.estimate;
                csd_final.(ptypes{i_cond}){i_shank} = k.csdEst;
        end
        
    end
end

% Try to compute the 2D CSD
if strcmpi(METHOD, 'kcsd');
    elPos_y = repmat([0:0.100:.700]', 1, 4);
    elPos_y = elPos_y(:);
    elPos_x = repmat([0.100:0.400:1.30], 8, 1);
    elPos_x = elPos_x(:);
    elPos = [elPos_x, elPos_y];
    for i_cond = 1:numel(ptypes)
        pots = [];
        for i_shank = 1:4
            pots = cat(1, pots, lfp_final.(ptypes{i_cond}){i_shank});
        end
        
        [X, Y] = meshgrid([0:0.010:1.4], [0:0.010:0.750]);
        
        k = kcsd2d(elPos, pots, 'X', X, 'Y', Y);
        k.plot_CSD;
        
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
        plotDat = csd_final.(ptypes{i_cond}){i_shank};
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
    Nsamps = ceil(0.015 .* sampFreq_lfp);
    idx = [t_zero_idx - Nsamps : t_zero_idx + Nsamps];
    
    % concatenate across TF conditions
    plotDat = [];
    for i_cond = 1:numel(ptypes)
        plotDat = cat(3, plotDat, csd_final.(ptypes{i_cond}){i_shank}(:,idx));
    end
    
    % average across TFs
    plotDat = mean(plotDat, 3);
    
    % add row/col that will be deleted by pcolor
    plotDat = [plotDat(1,:); plotDat];
    plotDat = [plotDat, plotDat(:,1)];
    
    
    % do the plotting
    s=subplot(4,1,i_shank);
    newN = size(plotDat,2);
    tt_sec = ((0:newN-1)-Nsamps) ./ sampFreq_lfp;
    pos = linspace(0,8,size(plotDat, 1));
    p=pcolor(tt_sec, pos, flipud(plotDat));
    shading interp;
    s.YTickLabel = flipud(s.YTickLabel);
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


for i_mouse = 1:numel(uniqueMice)
    
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
        i_bad = strcmpi(ptypes, 'preTime') | strcmpi(ptypes, 'postTime');
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
                
                % put the data into the appropriate position
                tmp_spiketimes = tmp_spikeTimes.(ptypes{i_type}).tt{i_ch};
                tmp_N = tmp_spikeTimes.(ptypes{i_type}).Ntrials{i_ch};
                mu_shank.(ptypes_tfonly{i_type}){i_pos}.spiketimes{site, shank} = tmp_spiketimes;
                mu_shank.(ptypes_tfonly{i_type}){i_pos}.Ntrials{site, shank} = tmp_N;
                
            end
            
        end
        
    end
    
    
    %
    % Plotting: PSTHs
    %
    BINWIDTH = 0.000500;
    STANDARDIZE_YLIMS = false;
    startTime = -dat{idx}.spikeTimes.preTime;
    hFig = figure;
    hFig.Units = 'normalized';
    hFig.Position = [0 0 1 1];
    hTabGroup = uitabgroup(hFig);
    for i_type = 1:numel(ptypes_tfonly)
        
        l_noData = cellfun(@isempty, dat{idx}.lfpsnips.(ptypes{i_type}));
        Ntimepoints = cellfun(@(x) size(x, 2), dat{idx}.lfpsnips.(ptypes{i_type}));
        Ntimepoints = unique(Ntimepoints(~l_noData));
        sampFreq_nsx = dat{idx}.info.sampFreq_lfp;
        totalTime = (Ntimepoints ./ sampFreq_nsx);
        endTime = totalTime + startTime;
        bins = startTime : BINWIDTH : endTime;
        
        assert(numel(Ntimepoints) == 1, 'ERROR: time is not consistent')
        
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
                        continue
                    end
                    
                    counts= histc(pltDat, bins);
                    counts = (counts./Ntrials) ./ BINWIDTH;
                    
                    plotCol = ((i_pos - 1) .* 4) + i_shank;
                    plotRow = i_site;
                    pltidx = sub2ind([Ncols, Nrows], plotCol, plotRow); % transpose row and col to get the correct linear index for plots
                    hBar(end+1) = subplot(Nrows, Ncols, pltidx, 'align');
                    bar(bins*1000, counts, 0.8, 'histc');
                    axis tight
                    if i_site == 1
                        title(sprintf('area: %s', mu_shank_info{i_pos}.brainAreas{i_shank}))
                    end
                end
            end
        end
        
        % standardize the ylims
        if STANDARDIZE_YLIMS
            ymax = 0;
            for i_h = 1:numel(hBar)
                yvals = get(hBar(i_h), 'ylim');
                ymax = max([ymax, yvals(2)]);
            end
            for i_h = 1:numel(hBar)
                set(hBar(i_h), 'ylim', [0 ymax])
            end
        end
        
    end
    
    %
    % Plotting: PSTHs of the first pulse (avg across TF)
    %
    BINWIDTH = 200e-6;
    STANDARDIZE_YLIMS = false;
    bins = -0.001 : BINWIDTH : 0.015;
    
    assert(numel(Ntimepoints) == 1, 'ERROR: time is not consistent')
    
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
                for i_type = 1:numel(ptypes_tfonly)
                    
                    pltDat = mu_shank.(ptypes_tfonly{i_type}){i_pos}.spiketimes{i_site, i_shank};
                    Ntrials = mu_shank.(ptypes_tfonly{i_type}){i_pos}.Ntrials{i_site, i_shank};
                    
                    if isempty(pltDat);
                        continue
                    end
                    
                    tmp_counts= histc(pltDat, bins);
                    tmp_counts = (tmp_counts./Ntrials) ./ BINWIDTH;
                    counts = cat(1, counts, tmp_counts);
                end
                
                counts = mean(counts, 1);
                
                plotCol = ((i_pos - 1) .* 4) + i_shank;
                plotRow = i_site;
                pltidx = sub2ind([Ncols, Nrows], plotCol, plotRow); % transpose row and col to get the correct linear index for plots
                hBar(end+1) = subplot(Nrows, Ncols, pltidx, 'align');
                bar(bins*1000, counts, 0.8, 'histc');
                axis tight
                if i_site == 1
                    title(sprintf('area: %s', mu_shank_info{i_pos}.brainAreas{i_shank}))
                end
                
            end
        end
    end
    % standardize the ylims
    if STANDARDIZE_YLIMS
        ymax = 0;
        for i_h = 1:numel(hBar)
            yvals = get(hBar(i_h), 'ylim');
            ymax = max([ymax, yvals(2)]);
        end
        for i_h = 1:numel(hBar)
            set(hBar(i_h), 'ylim', [0 ymax])
        end
    end
    
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



