%% NOTES

% Does Chronos cause artifical depression? Look for SOM cells that should
% evoke facilitation. Look for any instances of good facilitation using
% Chronos.

% Does oChIEF cause wonky STP? Strange facilitation followed by strong
% depression followed by rebound

% Do the oChIEF fits improve if I add an extra facilitation term?

% Do the cross validation trials work better if I use the recovery trains
% to train the model and then predict the RITs

% Should I use only a few RITs (maybe three) but present them multiple
% times?


% Try to look a the data in terms of cell type. Add DC current steps to the
% QC plots

% Are the poor fits associated with noisy Rs or 1st pulse amp? Delete the
% bad pulses and try again... make sure to save the original workbook in
% case it's not useful.

% read over code to make sure fitting and cross validation is working.



%% SPECIFY WHICH EXPERIMENTS SHOULD CONTRIBUTE, LOAD THE DATA

fin

%%%% DEFINE THE ANALYSIS PARAMS %%%%

params.pretime.vclamp = 0.002;     % seconds before pulse onset
params.posttime.vclamp = 0.015;    % seconds after pulse onset
params.pretime.dcsteps = 0.100;    % seconds before current onset
params.posttime.dcsteps = 0.300;   % seconds after current offset 
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% load in the workbook that contains all the experimental data
wb_path = [GL_DOCUPATH, 'Other_workbooks', filesep, 'wholeCellSTPCellList.xlsx'];
[~,~,wb_expt] = xlsread(wb_path, 1);

% generate the header index info
for i_atrib = 1:size(wb_expt,2)
    fldname = wb_expt{1,i_atrib};
    fldname(isspace(fldname)) = [];
    hidx.(fldname) = i_atrib;
end
    
% now that the header is formed, delete the first row.
wb_expt(1,:) = [];

% convert the file names into fully qualified paths so that par-for can run
% without calling a global
iclamp_fpath = cellfun(@(x,y) strcat(GL_DATPATH, x, filesep, 'Physiology', filesep, y, '.abf'), wb_expt(:,hidx.MouseName), wb_expt(:, hidx.ABFDCsteps), 'uniformoutput', false);
vclamp_fpath = cellfun(@(x,y) strcat(GL_DATPATH, x, filesep, 'Physiology', filesep, y, '.abf'), wb_expt(:,hidx.MouseName), wb_expt(:, hidx.ABFOptostimVclamp), 'uniformoutput', false);
wb_expt(:,hidx.ABFDCsteps) = iclamp_fpath;
wb_expt(:,hidx.ABFOptostimVclamp) = vclamp_fpath;                 
                   
% make each row it's own cell array so that it can be passed as a single
% argument to the function that does the major unpacking
attributes = {};
for i_ex = 1:size(wb_expt,1)
    attributes{i_ex,1} = wb_expt(i_ex,:);
end


%
% LOAD THE DATA FILES
%
dat = {};
Nexpts = numel(attributes);

pool = gcp('nocreate');
if isempty(pool)
    pool = parpool(21);
end

parfor i_ex = 1:Nexpts
    dat{i_ex} = wcstp_compile_data(attributes{i_ex}, hidx, params);
end

fprintf('All done importing data\n')

%% QULAITY CONTROL PLOTS

close all

for i_ex = 1:numel(dat)

    f = figure;
    f.Name = sprintf('Mouse %s, site %s', dat{i_ex}.info.mouseName, dat{i_ex}.info.siteNum);
    f.Position = [332  96 1259 665];
    
    for i_ch = 1:2
        
        fldname = sprintf('HS%dvalid', i_ch);
        isvalid = str2double(wb_expt(i_ex, hidx.(fldname)));
        if ~isvalid
            continue
        end
        
        % plot the current step data set to help identify cell types
        if isfield(dat{i_ex}, 'dcsteps')
            Vm = dat{i_ex}.dcsteps.Vm_raw{i_ch};
            Icmd = dat{i_ex}.dcsteps.Icmd{i_ch};
            N = size(Vm,2);
            tt = [0:N-1] ./ dat{i_ex}.info.sampRate.iclamp;
            if i_ch == 1; col = 1; else col = 3; end
            
            if any(Icmd<0)
                pltidx = sub2ind([4,3], col, 3);
                ha = subplot(3,4,pltidx);
                plot(tt, Vm(Icmd<0,:)');
                axis tight
                ha.Box = 'off';
                ha.TickDir = 'out';
                xlabel('Seconds')
                ylabel('mV')
            end
            
            vmidx = find(Icmd>0, 1, 'first');
            if any(vmidx);
                stopidx = min([vmidx+3, size(Vm,1)]);
                pltidx = sub2ind([4,3], col, 2);
                ha = subplot(3,4,pltidx);
                plot(tt, Vm(vmidx:stopidx,:)');
                axis tight
                ha.Box = 'off';
                ha.TickDir = 'out';
                ylabel('mV')
            end
            
            if stopidx < size(Vm,1)
                pltidx = sub2ind([4,3], col, 1);
                ha = subplot(3,4,pltidx);
                plot(tt, Vm(stopidx+1:end,:)');
                axis tight
                ha.Box = 'off';
                ha.TickDir = 'out';
                ylabel('mV')
            end
        end
        
        % figure out the subplot column number for the following axes
        if i_ch == 1; col = 2; else col = 4; end
        
        % series resistance
        if ~all(isnan(dat{i_ex}.qc.Rs{i_ch}))
            pltidx = sub2ind([4,3], col, 1);
            subplot(3,4,pltidx)
            tmp = squeeze(dat{i_ex}.qc.Rs{i_ch});
            plot(tmp)
            ylim([min([0,min(tmp(:))]), max(tmp(:))*1.5])
            ylabel('R_{s} (MOhm)')
            xlabel('trial number')
            title(sprintf('Channel %d', i_ch))
        end
        
        % verr
        if ~all(isnan(dat{i_ex}.qc.verr{i_ch}))            
            pltidx = sub2ind([4,3],col,2);
            subplot(3,4,pltidx)
            tmp = squeeze(dat{i_ex}.qc.verr{i_ch});
            plot(tmp)
            ylim([min([0,min(tmp(:))]), max(tmp(:))*1.5])
            ylabel('SS Verr (mV)')
            xlabel('trial number')
            
        end
        
        % p1amps
        if ~all(isnan(dat{i_ex}.qc.p1amp{i_ch}))
            pltidx = sub2ind([4,3], col, 3);
            subplot(3,4,pltidx), hold on,
            tmp = squeeze(dat{i_ex}.qc.p1amp{i_ch});
            tmp_norm = squeeze(dat{i_ex}.qc.p1amp_norm{i_ch});
            plot(tmp)
            plot(tmp_norm, 'r')
            ylim([min([0,min(tmp(:))]), max(tmp(:))*1.5])
            ylabel('P1 Amp')
            xlabel('trial number')
            %set(gca, 'yscale', 'log', 'ylim', [0.3333, 3])
        end
            
        
    end
    drawnow
    
end



%% STP SUMMARY FOR EACH RECORDING

close all
NORMAMPS = true;
PLOTRAW = true;

for i_ex = 1:numel(dat)
    
    f = figure;
    f.Name = sprintf('Mouse %s, site %s', dat{i_ex}.info.mouseName, dat{i_ex}.info.siteNum);
    f.Units = 'normalized';
    f.Position = [0.1363    0.0301    0.7535    0.8762];
    
    for i_ch = 1:2
        
        conds = fieldnames(dat{i_ex}.expt);
        Nconds = numel(conds);
        for i_cond = 1:Nconds
            if isempty(dat{i_ex}.expt.(conds{i_cond}).stats.EPSCamp{i_ch})
                continue
            end
            
            subplot(Nconds,2, 2.*(i_cond-1) + i_ch);
            
            tmp = dat{i_ex}.expt.(conds{i_cond}).stats.EPSCamp{i_ch};
            
            if NORMAMPS
                p1 = squeeze(dat{i_ex}.qc.p1amp{i_ch});
                normfact = squeeze(dat{i_ex}.qc.p1amp_norm{i_ch});
                realTrlNums = dat{i_ex}.expt.(conds{i_cond}).realTrialNum{i_ch};
                normfact = normfact(realTrlNums);
                normfact = permute(normfact, [3,1,2]);
                tmp = bsxfun(@rdivide, tmp, normfact);
            end
            
            tt = dat{i_ex}.expt.(conds{i_cond}).pOnTimes;
            if PLOTRAW
                plot(tt, permute(tmp, [1,3,2]), '-')
            else
                xbar = mean(tmp,3);
                sem = stderr(tmp,3);
                my_errorbar(tt, xbar, sem, 'ok', 'markersize', 2, 'linewidth', 1);
            end
            
            xlim([0, dat{i_ex}.info.sweepLength.vclamp]);
            if i_cond < Nconds
                set(gca, 'xticklabel', []);
            end
             
            if ~NORMAMPS
                set(gca, 'yticklabel', []);
            end
            
        end
    end
    drawnow
    
end




%% ESTIMATE TIME CONSTANTS OF SHORT TERM PLASTICITY

TRAINSET = 'rit';  % could be 'rit', 'recovery'
PLOTTRAININGDATA = true;
NORMALIZEDATA = true;
FITAVERAGEDATA = true;


% write an anyonomous helper function to find the training data
clear isTrainingSet
switch TRAINSET
    case 'rit'
        isTrainingSet = @(x) strncmp(x, 'RITv', 4);
    case 'recovery'
        isTrainingSet = @(x) ~strncmp(x, 'RITv', 4);
end

% aggregate some things across experiments
pop.R2_raw = {};
pop.R2_crossval = {};
pop.fitparams = {};

for i_ex = 1:numel(dat)
   clc
   hf = figure;
   hf.Position = [290 -66 1325 825];
   hf.Name = sprintf('Mouse %s, site %s, opsin: %s.  Train with: %s, Plot training set: %d',...
                      dat{i_ex}.info.mouseName, dat{i_ex}.info.siteNum, dat{i_ex}.info.opsin, TRAINSET, PLOTTRAININGDATA);
   chempty = false(1,2);
   for i_ch = 1:2
       
       % determine if there are data to fit
       fldname = sprintf('HS%dvalid', i_ch);
       isvalid = str2double(wb_expt(i_ex, hidx.(fldname)));
       if ~isvalid
           chempty(i_ch) = true;
           if i_ch == 1
               continue
           elseif all(chempty)
               close(hf)
               continue
           end
       end
       
       % make a mini dataset that's composed of only pOnTimes, EPSC amps,
       % and p1Amps
       [pOnTimes, rawAmps, p1Amps, raw_xbar, raw_sem] = deal({});
       condnames = fieldnames(dat{i_ex}.expt);
       l_trainingSet = isTrainingSet(condnames);
       for i_cond = 1:numel(condnames)
           
           % grab pOnTimes
           pOnTimes{i_cond} = dat{i_ex}.expt.(condnames{i_cond}).pOnTimes;
           rawAmps{i_cond} = dat{i_ex}.expt.(condnames{i_cond}).stats.EPSCamp{i_ch};
           
           % normalize to the p1amp_norm if need be
           trlNums = dat{i_ex}.expt.(condnames{i_cond}).realTrialNum{i_ch};
           normfacts = dat{i_ex}.qc.p1amp_norm{i_ch}(trlNums);
           normfacts = permute(normfacts, [1,3,2]);
           if NORMALIZEDATA
               
               rawAmps{i_cond} = bsxfun(@rdivide, rawAmps{i_cond}, normfacts);
               raw_sem{i_cond} = stderr(rawAmps{i_cond}, 3);
               raw_xbar{i_cond} = mean(rawAmps{i_cond},3);
               if FITAVERAGEDATA
                   p1Amps{i_cond} = 1;
               else
                   p1Amps{i_cond} = ones(1, size(rawAmps{i_cond},3));
               end
               
           else
               
               % leave 'rawAmps' unchanged, but update these:
               raw_sem{i_cond} = stderr(rawAmps{i_cond}, 3);
               raw_xbar{i_cond} = mean(rawAmps{i_cond},3);
               if FITAVERAGEDATA
                   p1Amps{i_cond} = raw_xbar{i_cond}(1);
               else
                   p1Amps{i_cond} = squeeze(normfacts);
               end
               
           end
       end
       
       % allocate the training data
       [training_amps, training_pOnTimes, training_p1Amps] = deal({});
       if FITAVERAGEDATA
           training_amps = raw_xbar(l_trainingSet);
       else
           training_amps = rawAmps(l_trainingSet);
       end
       training_pOnTimes = pOnTimes(l_trainingSet);
       training_p1Amps = p1Amps(l_trainingSet);
       
       
       % if there were no training data, then move along,
       if isempty(training_amps)
           chempty(i_ch) = true;
           if i_ch == 1
               continue
           elseif all(chempty)
               close(hf)
               continue
           end
       end
       NtrainingSets = sum(l_trainingSet);
       
       % look for instances where there are pOnTimes but no data (could
       % happen if some sweeps get deleted from one HS but not the other.
       l_empty = cellfun(@isempty, training_amps);
       training_amps(l_empty) = [];
       training_pOnTimes(l_empty) = [];
       training_p1Amps(l_empty) = [];
       
       % fit the RIT data
       [d, f, dTau, fTau] = deal([]); %#ok<*ASGLU>
       [d, f, dTau, fTau] = fitTau2STP(training_amps, training_pOnTimes, training_p1Amps, 'multistart');
       
       % predict all the data
       pred = {};
       for i_cond = 1:numel(condnames)
           if ~isempty(rawAmps{i_cond})
               A0 = mean(p1Amps{i_cond});
               pred{i_cond} = predictPSCfromTau(pOnTimes{i_cond}, d, dTau, f, fTau, A0);
           else
               pred{i_cond} = [];
           end
       end
       
       
       % plot the cross validation data set, and the prediction to the
       % corss validation
       figure(hf)
       if PLOTTRAININGDATA
           l_condsToPlot = isTrainingSet(condnames);
       else
           l_condsToPlot = ~isTrainingSet(condnames);
       end
       idx_condsToPlot = find(l_condsToPlot);
       xlims = [inf -inf];
       ylims = [inf -inf];
       hs = [];
       for i_cond = 1:numel(idx_condsToPlot)
           
           typeIdx = idx_condsToPlot(i_cond);
           if isempty(raw_xbar{typeIdx}); continue; end
           
           pltIdx = sub2ind([4, numel(idx_condsToPlot)], i_ch+1, i_cond);
           hs(i_cond) = subplot(numel(idx_condsToPlot), 4, pltIdx); hold on,
           
           xx = pOnTimes{typeIdx};
           my_errorbar(xx, raw_xbar{typeIdx}, raw_sem{typeIdx}, 'k');
           plot(xx, pred{typeIdx}, 'r', 'linewidth', 2)
           xlims(1) = min([min(xx), xlims(1)]);
           xlims(2) = max([max(xx), xlims(2)]);
           yvals = get(gca, 'ylim');
           ylims(1) = min([yvals(1), ylims(1)]);
           ylims(2) = max([yvals(2), ylims(2)]);
       end
       set(hs, 'XLim', xlims, 'YLim', ylims)
       
       
       
       % make a scatter plot of all predicted and actual PSC amps
       hs = [];
       training_raw = [];
       training_pred = [];
       crossval_raw = [];
       crossval_pred = [];
       for i_cond = 1:numel(condnames)
           if isTrainingSet(condnames{i_cond}); pltclr = 'k';else pltclr = 'r';end
           if i_ch == 1; pltcol=1; else pltcol=4; end
           
           
           
           if isempty(raw_xbar{i_cond}); continue; end
           
           tmp_raw = raw_xbar{i_cond};
           tmp_pred = pred{i_cond};
           
           tmp_raw = tmp_raw(:);
           tmp_pred = tmp_pred(:);
           assert(all(size(tmp_pred) == size(tmp_raw)))
           
           pltIdx = sub2ind([4, 3], pltcol, 1);
           hs = subplot(3, 4, pltIdx); hold on,
           plot(tmp_raw, tmp_pred, '.', 'color', pltclr)
           axis tight
           
           % concatenate the training data
           if isTrainingSet(condnames{i_cond})
               training_raw = cat(1, training_raw, tmp_raw(:));
               training_pred = cat(1, training_pred, tmp_pred(:));
           end
           
           % concatenate the cross-validation data
           if ~isTrainingSet(condnames{i_cond})
               crossval_raw = cat(1, crossval_raw, tmp_raw(:));
               crossval_pred = cat(1, crossval_pred, tmp_pred(:));
           end
       end
       maxval = max([hs.XLim, hs.YLim]);
       minval = min([hs.XLim, hs.YLim]);
       plot([minval, maxval], [minval, maxval], 'k--')
       %hs.XScale = 'log';
       %hs.YScale = 'log';
       xlabel('raw EPSC amp')
       ylabel('pred amp')
       title(sprintf('num RITs fit: %d', NtrainingSets))
       
       pltIdx = sub2ind([4, 3], pltcol, 2);
       subplot(3,4,pltIdx), hold on,
       resid = training_pred - training_raw;
       histogram(resid)
       plot(mean(resid), 10, 'rv', 'markerfacecolor', 'r')
       R2 = 1 - (sum(resid.^2) ./ sum((training_raw - mean(training_raw)).^2));
       xlabel('pred-real')
       
       % calculate the suffle corrected R2.
       N = numel(training_pred);
       iters = 5000;
       shuffle_inds = nan(N, iters);
       for i_itr = 1:iters
           shuffle_inds(:,i_itr) = randperm(N)';
       end
       shuffle_pred = training_pred(shuffle_inds);
       resid = bsxfun(@minus, shuffle_pred, training_raw);
       SS_resid = sum(resid.^2, 1);
       SS_raw = sum((training_raw-mean(training_raw)).^2);
       R2_shuffle = 1 - (SS_resid ./ SS_raw);
       title(sprintf('R2: %.3f, R2_shuff: %.3f', R2, mean(R2_shuffle)))
       
       
       % plot cross validation stuff
       if ~isempty(crossval_raw)
           pltIdx = sub2ind([4, 3], pltcol, 3);
           subplot(3,4,pltIdx), hold on,
           resid = crossval_pred - crossval_raw;
           histogram(resid);
           plot(mean(resid), 5, 'rv', 'markerfacecolor', 'r')
           R2 = 1 - (sum(resid.^2) ./ sum((crossval_raw-mean(crossval_raw)).^2));
           xlabel('cross-valid (pred-real)')
           
           N = numel(crossval_pred);
           shuffle_inds = nan(N, iters);
           for i_itr = 1:iters
               shuffle_inds(:,i_itr) = randperm(N)';
           end
           shuffle_pred = crossval_pred(shuffle_inds);
           resid = bsxfun(@minus, shuffle_pred, crossval_raw);
           SS_resid = sum(resid.^2, 1);
           SS_raw = sum((crossval_raw-mean(crossval_raw)).^2);
           R2_shuffle = 1 - (SS_resid ./ SS_raw);
           title(sprintf('R2: %.3f, R2_shuff: %.3f', R2, mean(R2_shuffle)))
       end
       
   end
   drawnow
   
   
end


%% PAIRED PULSE PLASTICITY MANIFOLDS (DATA COLLECTION)

% loop through the experiments. Pull out the trains data. Ignore the
% recovery train (if present) and aggregate across recovery conditions.
pprpop = [];
pprpop.TFsAllExpts = [];
pprpop.MaxNPulses = 0;
for i_ex = 1:numel(dat)
    
    % find the normal trains. Assume the field name is NOT 'ritv'
    condnames = fieldnames(dat{i_ex}.expt);
    l_trains = ~strncmp(condnames, 'RITv', 4);
    if sum(l_trains)==0; continue; end % no trains data
    trainParams = cellfun(@(x) dat{i_ex}.expt.(condnames{x}).tdict, mat2cell(find(l_trains), ones(sum(l_trains),1), 1), 'uniformoutput', false);
    trainParams = cat(1, trainParams{:});
    
    % make sure the pulse amplitude and width were identical across ttypes
    assert(numel(unique(trainParams(:,1)))==1, 'ERROR: more than one pulse amplitude')
    assert(numel(unique(trainParams(:,2)))==1, 'ERROR: more than one pulse width');
    
    % identify the unique TF conditions for this experiment, and update the
    % running log of TFs used across all experiments
    uniqueTFs = unique(trainParams(:,3));
    tmp = cat(1, pprpop.TFsAllExpts, uniqueTFs);
    pprpop.TFsAllExpts = unique(tmp);
    
    % store the stim params for all stim types, which will be useful for
    % indexing later.
    allStimParams = cellfun(@(x) dat{i_ex}.expt.(condnames{x}).tdict, mat2cell((1:numel(l_trains))', ones(numel(l_trains),1), 1), 'uniformoutput', false);
    allStimParams = cat(1, allStimParams{:});
    
    
    % aggregate data across within TF conditons (Separately for each
    % recording channel
    for i_tf = 1:numel(uniqueTFs);
        
        tfidx = find(trainParams(:,3) == uniqueTFs(i_tf)); % condnames that contain the train with a particular TF
        
        for i_ch = 1:2
            
            % check to make sure this neuron was defined
            fldname = sprintf('HS%dvalid', i_ch);
            isvalid = str2double(wb_expt(i_ex, hidx.(fldname)));
            if ~isvalid
                pprpop.dat{i_ex}.xbar{i_tf}{i_ch} = [];
                pprpop.tfs{i_ex}{i_ch} = [];
                pprpop.recoverypulse{i_ex}{i_tf}{i_ch} = [];
                continue
            end
            
            % iterate over the trains with the same freq. 
            catdat = [];
            for i_cond = 1:numel(tfidx)
                
                % pull out the data
                condIdx = ismember(allStimParams, trainParams(tfidx(i_cond),:), 'rows');
                tmpdat = dat{i_ex}.expt.(condnames{condIdx}).stats.EPSCamp{i_ch}; % [Npulses, 1, Nsweeps]
           
                % normalize by the p1Amp_norm
                realTrlNums = dat{i_ex}.expt.(condnames{condIdx}).realTrialNum{i_ch};
                p1Amp_norm = dat{i_ex}.qc.p1amp_norm{i_ch}(realTrlNums);
                assert(~any(isnan(p1Amp_norm)), 'ERROR: scale factor is a nan');
                p1Amp_norm = permute(p1Amp_norm, [1,3,2]);
                %tmpdat = bsxfun(@rdivide, tmpdat, p1Amp_norm); % [Npulses, 1, Nsweeps]
                
                % store in a matrix for averaging later
                catdat = cat(3, catdat, tmpdat);
                
            end
            
            % check to make sure there are data for these conditions. Even
            % though this recording channel should be defined (See above),
            % it's possible there are no data due to deletion of single
            % sweeps
            if isempty(catdat)
                pprpop.dat{i_ex}.xbar{i_tf}{i_ch} = [];
                pprpop.tfs{i_ex}{i_ch} = [];
                pprpop.recoverypulse{i_ex}{i_tf}{i_ch} = [];
            else
                
                % now average across trails, and re-normalize to the first
                % pulse. Store in the population data structure.
                avg = mean(catdat,3);
                ppr = avg./avg(1);
                
                % delete the last pulse (if it's a recovery pulse, but
                % store the recovery pulse in a different field of the
                % population structure
                isrecovery = trainParams(tfidx(i_cond),4) > 0;
                if isrecovery
                    pprpop.recoverypulse{i_ex}{i_tf}{i_ch} = ppr(end);
                    ppr(end) = [];
                end
                
                pprpop.dat{i_ex}.xbar{i_tf}{i_ch} = ppr;
                pprpop.tfs{i_ex}{i_ch} = uniqueTFs;
                pprpop.MaxNPulses = max([pprpop.MaxNPulses, numel(ppr)]);
                
            end
            
        end
        
    end
    
end


%% PAIRED PULSE PLASTICITY MANIFOLDS (PLOTS)

% define a set of attributes for each lineseries (or manifold) in the plot
% {CellType,  BrainArea,  OpsinType}
plotgroups = {
              'PY',    'al', 'chronos';...
              'PY',    'al', 'chief';...
              };

groupdata = repmat({[]}, 1, size(plotgroups, 1)); % should only have N cells, where N = size(plotgroups, 1). Each cell has a matrix with a cononicalGrid:
canonicalGrid = nan(pprpop.MaxNPulses, numel(pprpop.TFsAllExpts)); % need to get the maxNpulses into the pprpop struct
allTFs = pprpop.TFsAllExpts;

% iterate over the experiments. For each recording channel, determine what
% the attributes are, and place the data in the correct ploting group.
for i_ex = 1:numel(dat)
    for i_ch = 1:2
        % check to make sure this neuron was defined
        HSname = sprintf('HS%dvalid', i_ch);
        isvalid = str2double(wb_expt(i_ex, hidx.(HSname)));
        if ~isvalid
            continue
        end
        
        % check the attributes
        wbcol_CellType = hidx.(sprintf('Celltype%d', i_ch));
        wbcol_BrainArea = hidx.brainarea;
        wbcol_OpsinType = hidx.opsin;
        ch_attribs = wb_expt(i_ex, [wbcol_CellType, wbcol_BrainArea, wbcol_OpsinType]);
        l_nan = cellfun(@(x) all(isnan(x)), ch_attribs);
        ch_attribs(l_nan) = cellfun(@num2str, ch_attribs(l_nan), 'uniformoutput', false);
        
        l_cellType_match = cellfun(@(x) ~isempty(regexpi(ch_attribs{1}, x)), plotgroups(:,1)) | strcmpi(plotgroups(:,1), 'any');
        l_brainArea_match = cellfun(@(x) ~isempty(regexpi(ch_attribs{2}, x)), plotgroups(:,2)) | strcmpi(plotgroups(:,2), 'any');
        l_opsinMatch = cellfun(@(x) ~isempty(regexpi(ch_attribs{3}, x)), plotgroups(:,3)) | strcmpi(plotgroups(:,3), 'any');
        
        group_idx = sum([l_cellType_match, l_brainArea_match, l_opsinMatch], 2) == 3;
        assert(sum(group_idx)<=1, 'ERROR: found too many group indicies')
        if sum(group_idx) == 0; continue; end
        
        % add data to the appropriate group data array
        tmpgrid = canonicalGrid;
        ch_tfs = pprpop.tfs{i_ex}{i_ch};
        for i_tf = 1:numel(ch_tfs)
            grid_col_idx = allTFs == ch_tfs(i_tf);
            tmpdat = pprpop.dat{i_ex}.xbar{i_tf}{i_ch};
            Npulses = numel(tmpdat);
            tmpgrid(1:Npulses, grid_col_idx) = tmpdat(:);
        end
        
        if all(isnan(tmpgrid(:))); continue; end
        
        % aggregate the data
        groupdata{group_idx} = cat(3, groupdata{group_idx}, tmpgrid);
        
    end
end


% Plot average manfold. any global (above) that has multiple elements will
% be ploted against eachother, but the other globals will be held constant.
hf = figure;
hold on,
plotcolors = {'r', 'b', 'g'}
for i_group = 1:numel(groupdata)
    
    
    grid_average = nanmean(groupdata{i_group},3);
    grid_N = sum(~isnan(groupdata{i_group}),3)
    
    Y = 1:size(groupdata{i_group},1);
    X = 1:numel(allTFs);
    
    l_nan_tfs = all(isnan(grid_average), 1);
    
    hs = surf(X(:,~l_nan_tfs), Y, flipud(grid_average(:,~l_nan_tfs)));
    hs.EdgeColor = 'k';
    hs.EdgeAlpha = 0.5;
    hs.FaceColor = plotcolors{i_group};
    hs.FaceAlpha = 1;
    
end
set(gca, 'zscale', 'log', 'view', [-43    16])
set(gca, 'ytick', Y, 'yticklabel', cellfun(@(x) num2str(x), num2cell(rot90(Y,2)), 'uniformoutput', false))
set(gca, 'xtick', X, 'xticklabel', cellfun(@(x) num2str(x), num2cell(allTFs), 'uniformoutput', false))
xlabel('Temporal Frequency')
ylabel('Pulse Number')
zlabel('norm amp')




%% POST-TETANIC POTENTIATION

% plot recovery pulse amplitude (normalized) vs. train frequency. Do this
% separately for cell types, brain areas, opsins, etc...



