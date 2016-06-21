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


% decide what experiment to run
EXPTTYPE = 2;
switch EXPTTYPE
    case 1
        EXPTTYPE = 'all';
    case 2
        EXPTTYPE = 'Manifold';
    case 3
        EXPTTYPE = 'RIT_test';
end



%%%% DEFINE THE ANALYSIS PARAMS %%%%

params.pretime.vclamp = 0.002;     % seconds before pulse onset
params.posttime.vclamp = 0.015;    % seconds after pulse onset
params.pretime.dcsteps = 0.100;    % seconds before current onset
params.posttime.dcsteps = 0.300;   % seconds after current offset 
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% figure out which experiments should be analyzed
wb_path = [GL_DOCUPATH, 'Other_workbooks', filesep, 'wholeCellSTPCellList.xlsx'];
[~, ~, wb_expt] = xlsread(wb_path, 2);
if ~strcmpi(EXPTTYPE, 'all')
    header_idx = strcmpi(EXPTTYPE, wb_expt(1,:));
    assert(sum(header_idx) == 1)
    l_to_analyze = cellfun(@(x) numel(x)==1 && x==1, wb_expt(:, header_idx));
else
    l_to_analyze = ones(size(wb_expt,1), 1);
end
l_to_analyze(1) = []; % delete the header row
expt_idx = find(l_to_analyze);

% load in the workbook that contains all the experimental information
[~,~,wb_info] = xlsread(wb_path, 1);

% generate the header index info
for i_atrib = 1:size(wb_info,2)
    fldname = wb_info{1,i_atrib};
    fldname(isspace(fldname)) = [];
    hidx.(fldname) = i_atrib;
end
    
% now that the header is formed, delete the first row.
wb_info(1,:) = [];

% convert the file names into fully qualified paths so that par-for can run
% without calling a global
iclamp_fpath = cellfun(@(x,y) strcat(GL_DATPATH, x, filesep, 'Physiology', filesep, y, '.abf'), wb_info(:,hidx.MouseName), wb_info(:, hidx.ABFDCsteps), 'uniformoutput', false);
vclamp_fpath = cellfun(@(x,y) strcat(GL_DATPATH, x, filesep, 'Physiology', filesep, y, '.abf'), wb_info(:,hidx.MouseName), wb_info(:, hidx.ABFOptostimVclamp), 'uniformoutput', false);
wb_info(:,hidx.ABFDCsteps) = iclamp_fpath;
wb_info(:,hidx.ABFOptostimVclamp) = vclamp_fpath;                 
                   
% make each row it's own cell array so that it can be passed as a single
% argument to the function that does the major unpacking
attributes = {};
for i_ex = 1:numel(expt_idx)
    attributes{i_ex,1} = wb_info(expt_idx(i_ex),:);
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
        
        isvalid = dat{i_ex}.info.HS_is_valid(i_ch);
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

TRAINSET = 'recovery';  % could be 'rit', 'recovery', 'all'
PLOTTRAININGDATA = true;
NORMALIZEDATA = true;
FITAVERAGEDATA = true;
FITRECOVERYPULSE = true;


% write an anyonomous helper function to find the training data
clear isTrainingSet
switch TRAINSET
    case 'rit'
        isTrainingSet = @(x) strncmp(x, 'RITv', 4);
    case 'recovery'
        isTrainingSet = @(x) ~strncmp(x, 'RITv', 4);
    case 'all'
        isTrainingSet = @(x) true(size(x));
end


for i_ex = 1:numel(dat)
   clc
   hf = figure;
   hf.Units = 'Normalized';
   hf.Position = [0.0385 0.1759 0.6630 0.7481];
   hf.Name = sprintf('Mouse %s, site %s, opsin: %s.  Train with: %s, Plot training set: %d',...
                      dat{i_ex}.info.mouseName, dat{i_ex}.info.siteNum, dat{i_ex}.info.opsin, TRAINSET, PLOTTRAININGDATA);
   chempty = false(1,2);
   for i_ch = 1:2
       
       % determine if there are data to fit
       isvalid = dat{i_ex}.info.HS_is_valid(i_ch);
       if ~isvalid
           chempty(i_ch) = true;
           if  all(chempty)
               close(hf)
           end
           continue
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
           
           % delete the recovery pulse if need be
           if ~FITRECOVERYPULSE
               has_recov_pulse = dat{i_ex}.expt.(condnames{i_cond}).tdict(4) > 0;
               if has_recov_pulse
                   pOnTimes{i_cond}(end) = [];
                   rawAmps{i_cond} = rawAmps{i_cond}(1:end-1,1,:);
               end
           end
           
           
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
       
       % fit the RIT data, but only if the fit params do not already exist
       [d, f, dTau, fTau] = deal([]); %#ok<*ASGLU>
       MAKENEWFITS = true;
       if [isfield(dat{i_ex}, 'stpfits') ...
               && isfield(dat{i_ex}.stpfits, 'modelParams') ...
               && numel(dat{i_ex}.stpfits.modelParams)>= i_ch ...
               && ~isempty(dat{i_ex}.stpfits.modelParams{i_ch})] % end of if conditional
           
           % if you've made it here, then the old params are potentially
           % still appliciable, but I should check to make sure that they
           % were fit using the same technique
           matches = [strcmpi(dat{i_ex}.stpfits.trainingSet, TRAINSET); ...
                      dat{i_ex}.stpfits.normalizeData == NORMALIZEDATA; ...
                      dat{i_ex}.stpfits.fitWithAvg == FITAVERAGEDATA; ...
                      dat{i_ex}.stpfits.fitRecovPulse == FITRECOVERYPULSE];
           
          if all(matches)
              d = dat{i_ex}.stpfits.modelParams{i_ch}(1:2);
              f = dat{i_ex}.stpfits.modelParams{i_ch}(3);
              dTau = dat{i_ex}.stpfits.modelParams{i_ch}(4:5);
              fTau = dat{i_ex}.stpfits.modelParams{i_ch}(6);

              MAKENEWFITS = false;
          end
       end
       if MAKENEWFITS
           [d, f, dTau, fTau] = fitTau2STP(training_amps, training_pOnTimes, training_p1Amps, 'multistart');
       end
       
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
       
       
       % plot the training or cross validation data set, and the prediction
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
       if ~isempty(hs) && sum(hs)>0
           set(hs(hs~=0), 'XLim', xlims, 'YLim', ylims)
       end
       
       
       
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
       R2_train = 1 - (sum(resid.^2) ./ sum((training_raw - mean(training_raw)).^2));
       xlabel('pred-real')
       title(sprintf('R2 = %.2f', R2_train));
       
       
       
       % plot cross validation stuff
       R2_crossvalid = [];
       if ~isempty(crossval_raw)
           pltIdx = sub2ind([4, 3], pltcol, 3);
           subplot(3,4,pltIdx), hold on,
           resid = crossval_pred - crossval_raw;
           histogram(resid);
           plot(mean(resid), 5, 'rv', 'markerfacecolor', 'r')
           R2_crossvalid = 1 - (sum(resid.^2) ./ sum((crossval_raw-mean(crossval_raw)).^2));
           xlabel('cross-valid (pred-real)')
           title(sprintf('R2 = %.2f', R2_crossvalid));
       end
       
       
       % store some parameters in the dat array
       dat{i_ex}.stpfits.trainingSet = TRAINSET;
       dat{i_ex}.stpfits.normalizeData = NORMALIZEDATA;
       dat{i_ex}.stpfits.fitWithAvg = FITAVERAGEDATA;
       dat{i_ex}.stpfits.fitRecovPulse = FITRECOVERYPULSE;
       dat{i_ex}.stpfits.modelParams{i_ch} = [d, f, dTau, fTau];
       dat{i_ex}.stpfits.R2.training{i_ch} = R2_train;
       dat{i_ex}.stpfits.R2.crossvalid{i_ch} = R2_crossvalid;

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
    
    % store some metadata
    pprpop.info{i_ex} = dat{i_ex}.info;
    
    
    % aggregate data within TF conditons (Separately for each
    % recording channel
    for i_tf = 1:numel(uniqueTFs);
        
        tfidx = find(trainParams(:,3) == uniqueTFs(i_tf)); % condnames that contain the train with a particular TF
        
        for i_ch = 1:2
            
            % check to make sure this neuron was defined
            isvalid = dat{i_ex}.info.HS_is_valid(i_ch);
            if ~isvalid
                pprpop.dat{i_ex}.xbar{i_tf}{i_ch} = [];
                pprpop.tfs{i_ex}{i_ch} = [];
                continue
            end
            
            % iterate over the trains with the same freq. 
            catdat = [];
            for i_cond = 1:numel(tfidx)
                
                % pull out the data
                condIdx = ismember(allStimParams, trainParams(tfidx(i_cond),:), 'rows');
                assert(sum(condIdx)==1, 'ERROR, found zero or more than 1 instance of the trial type')
                tmpdat = dat{i_ex}.expt.(condnames{condIdx}).stats.EPSCamp{i_ch}; % [Npulses, 1, Nsweeps]
                if ~isempty(tmpdat)
                    % normalize by the p1Amp_norm
                    realTrlNums = dat{i_ex}.expt.(condnames{condIdx}).realTrialNum{i_ch};
                    p1Amp_norm = dat{i_ex}.qc.p1amp_norm{i_ch}(realTrlNums);
                    assert(~any(isnan(p1Amp_norm)), 'ERROR: scale factor is a nan');
                    p1Amp_norm = permute(p1Amp_norm, [1,3,2]);
                    tmpdat = bsxfun(@rdivide, tmpdat, p1Amp_norm); % [Npulses, 1, Nsweeps]
                    
                    % store in a matrix for averaging later
                    catdat = cat(3, catdat, tmpdat);
                end
            end
            
            % check to make sure there are data for these conditions. Even
            % though this recording channel should be defined (See above),
            % it's possible there are no data due to deletion of single
            % sweeps
            if isempty(catdat)
                pprpop.dat{i_ex}.xbar{i_tf}{i_ch} = [];
                pprpop.tfs{i_ex}{i_ch} = [];
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
                    ppr(end) = [];
                end
                
                pprpop.dat{i_ex}.xbar{i_tf}{i_ch} = ppr;
                pprpop.tfs{i_ex}{i_ch} = uniqueTFs;
                pprpop.MaxNPulses = max([pprpop.MaxNPulses, numel(ppr)]);
                
                
                
                % make a smooth manifold for this neuron.
                % Assume TF = 10 : 50;
                if isfield(dat{i_ex}, 'stpfits')
                    params = dat{i_ex}.stpfits.modelParams{i_ch};
                    isi_ms = fliplr([1000/50 : 1 : 1000/10]);
                    NumPulses = 10;
                    
                    
                    smoothManifold = nan(NumPulses, numel(isi_ms));
                    for i_isi = 1:numel(isi_ms)
                        A0 = 1;
                        tmp_pOntimes_ms = 0 : isi_ms(i_isi) : (isi_ms(i_isi)*NumPulses)-1;
                        tmp_pOntimes_sec = tmp_pOntimes_ms ./ 1000;
                        smoothManifold(:,i_isi) = predictPSCfromTau(tmp_pOntimes_sec, params(1:2), params(4:5), params(3), params(6), A0);
                    end
                    pprpop.smoothManifold{i_ex}{i_ch} = smoothManifold;
                    pprpop.smoothManifold_isi{i_ex}{i_ch} = isi_ms;
                end
                
            end
        end
        
    end
    
end


%% PAIRED PULSE PLASTICITY MANIFOLDS (PLOTS)
clc; close all

PLOT_SMOOTH_MANIFOLD = true;
PLOT_RAW_DATA = true;
PLOT_INTERP_RAW_DATA = false;
PLOT_INDIVIDUAL_DATASETS = false;

% define a set of attributes for each lineseries (or manifold) in the plot
% {CellType,  BrainArea,  OpsinType}
% Brain Area can be: 'AL', 'PM', 'AM', 'LM', 'any', 'med', 'lat'. CASE SENSITIVE
plotgroups = {
    'PY',    'LM', 'chronos';...
    'PY',    'AM', 'chronos';...
    };

groupdata_raw = repmat({[]}, 1, size(plotgroups, 1)); % should only have N cells, where N = size(plotgroups, 1). Each cell has a matrix with a cononicalGrid:
groupdata_smooth =  repmat({[]}, 1, size(plotgroups, 1));
groupexpinds =  repmat({[]}, 1, size(plotgroups, 1));
groupchinds = repmat({[]}, 1, size(plotgroups, 1));
canonicalGrid = nan(pprpop.MaxNPulses, numel(pprpop.TFsAllExpts)); % need to get the maxNpulses into the pprpop struct
allTFs = pprpop.TFsAllExpts;

% iterate over the experiments. For each recording channel, determine what
% the attributes are, and place the data in the correct ploting group.
for i_ex = 1:numel(dat)
    
    for i_ch = 1:2
        
        % check to make sure this neuron was defined
        isvalid = dat{i_ex}.info.HS_is_valid(i_ch);
        if ~isvalid
            continue
        end
        
        % check the attributes
        ch_attribs = {dat{i_ex}.info.cellType{i_ch}, upper(dat{i_ex}.info.brainArea), dat{i_ex}.info.opsin}; % force the brain area to be uppercase
        l_nan = cellfun(@(x) all(isnan(x)), ch_attribs);
        ch_attribs(l_nan) = cellfun(@num2str, ch_attribs(l_nan), 'uniformoutput', false);
        
        % is this expt and channel cooresponds to one of the plot_groups in
        % terms of cellType and opsin
        l_cellType_match = cellfun(@(x) ~isempty(regexpi(ch_attribs{1}, x)), plotgroups(:,1)) | strcmpi(plotgroups(:,1), 'any');
        l_opsinMatch = cellfun(@(x) ~isempty(regexpi(ch_attribs{3}, x)), plotgroups(:,3)) | strcmpi(plotgroups(:,3), 'any');
        
        % determine if this experiment has a brain area that corresponds to
        % one of the plot_groups. Start by adding a 'medial', 'lateral'
        % assignment to the brain area
        expt_area = ch_attribs{2};
        if ~isempty(regexp(expt_area, 'AM', 'once')) || ~isempty(regexp(expt_area, 'PM', 'once'))
            expt_area = [expt_area, ' med'];
        elseif ~isempty(regexp(expt_area, 'AL', 'once')) || ~isempty(regexp(expt_area, 'LM', 'once'))
            expt_area = [expt_area, ' lat'];
        end
        l_brainArea_match = cellfun(@(x) ~isempty(regexp(expt_area, x, 'once')), plotgroups(:,2)) | strcmpi(plotgroups(:,2), 'any');
        
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
        groupdata_raw{group_idx} = cat(3, groupdata_raw{group_idx}, tmpgrid);
        groupexpinds{group_idx} = cat(1, groupexpinds{group_idx}, i_ex);
        groupchinds{group_idx} = cat(1, groupchinds{group_idx}, i_ch);
        
        if PLOT_SMOOTH_MANIFOLD
            groupdata_smooth{group_idx} = cat(3, groupdata_smooth{group_idx}, pprpop.smoothManifold{i_ex}{i_ch});
        end
    end
end


% Plot average manfold. any global (above) that has multiple elements will
% be ploted against eachother, but the other globals will be held constant.
hf = figure;
hold on,
plotcolors = {'r', 'b', 'g'};
for i_group = 1:numel(groupdata_raw)
    
    if PLOT_RAW_DATA
        grid_average = nanmean(groupdata_raw{i_group},3);
        grid_N = sum(~isnan(groupdata_raw{i_group}),3)
        
        Y = 1:size(groupdata_raw{i_group},1);
        X = allTFs';
        
        l_nan_tfs = all(isnan(grid_average), 1);
        
        hs = surf(X(:,~l_nan_tfs), Y, flipud(grid_average(:,~l_nan_tfs)));
        hs.EdgeColor = plotcolors{i_group};
        hs.EdgeAlpha = 1;
        hs.LineWidth = 1.5;
        hs.FaceColor = plotcolors{i_group};
        hs.FaceAlpha = 0;
    end
    
    % plot the average smoothManifold
    if PLOT_SMOOTH_MANIFOLD
        grid_average = mean(groupdata_smooth{i_group}, 3);
        isi_ms = fliplr([1000/50 : 1 : 1000/10]);
        X = 1000./isi_ms;
        Y = 1:size(smoothManifold,1);
        hmod = surf(X,Y, flipud(grid_average));
        hmod.EdgeAlpha = 0;
        hmod.FaceColor = plotcolors{i_group};
        hmod.FaceAlpha = 0.5;
    end
    
end
set(gca, 'zscale', 'log', 'view', [-43    16])
set(gca, 'YTick', 1:10, 'YTickLabel', {'10','9','8','7','6','5','4','3','2','1'})
zmax = get(gca, 'zlim');
xlabel('Temporal Frequency')
ylabel('Pulse Number')
zlabel('norm amp')


% plot manifolds for each dataset individually
if PLOT_INDIVIDUAL_DATASETS
    for i_group = 1:numel(groupdata_raw)
        
        for i_examp = 1:size(groupdata_raw{i_group},3)
            f = figure;
            
            grid_average = groupdata_raw{i_group}(:,:,i_examp);
            
            Y = 1:size(groupdata_raw{i_group},1);
            X = allTFs';
            
            l_nan_tfs = all(isnan(grid_average), 1);
            
            hs = surf(X(:,~l_nan_tfs), Y, flipud(grid_average(:,~l_nan_tfs)));
            hs.EdgeColor = plotcolors{i_group};
            hs.EdgeAlpha = 1;
            hs.FaceColor = plotcolors{i_group};
            hs.FaceAlpha = 0;
            hs.LineWidth = 1.5;
            
            set(gca, 'zscale', 'log', 'view', [-43    16])
            set(gca, 'ytick', Y, 'yticklabel', cellfun(@(x) num2str(x), num2cell(rot90(Y,2)), 'uniformoutput', false))
            set(gca, 'xtick', X, 'xticklabel', cellfun(@(x) num2str(x), num2cell(allTFs), 'uniformoutput', false))
            xlabel('Temporal Frequency')
            ylabel('Pulse Number')
            zlabel('norm amp')
            
            i_ex = groupexpinds{i_group}(i_examp);
            i_ch = groupchinds{i_group}(i_examp);
            f.Name = sprintf('Mouse %s, site %s, HS%d', dat{i_ex}.info.mouseName, dat{i_ex}.info.siteNum, i_ch);
            
            % display the smooth manifold
            smoothManifold = pprpop.smoothManifold{i_ex}{i_ch};
            isi_ms = pprpop.smoothManifold_isi{i_ex}{i_ch};
            X = 1000./isi_ms;
            Y = 1:size(smoothManifold,1);
            hold on,
            hmod = surf(X,Y, flipud(smoothManifold));
            hmod.EdgeAlpha = 0;
            hmod.FaceAlpha = 0.5;
            
            
            
        end
        
    end
end


%% POST-TETANIC POTENTIATION (DATA COLLECTION)

% loop through the experiments. Pull out the trains data. Ignore the
% recovery train (if present) and aggregate across recovery conditions.
recovpop = [];
recovpop.TFsAllExpts = [];
recovpop.recoveryTimesAllExpts = [];
for i_ex = 1:numel(dat)
    
    % find the normal trains. Assume the field name is NOT 'ritv'
    condnames = fieldnames(dat{i_ex}.expt);
    l_trains = ~strncmp(condnames, 'RITv', 4);
    if sum(l_trains)==0; continue; end % no trains data
    trainParams = cellfun(@(x) dat{i_ex}.expt.(condnames{x}).tdict, num2cell(find(l_trains)), 'uniformoutput', false);
    trainParams = cat(1, trainParams{:});
    recovpop.trainParams{i_ex} = trainParams;
    
    % make sure the pulse amplitude and width were identical across ttypes
    assert(numel(unique(trainParams(:,1)))==1, 'ERROR: more than one pulse amplitude')
    assert(numel(unique(trainParams(:,2)))==1, 'ERROR: more than one pulse width');
    
    % store the stim params for all stim types, which will be useful for
    % indexing later.
    allStimParams = cellfun(@(x) dat{i_ex}.expt.(condnames{x}).tdict, num2cell((1:numel(l_trains))'), 'uniformoutput', false);
    allStimParams = cat(1, allStimParams{:});
    
    
    % aggregate data across TF conditons and recovery times (Separately for each
    % recording channel)
    
    for i_ch = 1:2
        
        % check to make sure this neuron was defined
        isvalid = dat{i_ex}.info.HS_is_valid(i_ch);
        
        for i_cond = 1:sum(l_trains)
            
            condIdx = ismember(allStimParams, trainParams(i_cond,:), 'rows');
            
            % check to make sure this is a recovery train
            isrecovery = allStimParams(condIdx,4) > 0;
            if ~isvalid || ~isrecovery
                recovpop.dat{i_ex}.recovAmp{i_ch}(i_cond,1) = NaN;
                continue
            end
            
            % pull out the data
            tmpdat = dat{i_ex}.expt.(condnames{condIdx}).stats.EPSCamp{i_ch}; % [Npulses, 1, Nsweeps]
            

            % check to make sure there are data for these conditions. Even
            % though this recording channel should be defined (See above),
            % it's possible there are no data due to deletion of single
            % sweeps
            if isempty(tmpdat)
                recovpop.dat{i_ex}.recovAmp{i_ch}(i_cond,1) = NaN;
                continue
            else
                
                % normalize by the p1Amp_norm
                realTrlNums = dat{i_ex}.expt.(condnames{condIdx}).realTrialNum{i_ch};
                p1Amp_norm = dat{i_ex}.qc.p1amp_norm{i_ch}(realTrlNums);
                assert(~any(isnan(p1Amp_norm)), 'ERROR: scale factor is a nan');
                p1Amp_norm = permute(p1Amp_norm, [1,3,2]);
                tmpdat = bsxfun(@rdivide, tmpdat, p1Amp_norm); % [Npulses, 1, Nsweeps]
                
                
                % now average across trails, and re-normalize to the first
                % pulse. Store in the population data structure.
                avg = mean(tmpdat,3);
                ppr = avg./avg(1);
                
                recovpop.dat{i_ex}.recovAmp{i_ch}(i_cond,1) = ppr(end);
                
            end
            
        end
        
    end
    
    %
    % update these things, but only if there were recovery data aggregated
    % for this experiment
    %
    dataPresent = cellfun(@(x) any(~isnan(x)), recovpop.dat{i_ex}.recovAmp);
    if any(dataPresent)
        
        % identify the unique TF conditions for this experiment, and update the
        % running log of TFs used across all experiments. This will be used by
        % the plotting routines to figure out the "canonical grid"
        uniqueTFs = unique(trainParams(:,3));
        tmp = cat(1, recovpop.TFsAllExpts, uniqueTFs);
        recovpop.TFsAllExpts = unique(tmp);
        
        % identify the unique recovery times for this experiment, and update the
        % running log of recov times used across all experiments. This will be used by
        % the plotting routines to figure out the "canonical grid"
        uniqueRecoveryTimes = unique(trainParams(:,4));
        tmp = cat(1, recovpop.recoveryTimesAllExpts, uniqueRecoveryTimes);
        recovpop.recoveryTimesAllExpts = unique(tmp);
        
    end
    
end



%% POST-TETANIC POTENTIATION (PLOTS)

% plot recovery pulse amplitude (normalized) vs. train frequency. Do this
% separately for cell types, brain areas, opsins, etc...
clc


% define a set of attributes for each lineseries (or manifold) in the plot
% {CellType,  BrainArea,  OpsinType}
plotgroups = {
              'PY',    'AL', 'chronos';...
              'PY',    'PM', 'chronos';...
              };

groupdata_raw = repmat({[]}, 1, size(plotgroups, 1)); % should only have N cells, where N = size(plotgroups, 1). Each cell has a matrix with a cononicalGrid:
canonicalGrid = nan(numel(recovpop.TFsAllExpts), numel(recovpop.recoveryTimesAllExpts)); % need to get the maxNpulses into the recovpop struct
allTFs = recovpop.TFsAllExpts;
allRecoveryTimes = recovpop.recoveryTimesAllExpts';

% iterate over the experiments. For each recording channel, determine what
% the attributes are, and place the data in the correct ploting group.
for i_ex = 1:numel(dat)
    
    for i_ch = 1:2
        % check to make sure this neuron was defined
        isvalid = dat{i_ex}.info.HS_is_valid(i_ch);
        if ~isvalid
            continue
        end
        
        % check the attributes
        ch_attribs = {dat{i_ex}.info.cellType{i_ch}, dat{i_ex}.info.brainArea, dat{i_ex}.info.opsin};
        l_nan = cellfun(@(x) all(isnan(x)), ch_attribs);
        ch_attribs(l_nan) = cellfun(@num2str, ch_attribs(l_nan), 'uniformoutput', false);
        
        % is this expt and channel cooresponds to one of the plot_groups in
        % terms of cellType and opsin
        l_cellType_match = cellfun(@(x) ~isempty(regexpi(ch_attribs{1}, x)), plotgroups(:,1)) | strcmpi(plotgroups(:,1), 'any');
        l_opsinMatch = cellfun(@(x) ~isempty(regexpi(ch_attribs{3}, x)), plotgroups(:,3)) | strcmpi(plotgroups(:,3), 'any');
        
        
        % determine if this experiment has a brain area that corresponds to
        % one of the plot_groups. Start by adding a 'medial', 'lateral'
        % assignment to the brain area
        expt_area = ch_attribs{2};
        if ~isempty(regexp(expt_area, 'AM', 'once')) || ~isempty(regexp(expt_area, 'PM', 'once'))
            expt_area = [expt_area, ' med'];
        elseif ~isempty(regexp(expt_area, 'AL', 'once')) || ~isempty(regexp(expt_area, 'LM', 'once'))
            expt_area = [expt_area, ' lat'];
        end
        l_brainArea_match = cellfun(@(x) ~isempty(regexp(expt_area, x, 'once')), plotgroups(:,2)) | strcmpi(plotgroups(:,2), 'any');
        
        group_idx = sum([l_cellType_match, l_brainArea_match, l_opsinMatch], 2) == 3;
        assert(sum(group_idx)<=1, 'ERROR: found too many group indicies')
        if sum(group_idx) == 0; continue; end

        
        % add data to the appropriate group data array
        tmpgrid = canonicalGrid;
        for i_cond = 1:size(recovpop.trainParams{i_ex}, 1);
            grid_row_idx = allTFs == recovpop.trainParams{i_ex}(i_cond, 3);
            grid_col_idx = allRecoveryTimes == recovpop.trainParams{i_ex}(i_cond, 4);
            tmpgrid(grid_row_idx, grid_col_idx) = recovpop.dat{i_ex}.recovAmp{i_ch}(i_cond);
        end
        
        if all(isnan(tmpgrid(:))); continue; end
        
        % aggregate the data
        groupdata_raw{group_idx} = cat(3, groupdata_raw{group_idx}, tmpgrid);
        
    end
end




% Plot average lineseries.
hf = figure;
hold on,
plotcolors = [1,0,0;...
              0,1,0;...
              0,0,1];
          
for i_group = 1:numel(groupdata_raw)
    
    grid_average = nanmean(groupdata_raw{i_group},3);
    if isempty(grid_average); continue; end
    
    N = sum(~isnan(groupdata_raw{i_group}),3)
    grid_sem = nanstd(groupdata_raw{i_group},[],3) ./ sqrt(N);
    
    groupcolors = repmat(plotcolors(i_group,:), size(grid_average,1), 1);
    ramp = linspace(0,0.8, size(grid_average,1))';
    l_off = sum(groupcolors, 1)==0;
    groupcolors(:, l_off) = [ramp, ramp];
    
    % loop over TF conds
    for i_tf = 1:size(grid_average, 1)
        l_nan = isnan(grid_average(i_tf,:));
        tmp_amps = grid_average(i_tf, ~l_nan);
        tmp_sem = grid_sem(i_tf, ~l_nan);
        tmp_recovTimes = allRecoveryTimes(~l_nan);
        my_errorbar(tmp_recovTimes, tmp_amps, tmp_sem, '-', 'color', groupcolors(i_tf,:), 'linewidth', 3);
    end
    
end




