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
    
    for i_ch = 1:2
        
        fldname = sprintf('HS%dvalid', i_ch);
        isvalid = str2double(wb_expt(i_ex, hidx.(fldname)));
        if ~isvalid
            continue
        end
        
        % series resistance
        if ~all(isnan(dat{i_ex}.qc.Rs{i_ch}))
            subplot(3,2,i_ch)
            tmp = squeeze(dat{i_ex}.qc.Rs{i_ch});
            plot(tmp)
            ylim([min([0,min(tmp(:))]), max(tmp(:))*1.5])
            ylabel('R_{s} (MOhm)')
            xlabel('trial number')
            title(sprintf('Channel %d', i_ch))
        end
        
        % verr
        if ~all(isnan(dat{i_ex}.qc.verr{i_ch}))            
            subplot(3,2,2+i_ch)
            tmp = squeeze(dat{i_ex}.qc.verr{i_ch});
            plot(tmp)
            ylim([min([0,min(tmp(:))]), max(tmp(:))*1.5])
            ylabel('SS Verr (mV)')
            xlabel('trial number')
            
        end
        
        % p1amps
        if ~all(isnan(dat{i_ex}.qc.p1amp{i_ch}))
            subplot(3,2,4+i_ch)
            tmp = squeeze(dat{i_ex}.qc.p1amp{i_ch});
            plot(tmp)
            ylim([min([0,min(tmp(:))]), max(tmp(:))*1.5])
            ylabel('P1 amplitude')
            xlabel('trial number')
        end
            
        
    end
    drawnow
    
end



%% STP SUMMARY FOR EACH RECORDING

close all

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
            
            % stem plot
            subplot(Nconds,2, 2.*(i_cond-1) + i_ch);
            tmp = dat{i_ex}.expt.(conds{i_cond}).stats.EPSCamp{i_ch};
            xbar = mean(tmp,3);
            sem = stderr(tmp,3);
            tt = dat{i_ex}.expt.(conds{i_cond}).pOnTimes;
            my_errorbar(tt, xbar, sem, 'ok', 'markersize', 3, 'linewidth', 1);
            xlim([0, dat{i_ex}.info.sweepLength.vclamp]);
            set(gca, 'yticklabel', []);
            if i_cond < Nconds
                set(gca, 'xticklabel', []);
            end

            
        end
    end
    drawnow
    
end




%% ESTIMATE TIME CONSTANTS OF SHORT TERM PLASTICITY


for i_ex = 1:numel(dat)
    
   hf = figure;
   hf.Position = [300 27 1325 957];
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
       
       % make a mini dataset that's composed of only the RIT data. I need the
       % pOnTimes, raw EPSC amps.
       [pOnTimes, rawAmps] = deal({});
       condnames = fieldnames(dat{i_ex}.expt);
       counter = 1;
       for i_cond = 1:numel(condnames)
           if strncmp(condnames{i_cond}, 'RITv', 4)
               pOnTimes{counter} = dat{i_ex}.expt.(condnames{i_cond}).pOnTimes;
               rawAmps{counter} = dat{i_ex}.expt.(condnames{i_cond}).stats.EPSCamp{i_ch};
               rawAmps{counter} = mean(rawAmps{counter},3); % mean across sweeps
               counter = counter+1;
           end
       end
       
       % if there were no RIT data, then move along,
       if isempty(rawAmps)
           chempty(i_ch) = true;
           if i_ch == 1
               continue
           elseif all(chempty)
               close(hf)
               continue
           end
       end
           
       
       % fit the RIT data
       [d, f, dTau, fTau] = fitTau2STP(rawAmps, pOnTimes, 'multistart');
       
       % predict all the data
       ttypes = fieldnames(dat{i_ex}.expt);
       [pred, raw_xbar, raw_sem] = deal({});
       for i_type = 1:numel(ttypes)
           pOnTimes = dat{i_ex}.expt.(ttypes{i_type}).pOnTimes;
           raw = dat{i_ex}.expt.(ttypes{i_type}).stats.EPSCamp{i_ch};
           raw_xbar{i_ch}{i_type} = mean(raw,3);
           raw_sem{i_ch}{i_type} = stderr(raw,3);
           A0 = raw_xbar{i_ch}{i_type}(1);
           pred{i_ch}{i_type} = predictPSCfromTau(pOnTimes, d, dTau, f, fTau, A0);
       end
       
       
       % plot the recovery trains (cross validation)
       figure(hf)
       trainNames = fieldnames(dat{i_ex}.expt);
       l_recov = ~(strncmp(trainNames, 'RIT', 3));
       recovIdx = find(l_recov);
       xlims = [inf -inf];
       ylims = [inf -inf];
       hs = [];
       for i_type = 1:numel(recovIdx)
           
           typeIdx = recovIdx(i_type);
           pltIdx = sub2ind([4, numel(recovIdx)], i_ch+1, i_type);
           hs(i_type) = subplot(numel(recovIdx), 4, pltIdx); hold on,
           
           xx = dat{i_ex}.expt.(trainNames{typeIdx}).pOnTimes;
           my_errorbar(xx, raw_xbar{i_ch}{typeIdx}, raw_sem{i_ch}{typeIdx}, 'k');
           plot(xx, pred{i_ch}{typeIdx}, 'r', 'linewidth', 2)
           xlims(1) = min([min(xx), xlims(1)]);
           xlims(2) = max([max(xx), xlims(2)]);
           yvals = get(gca, 'ylim');
           ylims(1) = min([yvals(1), ylims(1)]);
           ylims(2) = max([yvals(2), ylims(2)]);
       end
       set(hs, 'XLim', xlims, 'YLim', ylims)
       
       
       % make a scatter plot of all predicted and actual PSC amps
       hs = [];
       all_raw = [];
       all_pred = [];
       crossval_raw = [];
       crossval_pred = [];
       for i_type = 1:numel(trainNames)
           if strncmpi(trainNames{i_type}, 'rit', 3); pltclr = 'k';else pltclr = 'r';end
           if i_ch == 1; pltcol=1; else pltcol=4; end
           
           tmp_raw = dat{i_ex}.expt.(trainNames{i_type}).stats.EPSCamp{i_ch};
           tmp_pred = pred{i_ch}{i_type};
           tmp_pred = repmat(tmp_pred, size(tmp_raw,3), 1);
           tmp_raw = tmp_raw(:);
           assert(all(size(tmp_pred) == size(tmp_raw)))
           
           pltIdx = sub2ind([4, 3], pltcol, 1);
           hs = subplot(3, 4, pltIdx); hold on,
           plot(tmp_raw./tmp_raw(1), tmp_pred./tmp_pred(1), '.', 'color', pltclr)
           axis tight
           
           % concatenate all the data
           all_raw = cat(1, all_raw, tmp_raw(:));
           all_pred = cat(1, all_pred, tmp_pred(:));
           
           % concatenate the cross-validation data
           if ~strncmpi(trainNames{i_type}, 'rit', 3)
               crossval_raw = cat(1, crossval_raw, tmp_raw(:));
               crossval_pred = cat(1, crossval_pred, tmp_pred(:));
           end
       end
       maxval = max([hs.XLim, hs.YLim]);
       minval = min([hs.XLim, hs.YLim]);
       plot([minval, maxval], [minval, maxval], 'k--')
       hs.XScale = 'log';
       hs.YScale = 'log';
       xlabel('raw EPSC amp (norm)')
       ylabel('pred amp (norm)')
       
       pltIdx = sub2ind([4, 3], pltcol, 2);
       subplot(3,4,pltIdx), hold on,
       resid = all_pred - all_raw;
       histogram(resid)
       plot(mean(resid), 10, 'rv', 'markerfacecolor', 'r')
       R2 = 1 - (sum(resid.^2) ./ sum(all_raw.^2));
       xlabel('pred-real')
       
       % calculate the suffle corrected R2.
       N = numel(all_pred);
       iters = 5000;
       shuffle_inds = unidrnd(N, [N, iters]);
       shuffle_pred = all_pred(shuffle_inds);
       resid = bsxfun(@minus, shuffle_pred, all_raw);
       SS_resid = sum(resid.^2, 2);
       SS_raw = sum(all_raw.^2);
       R2_shuffle = 1 - bsxfun(@rdivide, SS_resid, SS_raw);
       title(sprintf('R2: %.3f, R2_shuff: %.3f', R2, mean(R2_shuffle)))
       
       
       % plot cross validation stuff
       pltIdx = sub2ind([4, 3], pltcol, 3);
       subplot(3,4,pltIdx), hold on,
       resid = crossval_pred - crossval_raw;
       histogram(resid);
       plot(mean(resid), 5, 'rv', 'markerfacecolor', 'r')
       R2 = 1 - (sum(resid.^2) ./ sum(all_raw.^2));
       xlabel('cross-valid (pred-real)')
       
       N = numel(crossval_pred);
       iters = 5000;
       shuffle_inds = unidrnd(N, [N, iters]);
       shuffle_pred = crossval_pred(shuffle_inds);
       resid = bsxfun(@minus, shuffle_pred, crossval_raw);
       SS_resid = sum(resid.^2, 2);
       SS_raw = sum(crossval_raw.^2);
       R2_shuffle = 1 - bsxfun(@rdivide, SS_resid, SS_raw);
       title(sprintf('R2: %.3f, R2_shuff: %.3f', R2, mean(R2_shuffle)))
       
   end
   drawnow
   
   
end








