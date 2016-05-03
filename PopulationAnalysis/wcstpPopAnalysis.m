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

% pool = gcp('nocreate');
% if isempty(pool)
%     pool = parpool(8);
% end

for i_ex = 1:Nexpts
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
    
    % make a mini dataset that's composed of only the RIT data, then
    % predict responses to recovery trains.
    ritdata = [];
    condnames = fieldnames(dat{i_ex}.expt);
    for i_cond = 1:numel(condnames)
        if strncmp(condnames{i_cond}, 'RITv', 4)
            ritdata.expt.(condnames{i_cond}) = dat{i_ex}.expt.(condnames{i_cond});
        end
    end
    
   for i_ch = 1:2
       
       fldname = sprintf('HS%dvalid', i_ch);
       isvalid = str2double(wb_expt(i_ex, hidx.(fldname)));
       if isempty(ritdata) || ~isvalid
           continue
       end
       
       % fit the RIT data
       [d, f, dTau, fTau] = fitTau2STP(ritdata, 'EPSCamp', i_ch, 'multistart');
       
       % predict all the data
       A0 = mean(dat{i_ex}.qc.p1amp{i_ch});
       ttypes = fieldnames(dat{i_ex}.expt);
       [pred, raw_xbar, raw_sem] = deal({});
       for i_type = 1:numel(ttypes)
           pOnTimes = dat{i_ex}.expt.(ttypes{i_type}).pOnTimes;
           raw = dat{i_ex}.expt.(ttypes{i_type}).stats.EPSCamp{i_ch};
           raw_xbar{i_ch}{i_type} = mean(raw,3);
           raw_sem{i_ch}{i_type} = stderr(raw,3);
           
           pred{i_ch}{i_type} = predictPSCfromTau(pOnTimes, d, dTau, f, fTau, A0);
       end
       
   end
end








