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
vclamp_fpath = cellfun(@(x,y) strcat(GL_DATPATH, x, filesep, 'Physiology', filesep, y, '.abf'), wb_expt(:,hidx.MouseName), wb_expt(:, hidx.ABFOptostim), 'uniformoutput', false);
wb_expt(:,hidx.ABFDCsteps) = iclamp_fpath;
wb_expt(:,hidx.ABFOptostim) = vclamp_fpath;                 
                   
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
    pool = parpool(min([32, Nexpts]));
end

parfor i_ex = 1:Nexpts
    dat{i_ex} = wcstp_compile_data(attributes{i_ex}, hidx, params);
end



%% QULAITY CONTROL PLOTS

close all

for i_ex = 1:numel(dat)
    
    if all(cellfun(@isempty, dat{i_ex}.qc.Rs))
        continue
    end
    
    f = figure;
    f.Name = sprintf('Mouse %s, site %s', dat{i_ex}.info.mouseName, dat{i_ex}.info.siteNum);
    
    for i_ch = 1:2
        
        if ~isempty(dat{i_ex}.qc.Rs{i_ch})
            
            % series resistance
            subplot(3,2,i_ch)
            tmp = dat{i_ex}.qc.Rs{i_ch};
            plot(tmp)
            ylim([min([0,min(tmp(:))]), max(tmp(:))*1.5])
            ylabel('R_{s} (MOhm)')
            xlabel('pulse number')
            title(sprintf('Channel %d', i_ch))
            
            % vhold
            subplot(3,2,2+i_ch)
            tmp = dat{i_ex}.qc.vhold{i_ch};
            plot(tmp)
            ylim([min([0,min(tmp(:))]), max(tmp(:))*1.5])
            ylabel('SS Verr (mV)')
            xlabel('pulse number')
            
            % p1amps
            subplot(3,2,4+i_ch)
            tmp = dat{i_ex}.qc.p1amp{i_ch};
            plot(tmp)
            ylim([min([0,min(tmp(:))]), max(tmp(:))*1.5])
            ylabel('P1 amplitude')
            xlabel('pulse number')
            
        end
    end
    drawnow
    
end



%% STP SUMMARY FOR EACH RECORDING

close all

for i_ex = 1:numel(dat)
    
    if all(cellfun(@isempty, dat{i_ex}.qc.Rs))
        continue
    end
    
    f = figure;
    f.Name = sprintf('Mouse %s, site %s', dat{i_ex}.info.mouseName, dat{i_ex}.info.siteNum);
    
    for i_ch = 1:2
        
        conds = fieldnames(dat{i_ex}.expt);
        Nconds = numel(conds);
        for i_cond = 1:Nconds
            if isempty(dat{i_ex}.expt.(conds{i_cond}).stats.EPSCamp{i_ch})
                continue
            end
            
            % stem plot
            subplot(Nconds,2, 2.*(i_cond-1) + i_ch)
            tmp = dat{i_ex}.expt.(conds{i_cond}).stats.EPSCamp{i_ch};
            xbar = mean(tmp,3);
            sem = stderr(tmp,3);
            tt = dat{i_ex}.expt.(conds{i_cond}).pOnTimes;
            my_errorbar(tt, xbar, sem, 'ok', 'markersize', 3, 'linewidth', 1)
            xlim([0, dat{i_ex}.info.sweepLength.vclamp]);
            set(gca, 'yticklabel', [])
            if i_cond < Nconds
                set(gca, 'xticklabel', [])
            end

            
        end
    end
    drawnow
    
end


