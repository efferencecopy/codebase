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

% pool = gcp('nocreate');
% if isempty(pool)
%     pool = parpool(16);
% end

for i_ex = 1:Nexpts
    dat{i_ex} = wcstp_compile_data(attributes{i_ex}, hidx, params);
end

fprintf('All done importing data\n')

%% QULAITY CONTROL PLOTS

close all

for i_ex = 1:numel(dat)
    
    if all(cellfun(@isempty, dat{i_ex}.qc.Rs))
        continue
    end
    
    f = figure;
    f.Name = sprintf('Mouse %s, site %s', dat{i_ex}.info.mouseName, dat{i_ex}.info.siteNum);
    
    for i_ch = 1:2
        
        % series resistance
        if ~isempty(dat{i_ex}.qc.Rs{i_ch})
            subplot(3,2,i_ch)
            tmp = dat{i_ex}.qc.Rs{i_ch};
            plot(tmp)
            ylim([min([0,min(tmp(:))]), max(tmp(:))*1.5])
            ylabel('R_{s} (MOhm)')
            xlabel('trial number')
            title(sprintf('Channel %d', i_ch))
        end
        
        % vhold
        if ~isempty(dat{i_ex}.qc.vhold{i_ch})            
            subplot(3,2,2+i_ch)
            tmp = dat{i_ex}.qc.vhold{i_ch};
            plot(tmp)
            ylim([min([0,min(tmp(:))]), max(tmp(:))*1.5])
            ylabel('SS Verr (mV)')
            xlabel('trial number')
            
        end
        
        % p1amps
        if ~isempty(dat{i_ex}.qc.p1amp{i_ch})
            subplot(3,2,4+i_ch)
            tmp = dat{i_ex}.qc.p1amp{i_ch};
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
    
    if all(cellfun(@isempty, dat{i_ex}.qc.Rs))
        continue
    end
    
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

clc
i_ex = 8;
channel = 2;
psctype = 'EPSCamp';

[d, f, dTau, fTau] = fitTau2STP(dat{i_ex}, psctype, channel, 'multistart');

ttypes = fieldnames(dat{i_ex}.expt);
for i_type = 1:numel(ttypes)
    pOnTimes = dat{i_ex}.expt.(ttypes{i_type}).pOnTimes;
    raw = dat{i_ex}.expt.(ttypes{i_type}).stats.(psctype){channel};
    
    A0 = mean(raw(1,1,:), 3);
    pred = predictPSCfromTau(pOnTimes, d, dTau, f, fTau, A0);

    figure
    plot(pred, '.b-'), hold on, stem(mean(raw(:,1,:),3), 'k');
end



%% TEST CODE FOR STP TIME CONSTANTS (1)

%
% compare fitted and 'real' parameters (without noise added)
%

clc, close all

% assume that pOnTimes exists. A cell array of pulse on times. Define some
% plasticity params, and then generate fake data. Add noise if desired

% setup some pulse times for the simulation
pOnTimes = dat{9}.expt.RITv8.pOnTimes;


% setup a bunch of different plasticity values
% version with facilitation and depression
d1_sim = [0.020:0.50:1];
d2_sim = [0.020:0.50:1];
tau_d1_sim = logspace(log10(0.001), log10(15), 5);
tau_d2_sim = logspace(log10(0.001), log10(15), 5);
f1_sim = logspace(log10(0.020), log10(15), 5);
tau_f1_sim = logspace(log10(0.001), log10(15), 5);
A0_sim = logspace(log10(10), log10(1000), 5);

[a,b,c,d,e,f,g] = ndgrid(d1_sim, d2_sim, tau_d1_sim, tau_d2_sim, f1_sim, tau_f1_sim, A0_sim);
W = [a(:),b(:),c(:),d(:),e(:),f(:),g(:)];
params = mat2cell(W, ones(size(W,1),1));
clear a b c d e f g W


for i_iter = 1:numel(params)
    fprintf('define: iter = %d\n', i_iter);
    psc_test = predictPSCfromTau(pOnTimes, params{i_iter}(1:2), params{i_iter}(3:4), params{i_iter}(5), params{i_iter}(6), params{i_iter}(7));
    if any(psc_test<0); error('found one below zero'); end
    testdat{i_iter}.expt.RITv1.stats.EPSCamp{1} = psc_test;
    testdat{i_iter}.expt.RITv1.pOnTimes = pOnTimes;
end


params_out = cell(numel(params),1);

parfor i_iter = 1:numel(params)
    fprintf('fit: iter = %d\n', i_iter);
    channel = 1;
    psctype = 'EPSCamp';
    [d_test, f_test, dTau_test, fTau_test] = fitTau2STP(testdat{i_iter}, psctype, channel, 'multistart');
    params_out{i_iter} = [d_test, dTau_test, f_test, fTau_test]
end


params_out = cat(1, params_out{:});
params = cat(1, params{:});


for i_param = 1:size(params,2)-1
    figure
    histogram(params(:,i_param)-params_out(:,i_param))
end







%% TEST CODE FOR STP TIME CONSTANTS (2)

%
% estimate the shape of the error function
%

d1_real = 0.4;
d2_real = 1;
tau_d1_real = 1;
tau_d2_real = 1;
f1_real = 2;
tau_f1_real = 2;
A0 = 700;

pOnTimes = dat{9}.expt.RITv8.pOnTimes;
pscreal = predictPSCfromTau(pOnTimes, [d1_real, d2_real], [tau_d1_real, tau_d2_real], f1_real, tau_f1_real, A0);

[a,b] = ndgrid(0.001:0.005:1, 0:0.050:3);
errvals = nan(size(a));
for idx = 1:numel(a)
        k_d = [a(idx), d2_real];
        tau_d = [tau_d1_real, tau_d2_real];
        k_f = b(idx);
        tau_f = tau_f1_real;

        pred = predictPSCfromTau(pOnTimes, k_d, tau_d, k_f, tau_f, A0);

        % pool the errors across pulse train types. Ingore the first pulse
        % (becuse the error is artifactually = zero). Do some error
        % checking along the way.
        err = pscreal(2:end) ./  pred(2:end);
        err = sum(abs(log10(err))); % big negative powers of ten are good 
        errvals(idx) = err;
end
       
h = surf(b, a, errvals);
h.EdgeAlpha = 0.1;

%% TEST CODE FOR STP TIME CONSTANTS (3)

%
% compare fitted and 'real' parameters with added noise, as a function of
% number of trails
%


clc, close all

% setup some pulse times for the simulation
pOnTimes = dat{9}.expt.RITv8.pOnTimes;
Npulses = numel(pOnTimes);

% setup a bunch of different plasticity values
% version with facilitation and depression
d1_real = 0.7;
d2_real = 0.95;
tau_d1_real = 0.4;
tau_d2_real = 2;
f1_real = 0.7;
tau_f1_real = 0.4;
A0 = 700;

% setup the noise and trial counts
sigma = A0 .* 0.10;
mu = 0;
Ntrials = 500;
sim_iters = 50;


pscreal = predictPSCfromTau(pOnTimes, [d1_real, d2_real], [tau_d1_real, tau_d2_real], f1_real, tau_f1_real, A0);
testdat = cell(sim_iters,1);

for i_iter = 1:sim_iters
    
    psc_iter = repmat(pscreal, [1,1,Ntrials]) + normrnd(mu, sigma, [Npulses, 1, Ntrials]);
    if any(psc_iter(:)<0); error('less than zero'); end
    testdat{i_iter}.expt.RITv1.stats.EPSCamp{1} = psc_iter;
    testdat{i_iter}.expt.RITv1.pOnTimes = pOnTimes;
end


params_out = cell(sim_iters,1);
parfor i_iter = 1:sim_iters
    fprintf('fit: iter = %d\n', i_iter);
    channel = 1;
    psctype = 'EPSCamp';
    [d_test, f_test, dTau_test, fTau_test] = fitTau2STP(testdat{i_iter}, psctype, channel, 'global');
    params_out{i_iter} = [d_test, dTau_test, f_test, fTau_test]
end



params_out = cat(1,params_out{:});
params_real = [d1_real, d2_real, tau_d1_real, tau_d2_real, f1_real, tau_f1_real];

figure
for i_param = 1:6
    subplot(2,3,i_param)
    histogram(params_out(:,i_param) - params_real(i_param))
end
    
    
    

