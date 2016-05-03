%% MAKE SOME PULSE TRAINS: RECOVERY, RIT, ENVELOPED-RIT

fin

PLOTFIGS = false;
CHANNEL = 1;
PSCTYPE = 'EPSCamp';
FITMETHOD = 'multistart';

SIGMAFACTOR = 10;  % noise is 1/sigmaFactor of the PSC amplitude
Niters = 100;
Ntrials = 25;

% simulation params
d1_sim = 0.80;
d2_sim = 0.98;
tau_d1_sim = 0.500;
tau_d2_sim = 5;
f1_sim = 0.4;
tau_f1_sim = 0.200;
A0_sim = 600;

[rit_fits, rit_fits_noiseless, train_fits, train_fits_noiseless] = deal(nan(Niters, 6));

for i_iter = 1:Niters;
    
    i_iter
    
    %
    % Define the params that influence the entire data file (trains and RITs)
    %
    params.si   = 50e-6;              % the sample INTERVAL (needs to be an iteger)
    params.swpDur = 10e4;             % The total duration of the sweep IN NUMBERS OF SAMPLES!!!!
    params.tStart = 0.100;            % the time of the first pulse
    params.pAmp = 1;                  % A vector of amplitudes for the pulse height [interleaved variable]
    params.pWidth = 350e-6;           % A vector of pulse widths (in seconds)  [interleaved variable]
    
    
    %
    % make the sweep templates for the pulse trains
    %
    params.type = 'trains';            % 'train', 'pulse'
    params.pFreq = [10, 20, 40];              % A vector of frequencies for the pulse train [interleaved variable]
    params.nPulses = 10;
    params.recovery = true;       %  true/false, should there be a recovery pulse?
    params.recoveryTime = [0.250, .500, 2, 3.8];   %  A vector of numbers corresponding to the recovery time in seconds [interleaved variable]
    params = makeSweepTemplates_trains(params); % templates are stored in params.templates_trains
    
    
    %
    % make the sweep template(s) for the Random impulse trains
    %
    params.ritFreq = [8];
    params.ritHiFreqCut = 58;  % ISIs faster than this will be cutout
    params.rit_Nversions = numel(params.templates_trains) .* Ntrials;
    params.ritUseEnvelope = true;
    params.ritEnvelopeFreq = [0.20];
    params = makeSweepTemplates_poiss(params); % templates are stored in params.templates_poiss
    
    
    
    % concatenate templates
    sweepTemplates = cat(2, params.templates_trains, params.templates_poiss);
    
    
    % GENERATE NOISELESS SIMULATED RESPONSES
    d = [d1_sim, d2_sim];
    dTau = [tau_d1_sim, tau_d2_sim];
    f = f1_sim;
    fTau = tau_f1_sim;
    A0 = A0_sim;
    
    % start with recovery trains
    for i_tr = 1:numel(params.templates_trains)
        
        % define pOnTimes
        tt = (0:numel(params.templates_trains{i_tr})-1) .* params.si;
        crossing = params.templates_trains{i_tr} > 0.25;
        crossing = [0; diff(crossing)] == 1;
        params.pOnTimes_trains{i_tr} = tt(crossing);
        
        % simulate the EPSCs
        params.psc_test_trains{i_tr} = predictPSCfromTau(params.pOnTimes_trains{i_tr}, d, dTau, f, fTau, A0);
        
        if PLOTFIGS
            figure
            plot(params.pOnTimes_trains{i_tr}, params.psc_test_trains{i_tr}, '-o')
            xlim([tt(1), tt(end)])
        end
        
    end
    for i_tr = 1:numel(params.templates_poiss)
        
        % define the pOnTimes
        tt = (0:numel(params.templates_poiss{i_tr})-1) .* params.si;
        crossing = params.templates_poiss{i_tr} > 0.25;
        crossing = [0; diff(crossing)] == 1;
        params.pOnTimes_poiss{i_tr} = tt(crossing);
        
        %simulate the EPSCs
        params.psc_test_poiss{i_tr} = predictPSCfromTau(params.pOnTimes_poiss{i_tr}, d, dTau, f, fTau, A0);
        
        if PLOTFIGS
            figure
            plot(params.pOnTimes_poiss{i_tr}, params.psc_test_poiss{i_tr}, '-o')
            xlim([tt(1), tt(end)])
        end
        
    end
    
    %
    % OPTIONAL PLOTS OFN THE TEMPLATES AND INSTANTANEOUS FREQ
    %
    if PLOTFIGS
        % start with recovery trains
        for i_tr = 1:numel(params.templates_trains)
            
            tt = (0:numel(params.templates_trains{i_tr})-1) .* params.si;
            
            f = figure;
            f.Position = [257   125   697   652];
            subplot(2,1,1)
            stem(params.pOnTimes_trains{i_tr}, ones(size(params.pOnTimes_trains{i_tr})))
            xlim([tt(1), tt(end)])
            box off
            
            subplot(2,1,2)
            isi = [0, diff(params.pOnTimes_trains{i_tr})];
            ifreq = 1./isi;
            plot(params.pOnTimes_trains{i_tr}, isi, '.-')
            xlim([tt(1), tt(end)])
            box off
            drawnow
            
        end
        
        
        for i_tr = 1:numel(params.templates_poiss)
            
            
            tt = (0:numel(params.templates_poiss{i_tr})-1) .* params.si;
            
            f = figure;
            f.Position = [257   125   697   652];
            subplot(2,1,1)
            stem(params.pOnTimes_poiss{i_tr}, ones(size(params.pOnTimes_poiss{i_tr})))
            xlim([tt(1), tt(end)])
            box off
            
            subplot(2,1,2)
            isi = [0, diff(params.pOnTimes_poiss{i_tr})];
            ifreq = 1./isi;
            plot(params.pOnTimes_poiss{i_tr}, isi, '.-')
            xlim([tt(1), tt(end)])
            box off
            drawnow
            
        end
    end
    
    
    %
    % CONFIRM FITTING CODE WORKS ON NOISELESS DATA
    %
    
    % things in common for both stimulus types
    Ntrains = numel(params.templates_trains); % be fair and only compare the same number of trials
    
    % construct a fake dataset from the trains
    testdat = [];
    for i_tr = 1:Ntrains
        ttype = sprintf('train_%d', i_tr);
        testdat.expt.(ttype).stats.EPSCamp{CHANNEL} = params.psc_test_trains{i_tr};
        testdat.expt.(ttype).pOnTimes = params.pOnTimes_trains{i_tr};
        testdat.expt.(ttype).raw.snips{CHANNEL} = [nan, nan]; % just needs to be defined inorder for fitTau2stp to run
    end
    
    [d_trains, f_trains, dTau_trains, fTau_trains] = fitTau2STP(testdat, PSCTYPE, CHANNEL, FITMETHOD);
    train_fits_noiseless(i_iter,:) = [d_trains, f_trains, dTau_trains, fTau_trains];
    
    
    % construct a fake dataset from the RIT stimuli
    testdat = [];
    for i_tr = 1:Ntrains;
        ttype = sprintf('RIT_%d', i_tr);
        testdat.expt.(ttype).stats.EPSCamp{CHANNEL} = params.psc_test_poiss{i_tr};
        testdat.expt.(ttype).pOnTimes = params.pOnTimes_poiss{i_tr};
        testdat.expt.(ttype).raw.snips{CHANNEL} = [nan, nan]; % just needs to be defined inorder for fitTau2stp to run
    end
    
    [d_rit, f_rit, dTau_rit, fTau_rit] = fitTau2STP(testdat, PSCTYPE, CHANNEL, FITMETHOD);
    rit_fits_noiseless(i_iter,:) = [d_rit, f_rit, dTau_rit, fTau_rit];
    
    %
    % ADD GAUSSIAN NOISE AND DETERMINE THE QUALITY OF FITS
    %
    
    % things in common for both stimulus types
    Ntrials = numel(params.templates_poiss) ./ numel(params.templates_trains)
    
    
    % construct a fake dataset from the trains, adding noise and making the
    % appropriate number of trials per stimulus
    testdat = [];
    for i_tr = 1:numel(params.templates_trains)
        ttype = sprintf('train_%d', i_tr);
        
        % pull out the noiseless pscs
        psc = params.psc_test_trains{i_tr};
        
        % make the noise have a fixed proportion of the (noiseless) versions
        sigma = psc ./ SIGMAFACTOR;
        psc = normrnd(repmat(psc, [1,1,Ntrials]), repmat(sigma, [1, 1, Ntrials]));
        if any(psc<0)
            psc(psc<0) = 0;
            fprintf('######### found one #######\n')
        end
        
        
        testdat.expt.(ttype).stats.EPSCamp{CHANNEL} = psc;
        testdat.expt.(ttype).pOnTimes = params.pOnTimes_trains{i_tr};
        testdat.expt.(ttype).raw.snips{CHANNEL} = [nan, nan]; % just needs to be defined inorder for fitTau2stp to run
    end
    
    [d_trains, f_trains, dTau_trains, fTau_trains] = fitTau2STP(testdat, PSCTYPE, CHANNEL, FITMETHOD);
    train_fits(i_iter,:) = [d_trains, f_trains, dTau_trains, fTau_trains];
    
    
    % construct a fake dataset from the RIT stimuli
    testdat = [];
    for i_tr = 1:numel(params.psc_test_poiss);
        
        ttype = sprintf('RIT_%d', i_tr);
        
        psc = params.psc_test_poiss{i_tr};
        sigma = psc ./ SIGMAFACTOR;
        noise = normrnd(psc, sigma);
        psc(psc<0) = 0;
        
        testdat.expt.(ttype).stats.EPSCamp{CHANNEL} = psc;
        testdat.expt.(ttype).pOnTimes = params.pOnTimes_poiss{i_tr};
        testdat.expt.(ttype).raw.snips{CHANNEL} = [nan, nan]; % just needs to be defined inorder for fitTau2stp to run
    end
    
    [d_rit, f_rit, dTau_rit, fTau_rit] = fitTau2STP(testdat, PSCTYPE, CHANNEL, FITMETHOD);
    rit_fits(i_iter,:) = [d_rit, f_rit, dTau_rit, fTau_rit];
    
end




%% PLOT THE RESULTS OF THE SIMULATION

% compile the noiseless and noise-added fitted values from trains and RITs
d_rit = rit_fits(:,1:2);
f_rit = rit_fits(:,3);
dTau_rit = rit_fits(:,4:5);
fTau_rit = rit_fits(:,6);

d_rit_noiseless = rit_fits_noiseless(:,1:2);
f_rit_noiseless = rit_fits_noiseless(:,3);
dTau_rit_noiseless = rit_fits_noiseless(:,4:5);
fTau_rit_noiseless = rit_fits_noiseless(:,6);

d_trains = train_fits(:,1:2);
f_trains = train_fits(:,3);
dTau_trains = train_fits(:,4:5);
fTau_trains = train_fits(:,6);

d_trains_noiseless = train_fits_noiseless(:,1:2);
f_trains_noiseless = train_fits_noiseless(:,3);
dTau_trains_noiseless = train_fits_noiseless(:,4:5);
fTau_trains_noiseless = train_fits_noiseless(:,6);

% fix the order of the depression constants (they can get fit in any order)
[~, minidx] = min(d_rit,[],2);
flipidx = minidx==2;
d_rit(flipidx,:) = fliplr(d_rit(flipidx,:));
dTau_rit(flipidx,:) = fliplr(dTau_rit(flipidx,:));

[~, minidx] = min(d_rit_noiseless, [], 2);
flipidx = minidx==2;
d_rit_noiseless(flipidx,:) = fliplr(d_rit_noiseless(flipidx,:));
dTau_rit_noiseless(flipidx,:) = fliplr(dTau_rit_noiseless(flipidx,:));

[~, minidx] = min(d_trains, [], 2);
flipidx = minidx==2;
d_trains(flipidx,:) = fliplr(d_trains(flipidx,:));
dTau_trains(flipidx,:) = fliplr(dTau_trains(flipidx,:));

[~, minidx] = min(d_trains_noiseless, [], 2);
flipidx = minidx==2;
d_trains_noiseless(flipidx,:) = fliplr(d_trains_noiseless(flipidx,:));
dTau_trains_noiseless(flipidx,:) = fliplr(dTau_trains_noiseless(flipidx,:));

%
% make two plots comparing the fitted values. Start with the noiseless fits
%
pclr = 'k';
hf = figure;
hf.Position = [221         -16        1033         683];
subplot(2,3,1), hold on,
plot(d_trains_noiseless(:,1), d_rit_noiseless(:,1), '.', 'color', pclr)
plot(mean(d_trains_noiseless(:,1)), mean(d_rit_noiseless(:,1)), 's', 'color', pclr, 'markerfacecolor', pclr, 'markersize', 8)
plot(d1_sim, d1_sim, 'ro', 'markerfacecolor', 'r', 'markersize', 8)
title('D1')

subplot(2,3,2), hold on,
plot(d_trains_noiseless(:,2), d_rit_noiseless(:,2), '.', 'color', pclr)
plot(mean(d_trains_noiseless(:,2)), mean(d_rit_noiseless(:,2)), 's', 'color', pclr, 'markerfacecolor', pclr, 'markersize', 8)
plot(d2_sim, d2_sim, 'ro', 'markerfacecolor', 'r', 'markersize', 8)
title('D2')

subplot(2,3,3), hold on,
plot(f_trains_noiseless, f_rit_noiseless, '.', 'color', pclr)
plot(mean(f_trains_noiseless), mean(f_rit_noiseless), 's', 'color', pclr, 'markerfacecolor', pclr, 'markersize', 8)
plot(f1_sim, f1_sim, 'ro', 'markerfacecolor', 'r', 'markersize', 8)
title('F1')

subplot(2,3,4), hold on,
plot(dTau_trains_noiseless(:,1), dTau_rit_noiseless(:,1), '.', 'color', pclr)
plot(mean(dTau_trains_noiseless(:,1)), mean(dTau_rit_noiseless(:,1)), 's', 'color', pclr, 'markerfacecolor', pclr, 'markersize', 8)
plot(tau_d1_sim, tau_d1_sim, 'ro', 'markerfacecolor', 'r', 'markersize', 8)
title('Tau D1')

subplot(2,3,5), hold on,
plot(dTau_trains_noiseless(:,2), dTau_rit_noiseless(:,2), '.', 'color', pclr)
plot(mean(dTau_trains_noiseless(:,2)), mean(dTau_rit_noiseless(:,2)), 's', 'color', pclr, 'markerfacecolor', pclr, 'markersize', 8)
plot(tau_d2_sim, tau_d2_sim, 'ro', 'markerfacecolor', 'r', 'markersize', 8)
title('Tau D2')

subplot(2,3,6), hold on,
plot(fTau_trains_noiseless, fTau_rit_noiseless, '.', 'color', pclr)
plot(mean(fTau_trains_noiseless), mean(fTau_rit_noiseless), 's', 'color', pclr, 'markerfacecolor', pclr, 'markersize', 8)
plot(tau_f1_sim, tau_f1_sim, 'ro', 'markerfacecolor', 'r', 'markersize', 8)
title('Tau F1')

for i_plt = 1:6;
    subplot(2,3,i_plt)
    minvals = min([get(gca, 'xlim'), get(gca, 'ylim')]);
    maxvals = max([get(gca, 'xlim'), get(gca, 'ylim')]);
    plot([minvals; maxvals], [minvals; maxvals], 'k:')
    axis tight
    xlabel('Recovery Trains')
    ylabel('RITs')
end



%
%  Now plot the simulation with added noise
%
f = figure;
f.Position = [221         -16        1033         683];
subplot(2,3,1), hold on,
plot(d_trains(:,1), d_rit(:,1), '.', 'color', pclr)
plot(mean(d_trains(:,1)), mean(d_rit(:,1)), 's', 'color', pclr, 'markerfacecolor', pclr, 'markersize', 8)
plot(d1_sim, d1_sim, 'ro', 'markerfacecolor', 'r', 'markersize', 8)
title('D1')

subplot(2,3,2), hold on,
plot(d_trains(:,2), d_rit(:,2), '.', 'color', pclr)
plot(mean(d_trains(:,2)), mean(d_rit(:,2)), 's', 'color', pclr, 'markerfacecolor', pclr, 'markersize', 8)
plot(d2_sim, d2_sim, 'ro', 'markerfacecolor', 'r', 'markersize', 8)
title('D2')

subplot(2,3,3), hold on,
plot(f_trains, f_rit, '.', 'color', pclr)
plot(mean(f_trains), mean(f_rit), 's', 'color', pclr, 'markerfacecolor', pclr, 'markersize', 8)
plot(f1_sim, f1_sim, 'ro', 'markerfacecolor', 'r', 'markersize', 8)
title('F1')

subplot(2,3,4), hold on,
plot(dTau_trains(:,1), dTau_rit(:,1), '.', 'color', pclr)
plot(mean(dTau_trains(:,1)), mean(dTau_rit(:,1)), 's', 'color', pclr, 'markerfacecolor', pclr, 'markersize', 8)
plot(tau_d1_sim, tau_d1_sim, 'ro', 'markerfacecolor', 'r', 'markersize', 8)
title('Tau D1')

subplot(2,3,5), hold on,
plot(dTau_trains(:,2), dTau_rit(:,2), '.', 'color', pclr)
plot(mean(dTau_trains(:,2)), mean(dTau_rit(:,2)), 's', 'color', pclr, 'markerfacecolor', pclr, 'markersize', 8)
plot(tau_d2_sim, tau_d2_sim, 'ro', 'markerfacecolor', 'r', 'markersize', 8)
title('Tau D2')

subplot(2,3,6), hold on,
plot(fTau_trains, fTau_rit, '.', 'color', pclr)
plot(mean(fTau_trains), mean(fTau_rit), 's', 'color', pclr, 'markerfacecolor', pclr, 'markersize', 8)
plot(tau_f1_sim, tau_f1_sim, 'ro', 'markerfacecolor', 'r', 'markersize', 8)
title('Tau F1')


for i_plt = 1:6;
    subplot(2,3,i_plt)
    minvals = min([get(gca, 'xlim'), get(gca, 'ylim')]);
    maxvals = max([get(gca, 'xlim'), get(gca, 'ylim')]);
    plot([minvals; maxvals], [minvals; maxvals], 'k:')
    axis tight
    xlabel('Recovery Trains')
    ylabel('RITs')
end


%% TEST CODE FOR STP TIME CONSTANTS (2)

%
% estimate the shape of the error function
%


pOnTimes = params.pOnTimes_poiss{1};
pscreal = predictPSCfromTau(pOnTimes, [d1_sim, d2_sim], [tau_d1_sim, tau_d2_sim], f1_sim, tau_f1_sim, A0);

[a,b] = ndgrid(0.001:0.005:1, 0:0.050:3);
errvals = nan(size(a));
for idx = 1:numel(a)
        k_d = [a(idx), d2_sim];
        tau_d = [tau_d1_sim, tau_d2_sim];
        k_f = b(idx);
        tau_f = tau_f1_sim;

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

