%% MAKE SOME PULSE TRAINS: RECOVERY, RIT, ENVELOPED-RIT

fin

PLOTFIGS = false;
NONSTATIONARY = true;
FITMETHOD = 'multistart';

SIGMAFACTOR = 5;  % noise is 1/sigmaFactor of the PSC amplitude. This is also the Coeff. of Var
Niters = 200;
Ntrials = 2;

% simulation params
d1_sim = 0.60;
d2_sim = 0.90;
tau_d1_sim = 0.100;
tau_d2_sim = 2;
f1_sim = 0.6;
tau_f1_sim = 0.350;
A0_sim = 600;

[rit_fits, rit_fits_noiseless, train_fits, train_fits_noiseless] = deal(nan(Niters, 6));
[R2_rit_pred_from_train_fits, R2_train_pred_from_rit_fits] = deal(nan(Niters, 1));

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
    
    %warning('only using 25 RIT versions')
    
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
        if NONSTATIONARY
            A0 = A0.*0.95;
        end
        
        params.psc_test_trains{i_tr} = predictPSCfromTau(params.pOnTimes_trains{i_tr}, d, dTau, f, fTau, A0);
        
        if PLOTFIGS
            figure
            plot(params.pOnTimes_trains{i_tr}, params.psc_test_trains{i_tr}, '-o')
            xlim([tt(1), tt(end)])
        end
        
    end
    
    A0 = A0_sim;
    for i_tr = 1:numel(params.templates_poiss)
        
        % define the pOnTimes
        tt = (0:numel(params.templates_poiss{i_tr})-1) .* params.si;
        crossing = params.templates_poiss{i_tr} > 0.25;
        crossing = [0; diff(crossing)] == 1;
        params.pOnTimes_poiss{i_tr} = tt(crossing);
        
        %simulate the EPSCs
        if NONSTATIONARY
            A0 = A0.*0.95;
        end
        
        params.psc_test_poiss{i_tr} = predictPSCfromTau(params.pOnTimes_poiss{i_tr}, d, dTau, f, fTau, A0);
        
        if PLOTFIGS
            figure
            plot(params.pOnTimes_poiss{i_tr}, params.psc_test_poiss{i_tr}, '-o')
            xlim([tt(1), tt(end)])
        end
        
    end
    
    %
    % OPTIONAL PLOTS OF THE TEMPLATES AND INSTANTANEOUS FREQ
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
    
    % construct a fake dataset from the trains I need the pOnTimes, raw EPSC amps.
    [pOnTimes, rawAmps] = deal({});
    for i_tr = 1:Ntrains
        pOnTimes{i_tr} =  params.pOnTimes_trains{i_tr};
        rawAmps{i_tr} = params.psc_test_trains{i_tr};
    end
    
    [d_trains, f_trains, dTau_trains, fTau_trains] = fitTau2STP(rawAmps, pOnTimes, FITMETHOD);
    train_fits_noiseless(i_iter,:) = [d_trains, f_trains, dTau_trains, fTau_trains];
    
    
    % construct a fake dataset from the RIT stimuli
    [pOnTimes, rawAmps] = deal({});
    for i_tr = 1:numel(params.psc_test_poiss);
        pOnTimes{i_tr} =  params.pOnTimes_poiss{i_tr};
        rawAmps{i_tr} = params.psc_test_poiss{i_tr};
    end
    
    [d_rit, f_rit, dTau_rit, fTau_rit] = fitTau2STP(rawAmps, pOnTimes, FITMETHOD);
    rit_fits_noiseless(i_iter,:) = [d_rit, f_rit, dTau_rit, fTau_rit];
    
    %
    % ADD GAUSSIAN NOISE AND DETERMINE THE QUALITY OF FITS
    %
    
    % things in common for both stimulus types
    %Ntrials = numel(params.templates_poiss) ./ numel(params.templates_trains)
    
    
    % construct a fake dataset from the trains, adding noise and making the
    % appropriate number of trials per stimulus
    [pOnTimes_trains, noisyAmps_trains] = deal({});
    for i_tr = 1:numel(params.psc_test_trains)
        
        % pull out the noiseless pscs
        psc = params.psc_test_trains{i_tr};
        
        % make the noise have a fixed proportion of the (noiseless) versions
        sigma = psc ./ SIGMAFACTOR;
        psc = normrnd(repmat(psc, [1,1,Ntrials]), repmat(sigma, [1, 1, Ntrials]));
        assert(all(psc(:)>=0), '######### found one #######\n')
        
        noisyAmps_trains{i_tr} = psc;
        pOnTimes_trains{i_tr} = params.pOnTimes_trains{i_tr};
    end
    
    [d_trains, f_trains, dTau_trains, fTau_trains] = fitTau2STP(noisyAmps_trains, pOnTimes_trains, FITMETHOD);
    train_fits(i_iter,:) = [d_trains, f_trains, dTau_trains, fTau_trains];
    
    
    % construct a fake dataset from the RIT stimuli
    [pOnTimes_RIT, noisyAmps_rit] = deal({});
    for i_tr = 1:numel(params.psc_test_poiss);
        
        psc = params.psc_test_poiss{i_tr};
        sigma = psc ./ SIGMAFACTOR;
        psc = normrnd(psc, sigma);
        assert(all(psc(:)>=0), '######### found one #######\n')
        
        noisyAmps_rit{i_tr} = psc;
        pOnTimes_RIT{i_tr} = params.pOnTimes_poiss{i_tr};
    end
    
    [d_rit, f_rit, dTau_rit, fTau_rit] = fitTau2STP(noisyAmps_rit, pOnTimes_RIT, FITMETHOD);
    rit_fits(i_iter,:) = [d_rit, f_rit, dTau_rit, fTau_rit];
    
    %
    % use the fits from one stimulus type to cross validate the other type.
    % Start with using the train fits to cross validate the RIT
    %
    
    % use train fits, but the RIT pOnTime
    A0 = cellfun(@(x) x(1), noisyAmps_rit);
    A0 = mean(A0);
    rit_Pred_From_Train_Fits={};
    for i_tr = 1:numel(pOnTimes_RIT)
        rit_Pred_From_Train_Fits{i_tr} = predictPSCfromTau(pOnTimes_RIT{i_tr}, d_trains, dTau_trains, f_trains, fTau_trains, A0);
    end
    resid = cellfun(@(x,y) x-y, rit_Pred_From_Train_Fits, params.psc_test_poiss, 'uniformoutput', false);
    resid = cat(1, resid{:});
    SS_resid = sum(resid.^2);
    SS_raw = sum(cat(1, params.psc_test_poiss{:}).^2);
    R2_rit_pred_from_train_fits(i_iter) = 1 - (SS_resid ./ SS_raw);
    
    % use rit fits, but the Train pOnTime
    A0 = cellfun(@(x) x(1), noisyAmps_trains);
    A0 = mean(A0);
    train_Pred_From_RIT_Fits={};
    for i_tr = 1:numel(pOnTimes_trains)
        train_Pred_From_RIT_Fits{i_tr} = predictPSCfromTau(pOnTimes_trains{i_tr}, d_rit, dTau_rit, f_rit, fTau_rit, A0);
    end
    resid = cellfun(@(x,y) x-y, train_Pred_From_RIT_Fits, params.psc_test_trains, 'uniformoutput', false);
    resid = cat(1, resid{:});
    SS_resid = sum(resid.^2);
    SS_raw = sum(cat(1, params.psc_test_trains{:}).^2);
    R2_train_pred_from_rit_fits(i_iter) = 1 - (SS_resid ./ SS_raw);
    
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

%% LOOK AT PREDICTIONS AND CROSS VALIDATION

% train the fits using the RIT data, then plot the fits to all the data,
% and calculate Rsquared values.


% construct a fake, noisy dataset from the RIT stimuli
[pOnTimes, rawAmps_noisyRIT] = deal({});
for i_tr = 1:numel(params.psc_test_poiss);
    
    psc = params.psc_test_poiss{i_tr};
    sigma = psc ./ SIGMAFACTOR;
    psc = normrnd(psc, sigma);
    assert(all(psc>=0), 'found one')
    
    rawAmps_noisyRIT{i_tr} = psc;
    pOnTimes{i_tr} = params.pOnTimes_poiss{i_tr};
    A0_noisy(i_tr) = psc(1);
end
A0_noisy = mean(A0_noisy);

[d_rit, f_rit, dTau_rit, fTau_rit] = fitTau2STP(rawAmps_noisyRIT, pOnTimes, FITMETHOD);


% look to see how well the fitted params estimate the true, but noise data
for i_type = 1:numel(params.pOnTimes_poiss)
    
    % use the fitted params to make predictions about the simulated noisy
    % data
    pOnTimes = params.pOnTimes_poiss{i_type};
    rit_fitToNoisyData = predictPSCfromTau(pOnTimes, d_rit, dTau_rit, f_rit, fTau_rit, A0_noisy);
    
    % plot actual noisy data vs. predicted noisy data. 
    rit_realNoisyData = rawAmps_noisyRIT{i_type};
    
    figure, hold on,
    stem(1:numel(rit_fitToNoisyData), rit_realNoisyData, 'b')
    plot(1:numel(rit_fitToNoisyData), rit_fitToNoisyData, 'k')
    plot(1:numel(rit_fitToNoisyData), params.psc_test_poiss{i_type}, '--k')
    legend('noisy ''real'' data', 'fit to noisy data', 'noise-free data')
end


% compare to the noiseless, but un-fit trains data
[pred, raw_xbar] = deal({});
for i_type = 1:numel(params.pOnTimes_trains)
    pOnTimes = params.pOnTimes_trains{i_type};
    raw_xbar{i_type} = predictPSCfromTau(pOnTimes, [d1_sim, d2_sim], [tau_d1_sim, tau_d2_sim], f1_sim, tau_f1_sim, A0_sim);
    A0 = raw_xbar{i_type}(1);
    pred{i_type} = predictPSCfromTau(pOnTimes, d_rit, dTau_rit, f_rit, fTau_rit, A0);
end


% plot the recovery trains (cross validation)
hf = figure;
xlims = [inf -inf];
ylims = [inf -inf];
hs = [];
for i_type = 1:numel(params.pOnTimes_trains)
    i_ch = 1;
    pltIdx = sub2ind([2, numel(params.pOnTimes_trains)], i_ch+1, i_type);
    hs(i_type) = subplot(numel(params.pOnTimes_trains), 2, pltIdx); hold on,
    
    xx = params.pOnTimes_trains{i_type};
    plot(xx, raw_xbar{i_type}, 'ko');
    plot(xx, pred{i_type}, 'r', 'linewidth', 2)
    xlims(1) = min([min(xx), xlims(1)]);
    xlims(2) = max([max(xx), xlims(2)]);
    yvals = get(gca, 'ylim');
    ylims(1) = min([yvals(1), ylims(1)]);
    ylims(2) = max([yvals(2), ylims(2)]);
    if i_type == 1; title('Train on noisy RIT, pred noiseless trains'); end
end
set(hs, 'XLim', xlims, 'YLim', ylims)



% make a scatter plot of all predicted and actual PSC amps from the noisy
% RITs
hs = [];
all_raw = [];
all_pred = [];
crossval_raw = [];
crossval_pred = [];
for i_type = 1:numel(params.pOnTimes_poiss)
   
    tmp_raw = dat{i_ex}.expt.(trainNames{i_type}).stats.EPSCamp{i_ch};
    tmp_pred = pred{i_ch}{i_type};
    tmp_pred = repmat(tmp_pred, size(tmp_raw,3), 1);
    tmp_raw = tmp_raw(:);
    assert(all(size(tmp_pred) == size(tmp_raw)))
    
    pltIdx = sub2ind([4, 3], pltcol, 1);
    hs = subplot(3, 4, pltIdx); hold on,
    plot(tmp_raw./tmp_raw(1), tmp_pred./tmp_pred(1), 'k.')
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
