%% MAKE SOME PULSE TRAINS: RECOVERY, RIT, ENVELOPED-RIT

fin

% simulation params
Ntrials = 3;
d1_sim = 0.7;
d2_sim = 0.95;
tau_d1_sim = 0.500;
tau_d2_sim = 1;
f1_sim = 0.6;
tau_f1_sim = 0.5;
A0_sim = 600;

rit_fits = [];
train_fits = [];

for i_iter = 1:100;
    
    i_iter
%
% Define the params that influence the entire data file (trains and RITs)
%
params.si   = 50e-6;              % the sample INTERVAL (needs to be an iteger)
params.swpDur = 20e4;             % The total duration of the sweep IN NUMBERS OF SAMPLES!!!!
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
params.recoveryTime = [0.250, .500, 2, 8];   %  A vector of numbers corresponding to the recovery time in seconds [interleaved variable]
params = makeSweepTemplates_trains(params); % templates are stored in params.templates_trains


%
% make the sweep template(s) for the Random impulse trains without an envelope
%
params.ritFreq = [8];
params.ritHiFreqCut = 58;  % ISIs faster than this will be cutout
params.rit_Nversions = numel(params.templates_trains) .* Ntrials;
params.ritUseEnvelope = true;
params.ritEnvelopeFreq = [0.25];
params = makeSweepTemplates_poiss(params); % templates are stored in params.templates_poiss



% concatenate templates
sweepTemplates = cat(2, params.templates_trains, params.templates_poiss);


% GENERATE NOISELESS SIMULATED RESPONSES

PLOTFIG = true;

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
    
    if PLOTFIG
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
    
    if PLOTFIG
        figure
        plot(params.pOnTimes_poiss{i_tr}, params.psc_test_poiss{i_tr}, '-o')
        xlim([tt(1), tt(end)])
    end
    
end

% %% PLOT THE TEMPLATES AND INSTANTANEOUS FREQ
% 
% % start with recovery trains
% for i_tr = 1:numel(params.templates_trains)
%     
%     tt = (0:numel(params.templates_trains{i_tr})-1) .* params.si;
%     
%     f = figure;
%     f.Position = [257   125   697   652];
%     subplot(2,1,1)
%     stem(params.pOnTimes_trains{i_tr}, ones(size(params.pOnTimes_trains{i_tr})))
%     xlim([tt(1), tt(end)])
%     box off
%     
%     subplot(2,1,2)
%     isi = [0, diff(params.pOnTimes_trains{i_tr})];
%     ifreq = 1./isi;
%     plot(params.pOnTimes_trains{i_tr}, isi, '.-')
%     xlim([tt(1), tt(end)])
%     box off
%     drawnow
%     
% end
% 
% 
% for i_tr = 1:numel(params.templates_poiss)
%     
%     
%     tt = (0:numel(params.templates_poiss{i_tr})-1) .* params.si;
%     
%     f = figure;
%     f.Position = [257   125   697   652];
%     subplot(2,1,1)
%     stem(params.pOnTimes_poiss{i_tr}, ones(size(params.pOnTimes_poiss{i_tr})))
%     xlim([tt(1), tt(end)])
%     box off
%     
%     subplot(2,1,2)
%     isi = [0, diff(params.pOnTimes_poiss{i_tr})];
%     ifreq = 1./isi;
%     plot(params.pOnTimes_poiss{i_tr}, isi, '.-')
%     xlim([tt(1), tt(end)])
%     box off
%     drawnow
%     
% end

% 
% %% ESTIMATE TAUS TO CONFIRM FITTING CODE WORKS. NOISELESS DATA
% 
% % things in common for both stimulus types
% channel = 1;
% psctype = 'EPSCamp';
% method = 'multistart';
% Ntrains = numel(params.templates_trains); % be fair and only compare the same number of trials
% 
% % construct a fake dataset from the trains
% testdat = [];
% for i_tr = 1:Ntrains
%     ttype = sprintf('train_%d', i_tr);
%     testdat.expt.(ttype).stats.EPSCamp{channel} = params.psc_test_trains{i_tr};
%     testdat.expt.(ttype).pOnTimes = params.pOnTimes_trains{i_tr};
% end
% 
% [d_trains, f_trains, dTau_trains, fTau_trains] = fitTau2STP(testdat, psctype, channel, method);
% 
% 
% 
% % construct a fake dataset from the RIT stimuli
% testdat = [];
% for i_tr = 1:Ntrains;
%     ttype = sprintf('RIT_%d', i_tr);
%     testdat.expt.(ttype).stats.EPSCamp{channel} = params.psc_test_poiss{i_tr};
%     testdat.expt.(ttype).pOnTimes = params.pOnTimes_poiss{i_tr};
% end
% 
% [d_rit, f_rit, dTau_rit, fTau_rit] = fitTau2STP(testdat, psctype, channel, method);
% 
% 

% ESTIMATE TAUS TO CONFIRM FITTING CODE WORKS. GAUSSIAN NOISE

% things in common for both stimulus types
channel = 1;
psctype = 'EPSCamp';
method = 'multistart';
Ntrials = numel(params.templates_poiss) ./ numel(params.templates_trains);

sigmafactor = 10;

assert(rem(Ntrials, 1)==0);

% construct a fake dataset from the trains, adding noise and making the
% appropriate number of trials per stimulus
testdat = [];
for i_tr = 1:numel(params.templates_trains)
    ttype = sprintf('train_%d', i_tr);
    
    % pull out the noiseless pscs
    psc = params.psc_test_trains{i_tr};
    
    % make the noise have a fixed proportion of the (noiseless) versions
    sigma = psc ./ sigmafactor;
    psc = normrnd(repmat(psc, [1,1,Ntrials]), repmat(sigma, [1, 1, Ntrials]));
    assert(~any(psc(:)<0));
    
    testdat.expt.(ttype).stats.EPSCamp{channel} = psc;
    testdat.expt.(ttype).pOnTimes = params.pOnTimes_trains{i_tr};
end

[d_trains, f_trains, dTau_trains, fTau_trains] = fitTau2STP(testdat, psctype, channel, method);
train_fits(i_iter,:) = [d_trains, f_trains, dTau_trains, fTau_trains];


% construct a fake dataset from the RIT stimuli
testdat = [];
for i_tr = 1:numel(params.psc_test_poiss);
    
    ttype = sprintf('RIT_%d', i_tr);
    
    psc = params.psc_test_poiss{i_tr};
    sigma = psc ./ sigmafactor;
    noise = normrnd(psc, sigma);
    psc = noise + psc;
    assert(~any(psc(:)<0));
    
    testdat.expt.(ttype).stats.EPSCamp{channel} = psc;
    testdat.expt.(ttype).pOnTimes = params.pOnTimes_poiss{i_tr};
end

[d_rit, f_rit, dTau_rit, fTau_rit] = fitTau2STP(testdat, psctype, channel, method);
rit_fits(i_iter,:) = [d_rit, f_rit, dTau_rit, fTau_rit];

end









