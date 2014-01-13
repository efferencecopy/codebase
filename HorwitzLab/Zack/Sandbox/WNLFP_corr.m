% Absolute pin number bank A
% 99 96 93 89 86 83 81 77 74 71 70 67 64 60 58 55 52 48 45 42 39 36 33 29 26 23 21 17 14 11 7 4
% Relative pin number within bank A
% 6 16 26 5 15 25 32 12 22 31 2 11 21 1 8 18 28 7 17 27 4 14 24 3 13 23 30 10 20 29 9 19
% Absolute pin number bank B
% 98 95 92 88 85 82 79 76 73 69 66 63 61 57 54 51 50 47 44 38 35 32 28 25 22 19 16 13 10 9 6 3
% Relative pin number within bank B
% 9 19 29 10 20 30 5 15 25 6 16 26 31 11 21 32 1 12 22 7 17 27 8 18 28 3 13 23 2 4 14 24
% Absolute pin number bank C
% 97 94 90 87 84 80 78 75 72 68 65 62 59 56 53 49 46 43 41 37 34 31 30 27 24 20 18 15 12 8 5 2
% Relative pin number within bank C
% 14 24 4 13 23 3 10 20 30 9 19 29 6 16 26 5 15 25 32 12 22 31 2 11 21 1 8 18 28 7 17 27
%
% The 10-by-10 array is oriented like [   ]-   where "-" represents the wire, and
% the absolute numbering is flipud(reshape(1:100,[10 10])'). Pins 1, 40, 91, and 100 are inactive.

%% 1 - Attempting to plot correlations (between pins) as a function of pin separation and time lags
% Uses gun noise trials only

WN = wnobj(findfile('A051311002.nex'));
bank = 'B'; % the most lateral bank is "A", in-between bank "B", most medial bank "C"
GUNNOISE = 1; CONENOISE = 2;
noise_idxs = WN.trial(:,WN.idx.noise_type) == CONENOISE;

LFPsamprate = WN.sum.analog.storeRates{1};  % Assuming all channels are sampled at the same rate
LFPchans = find(WN.idx.lfp);
LFPstarttimes = cell2mat(WN.ras(:,strcmp(WN.sum.rasterCells, 'anlgStartTime')));
LFPstarttimes = LFPstarttimes(noise_idxs); % gun noise only

poststim_t = 0.2; % secs
prestim_t = 0.2; % secs
maxlag = ceil(LFPsamprate .* 0.2);
nprestimsamples = prestim_t*LFPsamprate;
npoststimsamples = poststim_t*LFPsamprate;
nsamples = nprestimsamples+npoststimsamples+1;

WNtrial = WN.trial(noise_idxs,:);
WNras = WN.ras(noise_idxs,:);
ntrials = size(WNtrial,1);

% % xcorr for the average PSTH
% zPSTHs = zeros(nsamples,length(LFPchans));
% for i = 1:ntrials
%     t_stimon = WNtrial(i,WN.idx.stim_on);
%     for chan = LFPchans
%         LFP = zscore(WNras{i,chan});
%         LFPtimes = LFPstarttimes(i)+(0:length(LFP)-1)/LFPsamprate; % sec
%         err = (LFPtimes-t_stimon).^2;
%         startidx = find(err == min(err),1);
%         if startidx > nprestimsamples && length(LFP)-startidx > npoststimsamples
%             LFPclip = LFP(startidx+(-nprestimsamples:npoststimsamples));
%             zPSTHs(:,chan == LFPchans) = zPSTHs(:,chan == LFPchans) + abs(LFPclip);
%         end
%     end
% end
% 
% [zcorrs,lags] = xcorr(zPSTHs,maxlag); % columns of zcorrs: [corr(lfp1,lfp1) corr(lfp1,lfp2) ... corr(lfp1,lfp32) corr(lfp2,lfp1) ... corr(lfp32,lfp32)]
% dists = distance_between_pins(pins_rel2abs(bank, 1:32)); % going to assume the order of AD channels maps directly to pin number
% % if this turns out not to be the case, then reorder the second argument to reflect reality
% % clear out auto correlations?
% Lautos = dists == 0;
% zcorrs(:,Lautos) = [];
% dists(Lautos) = [];
% 
% figure;
% h = surf(dists,lags,zcorrs);
% set(h,'edgecolor','none'); colormap(pmkmp([],'cubicl')); colorbar;
% xlabel('Normalized pin distance'); ylabel('Lag'); zlabel('Correlation');

% now trying to do xcorr on a trial by trial basis.
avgxcorr = zeros(2*maxlag+1, length(LFPchans)^2);
nValidTrials = 0;
for i = 1:ntrials
    trialLFP = zeros(nsamples, length(LFPchans));
    t_stimon = WNtrial(i,WN.idx.stim_on);
    for chan = LFPchans
        
        LFP = zscore(WNras{i,chan});
        LFPtimes = LFPstarttimes(i)+(0:length(LFP)-1)/LFPsamprate; % sec
        err = (LFPtimes-t_stimon).^2;
        startidx = find(err == min(err),1);        
        
        % look for oob errors
        nADSteps = 2^12;
        anlgChannelName = WN.sum.rasterCells{chan};
        anlgChannelIdx = strcmpi(WN.sum.analog.sigid, anlgChannelName);
        ADtoMV = WN.sum.analog.ADtoMV{anlgChannelIdx};
        maxVal = (ADtoMV * nADSteps) ./ 2;
        if any(abs(WNras{i,chan}) >= maxVal) % oobs
            continue
        end        
        
        if startidx > nprestimsamples && length(LFP)-startidx > npoststimsamples
            trialLFP(:,chan == LFPchans) = LFP(startidx+(-nprestimsamples:npoststimsamples));
            if chan == LFPchans(1)
                nValidTrials = nValidTrials + 1;
            end
        end
    end
    
    % now compute the xcorr
    avgxcorr = avgxcorr + xcorr(trialLFP, maxlag);
end

% divide by the number of trials
avgxcorr = avgxcorr / nValidTrials;
lags = (-maxlag : maxlag) / LFPsamprate;

dists = distance_between_pins(pins_rel2abs(bank, 1:32)); % going to assume the order of AD channels maps directly to pin number
% if this turns out not to be the case, then reorder the second argument to reflect reality

% % clear out auto correlations?
Lautos = dists == 0;
avgxcorr(:,Lautos) = [];
dists(Lautos) = [];

% % all the data
% figure, hold on,
% for a = 1:size(avgxcorr,2)
%     plot3(repmat(dists(a), numel(lags),1), lags, avgxcorr(:,a), 'k')
% end
% xlabel('Normalized pin distance'); ylabel('Lag'); zlabel('Correlation');


% now binning "distances"
uniqueDists = unique(dists);
xcorrAvgOverDist = zeros(numel(uniqueDists),numel(lags));
for a = 1:numel(uniqueDists)
    idx = dists == uniqueDists(a);
    xcorrAvgOverDist(a,:) = mean(avgxcorr(:, idx), 2);
end

figure;
h = surf(lags, uniqueDists, xcorrAvgOverDist);
set(h,'edgecolor','none'); colormap(pmkmp([],'cubicl')); colorbar;
xlabel('Lag (sec)'); ylabel('Pin distance'); zlabel('Correlation');
