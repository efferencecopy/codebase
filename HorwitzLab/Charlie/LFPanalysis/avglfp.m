function [lfp, timeVec, shufflfp] = avglfp(obj, lfpchan, tList, params)

%params should have the following fields:
%   .trig       => an index into the .trial field (for alligning the data)
%   .notch.f1   => the lower edge of the band stop
%   .notch.f2   => the upper edge of the band stop
%   .notch.deg  => the filter deg
%   .detrend.on => if set to 1, then detrend
%   .detrend.win => the window for linear removal of slow trends
%   .pretime    => seconds b/4 trig to consider
%   .posttime   => seconds after trig to consider
%   .bootstraps => the number of shuffle predictors to calculate
%   .average    => default is to average over trials. Set to zero otherwise
%   .zscore     => convert time domain LFP into z scores


% set the defaults (if need be)
if ~exist('params', 'var')
    %notch filter
    params.notch.deg = 3; %filter deg
    params.notch.f1 = 59; %low edge of the bandstop
    params.notch.f2 = 61; %upper edge of the bandstop
    
    %detrending
    params.detrend.on = 1;
    params.detrend.win = [.500 0.100];
    
    %pre/post time
    params.pretime = 0;
    params.posttime = 0.500;
    
    %miscalenous
    params.bootstraps = 0;
    params.average = 1;
    params.zscore = 0;
    
    %tailor the trig event to the paradigm type
    switch obj.sum.paradigmID
        case 210
            params.trig = obj.idx.flashOn;
        case 150
            params.trig = obj.idx.stimon;
        case 100
            params.trig = obj.idx.stim_on;
        otherwise
            error('unknown paradigm type')
    end
end
    
%set up the notch filter
anlgRate = obj.sum.analog.storeRates{lfpchan==find(obj.idx.lfp)};
Wn = [params.notch.f1*2/anlgRate, params.notch.f2*2/anlgRate];
[B, A] = butter(params.notch.deg, Wn, 'stop');




%iterate over trials and process the analog data
preSamples = params.pretime * anlgRate;
postSamples = params.posttime * anlgRate;
lfp = nan(length(tList), preSamples+postSamples+1);
[shufflfp{1:length(tList)}] = deal(nan(params.bootstraps, preSamples+postSamples+1));
for a = 1:length(tList)
    trigTime = obj.trial(tList(a), params.trig);
    anlgStart = obj.ras{tList(a), obj.idx.anlgStart};
    idxAtTrig = round((trigTime-anlgStart)*anlgRate + 1);
    triallfp = obj.ras{tList(a), lfpchan};
    
    %notch filter
    triallfp = filtfilt(B, A, triallfp);
    
    %detrend (if need be)
    if params.detrend.on
        triallfp = locdetrend(triallfp', anlgRate, params.detrend.win);
    end
    
    % Zscore if need be
    if params.zscore
        triallfp = (triallfp-mean(triallfp)) ./ std(triallfp);
    end
    
    if params.absval
        triallfp = abs(triallfp);
    end
    
    
    % select the appropriate samples and store in an array
    if (idxAtTrig-preSamples)>=1 && (idxAtTrig+postSamples)<=numel(triallfp);
        lfp(a,:) = triallfp(idxAtTrig-preSamples:idxAtTrig+postSamples);
    end
    
    %shuffle the lfp if need be
    if params.bootstraps
        for i = 1:params.bootstraps
            idx = round(unifrnd(preSamples+1, length(triallfp)-postSamples-1));
            shufflfp{a}(i,:) = triallfp(idx-preSamples:idx+postSamples);
        end
    end
   
end

%the main output args.
if params.average
    lfp = nanmean(lfp, 1);
end
timeVec = -params.pretime:1/anlgRate:params.posttime;

%the shuffled lfp
if params.bootstraps
    shufflfp = vertcat(shufflfp{:});
    shufflfp = mean(shufflfp,1);
else
    shufflfp = [];
end



