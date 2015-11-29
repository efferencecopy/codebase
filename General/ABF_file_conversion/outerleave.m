function tDict = outerleave(stimWF, sampRate)
%
%  EXAMPLE    tDict = outerleave(ax, ledIdx)
% 
%  tDict.vars     => {'pAmp', 'pWidth', 'tFreq', 'tRecov'}  a reference to colums in tDict.conds
%  tDict.conds    => matrix of values. each ROW corresponds to a unique type
%  tDict.trlList  => matrix of scalars that map each sweep onto a 'condition', one column for each channel recorded
%
%  ** Be advised that this function will round pulse widths to the nearest
%     100 usec, becuse there is jitter in blacrock's timing of events
%     relative to how clampex aquires signals.


% condition the inputs if the data come from a .abf file
if ~iscell(stimWF)
    stimWF = permute(stimWF, [1,3,2]);
    stimWF = mat2cell(stimWF, size(stimWF, 1), ones(1,size(stimWF,2)));
end


Nsweeps = numel(stimWF);
si = 1./sampRate; % the sample interval...

% iterate over the sweeps determining the pulse amplitude, width, freq
[pWidth, pAmp, tFreq] = deal(nan(Nsweeps, 1));
tRecov = zeros(Nsweeps, 1); % using a numeric value as the default b/c nans will make 'unique' bonk later in the function
for swp = 1:Nsweeps
    
    tmpWF = stimWF{swp};
    tmpWF = tmpWF(:);
    
    % amplitude first (which sets the threshold)
    pAmp(swp) = max(tmpWF);
    pAmp(swp) = round(pAmp(swp).*100) ./ 100; % round to the one hundreths place
    
    % now find pulse width and IPI
    thresh = pAmp(swp) .* 0.8;
    above = tmpWF > thresh;
    change = [0; diff(above)];
    xUp = change == 1;
    xDown = change == -1;
    
    tt = [0:numel(xUp)-1] ./ sampRate;
    pOnTimes = tt(xUp);
    pOffTimes = tt(xDown);
    
    
    tmpWidth = mean(pOffTimes-pOnTimes) - si; % need to subtract one due to the OBO error introduced by the thresholding procedure
    tmpWidth = tmpWidth .* 1e4; % in hundreds of usec.
    tmpWidth = round(tmpWidth) ./ 1e4;
    assert(si<=50e-6, 'ERROR: pulse widths may be unreliable')
    pWidth(swp) = tmpWidth;
    
    
    if numel(pOnTimes)>1
        tFreq(swp) =  1./(pOnTimes(2)-pOnTimes(1));
        tFreq(swp) = round(tFreq(swp)); % round to the tenth place
        
        % was there a recovery pulse?
        lastIPI = round([pOnTimes(end)-pOnTimes(end-1)]*1000); % in milliseconds
        if lastIPI > ((1./tFreq(swp))+0.025)*1000 % needs to be 25 ms longer than the typical IPI
            tRecov(swp) = lastIPI;
        end
        
    else
        tFreq(swp) = 0; % a single pulse
    end
    
end


% now determine the number of unique trial types
tDict.vars = {'pAmp', 'pWidth', 'tFreq', 'tRecov'};
tDict.conds = unique([pAmp, pWidth, tFreq, tRecov], 'rows');

Nconds = size(tDict.conds, 1);

tDict.trlList = nan(Nsweeps, 1);
for a = 1:Nconds
    tmp = tDict.conds(a,:);
    l_cond = ismember([pAmp, pWidth, tFreq, tRecov], tmp, 'rows');
    tDict.trlList(l_cond) = a;
end



