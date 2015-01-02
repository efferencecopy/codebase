function tDict = outerleave(ax, ledIdx)
% Approach:
%
% 1) unpack the .abf file as per usual
%
% 2) find the threshold crossing for all sweeps. based off the threshold
% crossings, determing the pulse width, amplitude, IPI. Look at all of the
% sweeps to make sure there are the appropriate number of unique
% conditions, and make sure that the digitized pulse widths are consistent
% with the number of samples that should be in a specific pulse...
%
% 3) leave the sweeps interleaved but create a dictionary that can be used
% to extract the data. There needs to be a correspondence between
% dictionary entries and trial types.
% 
%  tDict.vars     => {'pAmp', 'pWidth', 'tFreq'}  a reference to colums in tDict.conds
%  tDict.conds    => matrix of values. each ROW corresponds to a unique type
%  tDict.trlList  => matrix of scalars that map each sweep onto a 'condition', one column for each channel recorded

thresh = 0.025; % volts

tmpWF = ax.dat(:,ledIdx,:);
tmpWF = permute(tmpWF, [1,3,2]);
[~, Nsweeps] = size(tmpWF);

above = tmpWF > thresh;
change = [zeros(1,Nsweeps); diff(above, 1, 1)];
xUp = change == 1;
xDown = change == -1;
si = 1./ax.head.sampRate; % the sample interval...

% iterate over the sweeps determining the pulse amplitude
[pWidth, pAmp, tFreq] = deal(nan(Nsweeps, 1));
for swp = 1:Nsweeps
    
    pOnTimes = ax.tt(xUp(:,swp));
    pOffTimes = ax.tt(xDown(:,swp));
    
    pWidth(swp) = mean(pOffTimes-pOnTimes) - si; % need to subtract one due to the o.b.o.e introduced by the thresholding procedure
    
    pAmp(swp) = max(tmpWF(:,swp));
    pAmp(swp) = round(pAmp(swp).*100) ./ 100; % round to the one hundreths place
    
    if numel(pOnTimes)>1
        tFreq(swp) =  1./(pOnTimes(2)-pOnTimes(1));
        tFreq(swp) = round(tFreq(swp).*10)./10; % round to the tenth place
    else
        tFreq(swp) = 0; % a single pulse
    end
    
end


% now determine the number of unique trial types
tDict.vars = {'pAmp', 'pWidth', 'tFreq'};
tDict.conds = unique([pAmp, pWidth, tFreq], 'rows');
Nconds = size(tDict.conds, 1);

tDict.trlList = nan(Nsweeps, 1);
for a = 1:Nconds
    tmp = tDict.conds(a,:);
    l_cond = ismember([pAmp, pWidth, tFreq], tmp, 'rows');
    tDict.trlList(l_cond) = a;
end




