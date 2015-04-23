function tDict = outerleave(ax, ledIdx)
%
%  EXAMPLE    tDict = outerleave(ax, ledIdx)
% 
%  tDict.vars     => {'pAmp', 'pWidth', 'tFreq', 'tRecov'}  a reference to colums in tDict.conds
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

% iterate over the sweeps determining the pulse amplitude, width, freq
[pWidth, pAmp, tFreq] = deal(nan(Nsweeps, 1));
tRecov = zeros(Nsweeps, 1); % using a numeric value as the default b/c nans will make 'unique' bonk later in the function
for swp = 1:Nsweeps
    
    pOnTimes = ax.tt(xUp(:,swp));
    pOffTimes = ax.tt(xDown(:,swp));
    
    pWidth(swp) = mean(pOffTimes-pOnTimes) - si; % need to subtract one due to the OBO error introduced by the thresholding procedure
    
    pAmp(swp) = max(tmpWF(:,swp));
    pAmp(swp) = round(pAmp(swp).*100) ./ 100; % round to the one hundreths place
    
    if numel(pOnTimes)>1
        tFreq(swp) =  1./(pOnTimes(2)-pOnTimes(1));
        tFreq(swp) = round(tFreq(swp).*10)./10; % round to the tenth place
        
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




