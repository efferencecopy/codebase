function ci = boot2ci(monk, cell, expt, nStraps, confInterval)

% This funtion will accept neurometric and psychometric data and create
% confidence intervals for the threshold ratio(s) via nonparametric
% bootstrap. The first three input arguments are created by DTmocsUnpack.m.
% "confInterval" should be a percentage between 0 and 1;
%
% CAH 11/11/09


%create distributions of thresholds for all color/sfs combos. Do this for
%both the psychometric and neurometric functions.
nStraps = ceil(sqrt(nStraps));%b/c I'm going to use all pairwise comparisons;
pairs = fullfact([nStraps, nStraps]); 
alphaLevel = 1-confInterval; 
nSfs = length(expt.sfs);
nColors = size(expt.standColors,1);
for a = 1:nColors;
    for b = 1:nSfs;
        if isnan(cell.alpha(a,b))  %bail if the neurometric function is flat
            ci.up(a,b) = nan;
            ci.lo(a,b) = nan;
        else
            psyErrs = btstrapPsyFun(monk.alpha(a,b), monk.beta(a,b), monk.gamma(a,b),...
                expt.norms{a,b}, monk.nTrials{a,b}, nStraps, 100);
            neuroErrs = btstrapPsyFun(cell.alpha(a,b), cell.beta(a,b), cell.gamma(a,b),...
                expt.norms{a,b}, cell.nTrialsIn{a,b}, nStraps, 100);
            btTRs = neuroErrs.alpha(pairs(:,1)) ./ psyErrs.alpha(pairs(:,2)); %all pairwise comparisons
            
            %filter out the neurons that give rise to super high btstrp
            %estimates (i.e. 10^308).... which cause inf TRs.
            if any(isinf(btTRs))
                ci.up(a,b) = nan;
                ci.lo(a,b) = nan;
            else
                %determine the bias
                TRest = cell.alpha(a,b) ./ monk.alpha(a,b);
                totalStraps = nStraps^2;
                bias = sum(btTRs < TRest)./totalStraps; %proportion less than the estimate...
                bias = norminv(bias, 0, 1); %inv of normcdf evaluated at proportion < estimate
                
                %determine the acceleration
                accel = 0; %I don't know how to estimate this, so for now I'm setting it to zero.
                a_low = getAlphaLev(bias, accel, alphaLevel./2) .* 100;  %put in % b/w 1 and 100 for prctile.m
                a_high = getAlphaLev(bias, accel, 1-(alphaLevel./2)) .* 100;
                ci.up(a,b) = prctile(btTRs, a_high);
                ci.lo(a,b) = prctile(btTRs, a_low);
            end
        end
    end %sfs
end %colors
end% function


function BCalpha = getAlphaLev(bias, accel, alphaLevel)
    z_alpha = norminv(alphaLevel, 0, 1);
    phi = bias + ((bias + z_alpha) ./ (1 - accel.*(bias + z_alpha)));
    BCalpha = normcdf(phi, 0, 1);
end

