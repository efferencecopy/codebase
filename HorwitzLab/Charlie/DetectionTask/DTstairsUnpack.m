function [norms, thresh, cycperdeg, units] = DTstairsUnpack(stro)

if(stro.sum.exptParams.expt_meth ~= 2)
    error('this experiment does not use staircasing')
end

%***** CONSTANTS *********%
NTRIALS = 30;

% create the relavant indicies and recompute the relavant monitor
% calibration data
DTindicies;
bkgndrgb = [stro.sum.exptParams.bkgnd_r, stro.sum.exptParams.bkgnd_g, stro.sum.exptParams.bkgnd_b];
M = reshape(stro.sum.exptParams.m_mtx, 3, 3);
bkgndLMS = M * bkgndrgb';
x = 0:255; %the normal range of the gamma look up table
xx = linspace(0, 255, 2^16); %the desired quantization of the gammaTable
g1 = reshape(stro.sum.exptParams.gamma_table, 256, 3);
gammaTable = [spline(x, g1(:,1), xx)', spline(x, g1(:,2), xx)', spline(x, g1(:,3), xx)'];


% determine some relavent experimental characteristics
numTrials = size(stro.trial,1);
sptPeriods = unique(stro.trial(:,gaborLambdaInd));
nSptFreqs = size(sptPeriods, 1);
nColors = max(stro.trial(:,colorDirInd));


%convert from DAC gun space to cone contrasts. Do this seperately for the
%predicted and actual gun values.
predGunVals = [stro.trial(:,flashRInd), stro.trial(:,flashGInd), stro.trial(:,flashBInd)] + 1;
pred_rgb = [gammaTable(predGunVals(:,1), 1), gammaTable(predGunVals(:,2), 2), gammaTable(predGunVals(:,3), 3)];
pred_rgb = pred_rgb - repmat(bkgndrgb, numTrials, 1);
pred_lms = [(M * pred_rgb') ./ repmat(bkgndLMS, 1, numTrials)]';


actGunVals = [stro.trial(:,actflashRInd), stro.trial(:,actflashGInd), stro.trial(:,actflashBInd)];
act_rgb = [gammaTable(actGunVals(:,1), 1), gammaTable(actGunVals(:,2), 2), gammaTable(actGunVals(:,3), 3)];
act_rgb = act_rgb - repmat(bkgndrgb, numTrials, 1);
act_lms = [(M * act_rgb') ./ repmat(bkgndLMS, 1, numTrials)]';


for j = 1:nColors;
    for a = 1:nSptFreqs;
        %determine which trials have a particular color and sptFreq
        colorInd = stro.trial(:, colorDirInd) == j;
        sptFreqInd = stro.trial(:, gaborLambdaInd) == sptPeriods(a);
        tind = colorInd & sptFreqInd;
        
        if sum(tind); %if there are any of this trial type
            pred_tLMS = pred_lms(tind, :);
            pred_norms = sqrt(sum((pred_tLMS.^2), 2));
            pred_thresh = mean(pred_norms((end-(NTRIALS-1)):end));
            pred_unit = pred_tLMS ./ repmat(pred_norms, 1, 3);
            
            act_tLMS = act_lms(tind, :);
            act_norms = sqrt(sum((act_tLMS.^2), 2));
            act_thresh = mean(act_norms((end-(NTRIALS-1)):end));
            act_unit = act_tLMS ./ repmat(act_norms, 1, 3);
            
            norms{j,a,1} = act_norms;
            thresh(j,a,1) = act_thresh;
            units{j,a,1} = act_unit;
            
            norms{j,a,2} = pred_norms;
            thresh(j,a,2) = pred_thresh;
            units{j,a,2} = pred_unit;
            
            corrects = stro.trial(tind, correctInd);
            perf(j,a) = sum(corrects((end-(NTRIALS-1)):end)) ./ NTRIALS;
        end
    end
end


%what was the range of sptFreqs?
cycperdeg = 1./(sptPeriods ./ stro.sum.exptParams.pixperdeg);

%loop through the conditions and plot a figure of the trial by trial CC's
if (isunix)
    [up, down, left, right, esc] = deal([82], [81], [80], [79], [41]);
elseif (ispc)
    [up, down, left, right, esc] = deal([38], [40], [37], [39], [27]);
end

sptFreq = 1;
color = 1;
linecolors = ['k', 'r', 'b'];
h = figure;
while(1)
    [keyIsDown, secs, keyCode] = KbCheck();
    while ~keyIsDown
        [keyIsDown, secs, keyCode] = KbCheck();
        drawnow
    end
    
    if(keyCode(right) && (sptFreq<nSptFreqs))
        sptFreq = sptFreq+1;
    end
    if(keyCode(left) && (sptFreq>1))
        sptFreq = sptFreq - 1;
    end
    if(keyCode(up) && (color<nColors))
        color = color+1;
    end
    if(keyCode(down) && (color>1))
        color = color-1;
    end
    if(keyCode(esc))
        close(h);
        break
    end
    
    figure(h)
    clf
    hold on,
    plot(norms{color, sptFreq}, [linecolors(color) '.-'], 'LineWidth', 2);
    plot([length(norms{color, sptFreq})-NTRIALS+1, length(norms{color, sptFreq})], [thresh(color, sptFreq), thresh(color, sptFreq)], 'g', 'LineWidth', 2);
    hold off
    xlabel('trial number')
    ylabel('Cone Contrast')
    title(sprintf('Spt Freq No. %d of %d: %.2f cyc/deg \n Color No: %d, Perf: %.3f', sptFreq, nSptFreqs, cycperdeg(sptFreq), color, perf(color, sptFreq)));
    pause(0.2)
end



%% CONCERNS
%1) WHEN I TAKE THE DIFF(NORMS) I EXPECT TO SEE TWO VALUES (DECVAL, AND
%INCVAL) THIS IS NOT THE CASE. THERE ARE A GREAT MANY VALUES. SOME MAY
%CORRESPOND TO OUT OF RANGE CORRECTIONS, BUT I DON'T UNDERSTAND WHY THERE
%ARE SO MANY DIFFERENT VALUES!!!


% should I change diag to sum when computing norms?