function out = NTpreprocess(stro,RTHRESH,STBLTHRESH,SORTTHRESH)

% out = NTpreprocess(stro,RTHRESH,STBLTHRESH,SORTTHRESH)
%
% A function for preprocessing a NeuroThresh stro file.
% The first output is a matrix with the following columns:
% 1) color index (a unique integer identifying the color direction)
% 2) L component of color direction
% 3) M component of color direction
% 4) S component of color direction
% 5) Threshold (point at which the staircase converged (last contrast
% tested if out of gamut).
% 6) Level (in which iteration round the color direction was probed)
% 7) Out of gamut boolean (0 = threshold measurement in gamut, 1 = no
% threshold measured)
%
% Filtering is accomplished as follows:
% 1) Any color direction with zero stepsize (probe trials) are excluded.
% 2) Any color directions that were cut off prematurely (by closing the
% file) are excluded.
% 3) STBLTHRESH: If the number of contrast troughs (peaks) in the staircase that 
% exceed the final (fall below) the contrast value exceeds STBLTHRESH, 
% the color direction is excluded. So 0 is the strictest value.
% 4) If more than SORTTHRESH of the step directions (increase or decrease in contrast)
% disagree with the directions expected on the basis of posthoc spike
% sorting, the color direction is excluded.  (set to 0 by default)
% 5) If the rank correlation coefficient between contrast and spike count
% is less than RTHRESH, the color direction is excluded.  (Out of gamut
% points are not subject to this constraint).
% GDLH 12/15/09

if (nargin < 4)
    SORTTHRESH = 0;
end
% Getting a bunch of important stuff from the stro file.
spikeidx = strcmp(stro.sum.rasterCells(1,:),getSpikenum(stro,'first'));
coloridxs = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'coloridx'));
levels = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'levelidx'));
reversals = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'reversal'));
stepsize = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stepsize'));

lmsidxs = [find(strcmp(stro.sum.trialFields(1,:),'lcont'))...
    find(strcmp(stro.sum.trialFields(1,:),'mcont'))...
    find(strcmp(stro.sum.trialFields(1,:),'scont'))];

stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_on'));
stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_off'));

lms = stro.trial(:,lmsidxs);
spikerates = zeros(size(stro.trial,1),1);
for i = 1:size(stro.trial,1)
    spiketimes = stro.ras{i,spikeidx};
    nspikes = sum(spiketimes > stimon_t(i)+stro.sum.exptParams.latency/1000 & spiketimes < stimoff_t(i));
    spikerates(i) = nspikes./(stimoff_t(i)-stimon_t(i)-stro.sum.exptParams.latency/1000);
end

uniquecoloridxs = unique(coloridxs);

fundamentals = stro.sum.exptParams.fundamentals;
if isnan(fundamentals)
    out = [];
    return
end
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = stro.sum.exptParams.mon_spd;
mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals'*mon_spd;

% OK, got the relevant information from the stro file.  Now do the actual
% crunching.
out = [];
once = false;
for i = uniquecoloridxs'
    L = logical(coloridxs == i);
    if (any(stepsize(L) == 0))  % Skip color directions with 0 stepsize
        continue;
    end
    rev = reversals(L);
    ss = stepsize(L);
    lev = unique(levels(L));
    [u,s,v] = svd(lms(L,:));
    unitvector = v(:,1);
    contrasts = lms(L,:)*unitvector;
    if(sum(contrasts) < 0)
        unitvector = -unitvector;
        contrasts = -contrasts;
    end
    contrastpeaks = contrasts(rev == 1);
    contrasttroughs = contrasts(rev == -1);
    
    tmpdata = [];  % vector into which we put data needed for filtering trials
    
    %% tmpdata(1)
    % Number of high to low reversals that fall below the final contrast
    tmpdata(1) = sum(contrastpeaks < contrasts(end));
    
    %% tmpdata(2)
    % Number of low to high reversals that fall above the final contrast
    tmpdata(2) = sum(contrasttroughs > contrasts(end));
    
    %% tmpdata(3)
    % Did we go out of gamut in this color direction?
    % Did the experiment terminate during measurement along this color direction?
    GV1 = GamutViolation(contrasts(end)*(1+ss(end)), unitvector, M, stro.sum.exptParams.bkgndrgb);
    GV2 = GamutViolation(contrasts(end)*(1-ss(end)), unitvector, M, stro.sum.exptParams.bkgndrgb);
    thresh = contrasts(end);
    if (sum(abs(rev)) < stro.sum.exptParams.nreversals)
        if (GV1)
            tmpdata(3) = 1;
            thresh = thresh*(1+ss(end));
        elseif (GV2)
            tmpdata(3) = 1;
            thresh = thresh*(1-ss(end));
        else
            tmpdata(3) = nan;  % incomplete color direction
        end
    else
        tmpdata(3) = 0;  % Stayed within gamut
    end
    
    if isempty(ss(end) > min(stepsize(stepsize > 0)))
        continue
    end
    % GDLH ugly temporary hack  2/5/11
    % The problem is that in some datafiles there are staircases with two
    % adjacent reversals in the same direction, which should never happen.
    if (ss(end) > min(stepsize(stepsize > 0)) && tmpdata(3) == 0)
        tmpdata(3) = nan;
    end
    
    %% tmpdata(4)
    % Do the reversals implemented online agree with the placement of reversals on
    % the basis of spikes sorted post-hoc?
    Loverthresh =  spikerates(L) > stro.sum.exptParams.threshold;
    if (length(contrasts) > 1)
        tmpdata(4) = sum(xor(Loverthresh(1:end-1)', sign(diff(contrasts))'== -1));
        % Contrast trajectory not in agreement with that calculated from
        % post-hoc sorted spikes
    else
        tmpdata(4) = 0;
    end
    
    %% tmpdata(5)
    % The rank correlation coefficient between firing rate and contrast
    r = corr([spikerates(L), contrasts],'type','Spearman');
    if (tmpdata(3) == 0)  % If NOOG
        tmpdata(5) = r(1,2);
    else
        tmpdata(5) = 1;
    end
    
    if (sum(stepsize(L)) == 0)
        % This should never happen
        keyboard
    end

    if (ss(end) > min(stepsize(stepsize > 0)) && tmpdata(3) == 0)
        % This should never happen
        keyboard
    end
    
    Lviolation = [tmpdata(1) > STBLTHRESH, tmpdata(2) > STBLTHRESH, isnan(tmpdata(3)), tmpdata(4)>SORTTHRESH, tmpdata(5) < RTHRESH];
    if ~any(Lviolation)
        try
         out = [out; i unitvector' thresh lev tmpdata(3)];
        catch exception
            continue
        end
    end
end
%disp(['Fraction of color directions discarded: ',num2str(1-length(out)/length(uniquecoloridxs))]);
end