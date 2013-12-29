function stro = DTfilterquesttrials(stro, nskip, WINSIZE, LOWERPERFBOUND, UPPERPERFBOUND)
%
%   EXAMPLE:  stro = DTfilterquesttrials(stro, nskip, [WINSIZE],[LOWERPERFBOUND],[UPPERPERFBOUND]);
%
% Takes a stro structure from a DTspot Quest experiment and strips out
% trials that don't meet user-defined criteria.
% Much of this code was taking from DTquestUnpack.m.
% 
% nskip: how many trials to omit from the begining of each staircase. 
% DON'T USE THIS OPTION IF YOU'RE CALCULATING PSYCHOPHYSICAL THRESHOLDS
%
% This filter acts on each spatial frequency/color direction combination
% independently.  Conditions are rejected (and stripped from the output stro
% structure) if the number of corrects over the last WINSIZE trials exceeds
% UPPERPERFBOUND or is below LOWERPERFBOUND.  These arguments can be given
% as counts or proportions.
%
% Defaults:
%       WINSIZE: 15 trials.
%       LOWERPERFBOUND: 10 correct
%       UPPERPERFBOUND: 13 correct
% (Assuming correct~bino(15,.82) these bounds cut off 12% on the low end and 22%
% on the high end.)
%
% Set WINSIZE to Nan if you don't want to filter trials based on %correct
%
% CAH 03/09 GDLH 03/09
%
% Set 'nskip' to the string 'PaperDefault' to cut all conditions that are 
% not achromatic or L-M and have a spatialperiod of [206 58 16]
%

paperdefault = 0;
if (isempty(nskip))
    nskip = 0;
end
if ~exist('WINSIZE', 'var')
    WINSIZE = 15;
end
if ~exist('PERFRANGE', 'var')
    UPPERPERFBOUND = 13;
    LOWERPERFBOUND = 10;
else
   if (UPPERPERFBOUND < 1 && UPPERPERFBOUND ~= 0)
       UPPERPERFBOUND = UPPERPERFBOUND*WINSIZE;
   end
   if (LOWERPERFBOUND < 1 && LOWERPERFBOUND ~= 0)
       LOWERPERFBOUND = LOWERPERFBOUND*WINSIZE;
   end
end

if (strcmp(nskip,'PaperDefaults'))  % Defaults for the paper
   paperdefault = 1;
   nskip = 10;
   UPPERPERFBOUND = 0;
   LOWERPERFBOUND = 0;
   WINSIZE = 0;
end


% For paperdefaults = 1
standardperiods = [206 58 16];
standardcolors = mkbasis([1 1 1; 1 -1 0]')';

%unpack the necessary experimental params.
DTindicies;
colorDirs = reshape(stro.sum.exptParams.RF_colors, 3,3)';
colorDirs(sum(abs(colorDirs), 2)==0, :) = [];
colorDirs = mkbasis(colorDirs');
spatialPeriods = sort(unique(stro.trial(:,gaborLambdaInd)));
for clr = 1:size(colorDirs, 2);
    for spfr = 1:length(spatialPeriods);
        %pull out the appropriate trials
        l_sfs = stro.trial(:,gaborLambdaInd) == spatialPeriods(spfr);
        l_clr = stro.trial(:,colorDirInd) == clr;
        correct = stro.trial(l_sfs & l_clr, correctInd);
        reject = 0;
        if (paperdefault)
            if (~ismember(spatialPeriods(spfr),standardperiods))
                reject = 1;
               % disp(['Rejecting period: ',num2str(spatialPeriods(spfr))]);
            end
            if (~any(standardcolors*colorDirs(:,clr) > .99));
                reject = 1;
              %  disp(['Rejecting color: ',num2str(colorDirs(:,clr)')]);
            end
        end
        if (sum(l_sfs & l_clr) < WINSIZE)
            reject = 1;
        elseif (~isnan(WINSIZE))
            finalperf = sum(correct(end-WINSIZE+1:end));
            if (finalperf > UPPERPERFBOUND) || (finalperf < LOWERPERFBOUND)
                reject = 1;
            end
        end
        if (reject)
            stro.ras(l_sfs & l_clr,:) = [];
            stro.trial(l_sfs & l_clr,:) = [];
            stro.sum.absTrialNum(l_sfs & l_clr) = [];
        else
            % Getting rid of trials at the begining of the run
            % (as long as whole run isn't rejected)
            if (nskip > 0)
                idxs = find(l_sfs & l_clr);
                if (length(idxs) > nskip)
                    stro.ras(idxs(1:nskip),:) = [];
                    stro.trial(idxs(1:nskip),:) = [];
                    stro.sum.absTrialNum(idxs(1:nskip)) = [];
                end
            end
        end
    end
end


