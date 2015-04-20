%function [d, f, dTau, fTau] = fitTau2STP(group)

% game plan
%
%
%  use euler's method to determine the predicted response for a particular
%  set of depression and facilitation constants.
%
%  Technically, I only need to calculate the predictions once for each
%  unique spike train condition, and then calculate the residules for each
%  measured condition. 
%
%  I don't know if there's a faster way to fit a dynamical system. If there
%  is, that would be cool!
%
%  Varela 1997 used 2 depression factors and 1 facilitation factor. Each
%  factor consisted of a constant that was applied to each action
%  potential, and a recovery time constant that forced the factor back to
%  its resting state (a value of 1). 


% Step 1: import a file that's been pre-processed using the
% 'anlyMod_avgOuterleave' module. Make sure that it has the correct field.
% For example 'recovery' group is permissible. Bottom line, you're looking
% for a tdict structure.

% Step 2: iterate over channels. It would be nice if I can figure out a
% clever way of excluding junk channels.

% Step 2.1: For each trial type in the tdict, determine what the predicted
% pulse amplitude would be given a set of parameters. Do fminsearch to find
% the best fitting paramters.



i_fid = 1; % just looking at the first file for now
i_ch = 2; % just looking at the second channel for now (HS2)


% pull out the average current records and the pOn_times
raw_pa = group.avg.trace_pA{i_fid}(:,i_ch);

% preallocate the output based on the number of pulses presented:
psc_amp_pa = cellfun(@(x) x.* nan, group.tdict{i_fid}.pOnIdx(:,i_ch), 'uniformoutput', false);


% calculate the psc amplitude following each LED pulse
sampRate = group.head{i_fid}.sampRate;
preTime_idx = ceil(0.003 .* sampRate);
postTime_idx = ceil(0.010 .* sampRate);
artifactTimeout_idx = ceil(0.0015 .* sampRate);

Ntraintypes = size(group.tdict{i_fid}.pOnIdx, 1);
for i_train = 1:Ntraintypes
    
    Npulses = numel(group.tdict{i_fid}.pOnIdx{i_train});
    for i_pulse = 1:Npulses;
        pOn_idx = group.tdict{i_fid}.pOnIdx{i_train}(i_pulse);
        psc_idx = (pOn_idx+artifactTimeout_idx) : (pOn_idx+postTime_idx);
        baseline_idx = (pOn_idx-preTime_idx) : (pOn_idx+1);
        
        % the raw data
        psc_raw = raw_pa{i_train}(psc_idx);
        baseline_raw = raw_pa{i_train}(baseline_idx);
        
        % baseline subtract things
        psc_raw = psc_raw - mean(baseline_raw);
        
        % calculate the psc_amp
        [min_val, min_idx] = min(psc_raw);
        psc_amp_pa{i_train}(i_pulse) = abs(min_val);
        
    end
end










