function [d, f, dTau, fTau] = fitTau2STP(group)

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
        pOn_idx = group.tdict{i_fid}.pOnIdx{i_train, i_ch}(i_pulse);
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


% a quick figure as a sanity check
tmp = [psc_amp_pa{:}]';
tmp = bsxfun(@rdivide, tmp, tmp(:,1));
[x,y] = meshgrid(1:size(tmp,2),1:size(tmp,1));
figure
surf(x,y,tmp,'facealpha', 0.8, 'linewidth', 1)
set(gca, 'zscale', 'log')



%
% now go about finding the best fitting STP parameters.
%
%%%%%%%%%%%%%%%%%%%%

% initialize the dynamical variables and Tau's
nD = 2; % the number of depression factors
nF = 1; % the number of facilitation factors
[d, dTau] = deal(nan(1,nD));
[f, fTau] = deal(nan(1,nF));


% pull out the actual pulse times
pOn_idx = [group.tdict{i_fid}.pOnIdx{:, i_ch}]';
pOn_idx = bsxfun(@minus, pOn_idx, pOn_idx(:,1)); % relative to first pulse
pOn_time = pOn_idx ./ sampRate;

% make the raw data a matrix as well (instead of a cell array...)
psc_amp_pa = [psc_amp_pa{:}]';

% do the fminsearch
guesses = [0.5 0.8 0.030 0.075 1.05 0.250]; % [d1, d2, ?d1, ?d2, f1, ?f1]

[out, success, fval] = fminsearch(@(x) fittau_rms(x, psc_amp_pa, pOn_time) , guesses);








end %main function


function rms = fittau_rms(params, rawdata, ptimes)
    
    d1 = params(1);
    d2 = params(2);
    tau_d1 = params(3);
    tau_d2 = params(4);
    f1 = params(5);
    tau_f1 = params(6);
    
    % use forward euler's technique to numerically solve the dynamical
    % system and determine the predicted data give the spike times and
    % params
    
    pred = nan(size(rawdata));
    A0 = rawdata(:,1);
    pred(:,1) = A0;
    [D1, D2, F1] = deal([1]); % all dynamical variables start at one.
    for i_pulse = 2:size(ptimes,2);
        % update Ds and Fs
        ipi = ptimes(:,i_pulse) - ptimes(:,i_pulse-1);
        
        % let the system recover according to the time constants and the
        % asympototic values of D and F (all = 0)
        D1 = D1 + (1-D1) .* exp(-ipi./tau_d1);
        D2 = D2 + (1-D2) .* exp(-ipi./tau_d2);
        F1 = F1 - (F1-1) .* exp(-ipi./tau_f1);
        
        % now add the per-spike plasticity (d, f).
        D1 = D1 .* d1;
        D2 = D2 .* d2;
        F1 = F1 + f1
        
        % make a prediction
        pred(:,i_pulse) = A0 .* D1 .* D2 .* F1;
        
    end
    
    
end













