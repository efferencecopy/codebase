function [trialSnips, pulseSnips] = getTimeLockedLFP(lfp_data, blk, STIMTYPE, NSX)


stim_ch_idx = blk.sum.idx.(['stim_on_', NSX]);
tmpWFs = blk.ras(:, stim_ch_idx);
tmpWFs = cellfun(@(x) x./1000, tmpWFs, 'uniformoutput', false); % outerleave expects things to be in volts.
sampFreq_nsx = blk.sum.(NSX).MetaTags.SamplingFreq;
tdict = outerleave(tmpWFs, sampFreq_nsx);

switch STIMTYPE
    case 'sinusoid'
        tdict.conds = tdict.conds(1,:);
        tdict.trlList(:) = 1;
        preTime = 0.050;
        postTime = 0.010;
        
    case 'train'
        preTime = 0.050;
        postTime = 0.250;
end


% get the pulse onset/offset times
Ntrials = size(blk.ras, 1);
thresh = min(tdict.conds(:,1).*1000.*0.9); % in mV;
[pulseOffIdx, pulseOnIdx] = deal({});
for i_trl = 1:Ntrials
    
   trl_stim = blk.ras{i_trl, stim_ch_idx}; 
   aboveThresh = trl_stim > thresh;
   on_idx = [nan, diff(aboveThresh)] == 1;
   on_idx = find(on_idx);
   off_idx = [nan, diff(aboveThresh)] == -1;
   off_idx = find(off_idx);
   
   assert(numel(off_idx) == numel(on_idx), 'ERROR: mismatch in pulse number')
   
   pulseOnIdx{i_trl} = on_idx;
   pulseOffIdx{i_trl} = off_idx;
end

%
% pull out snippets of data from before/after each pulse in the pulse
% train. baseline subtract each of them from their pre-pulse baseline.
% Average across trials
pulseSnips = {};
trialSnips = {};
Nptypes = size(tdict.conds, 1);
for i_ptype = 1:Nptypes
    
    % give this pulse type a name
    p_amp = tdict.conds(i_ptype, 1);
    p_width = tdict.conds(i_ptype, 2).*1e6;
    p_tf = tdict.conds(i_ptype, 3);
    fldname = ['amp', num2str(p_amp), '_pw', num2str(p_width), '_tf', num2str(p_tf)];
    
    % basic stuff for this ptype
    trial_idx = tdict.trlList == i_ptype;
    trial_nums = find(trial_idx);
    Nchannels = size(lfp_data,2);
    Npulses = unique(cellfun(@numel, pulseOnIdx));
    assert(numel(Npulses) == 1, 'ERROR: mismatch in pulse numbers')
    
    % initialize the outputs
    trialSnips.(fldname) = repmat({[]}, 1, Nchannels);
    
    
    for i_trl = 1:numel(trial_nums)
        
        trlNum = trial_nums(i_trl);
        
        for i_ch = 1:Nchannels
            
            % set aside the trial snippets
            idx_on_first = pulseOnIdx{trlNum}(1);
            preTimeSamps = round((preTime .* sampFreq_nsx));
            idx_first = idx_on_first - preTimeSamps;
            ipi = 1/p_tf;
            idx_last = idx_on_first + round(((ipi .* Npulses)+postTime) .* sampFreq_nsx);
            
            snip = lfp_data{trlNum, i_ch}(idx_first : idx_last);
            snip = snip - mean(snip(1:preTimeSamps));
            trialSnips.(fldname){1, i_ch} = cat(1, trialSnips.(fldname){i_ch}, snip);
            trialSnips.preTime = preTime;
            trialSnips.postTime = postTime;
            
            
        end
        
    end
    
end
