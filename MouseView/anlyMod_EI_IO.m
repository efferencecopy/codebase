function params = anlyMod_EI_IO(params)

% the goal is to have E/I ratios for each pulse (single or train) as a
% function of pulse intensity/freq/width (all the things that can be
% interleaved)

% approach: calculate conductance trace, and peak-by-pulse for each data
% file within a "group". Figure out which ax-files correspond to
% excitation, and to inhibition. Then match the condition types accordingly
% (so that excitation and inhibition are compared for the same stimulus).
% The resulting data could be a list of stimulus conditions and a list of
% E/I ratios associated with each condition.

% remember that each channel can have a unique holding potential for E or
% I.


% Each pulse should have its own analysis window, from 3 ms after the
% pulse, to about 15 ms after the pulse. the upper limit should not encroch
% into the next pulse epoch (if present).



% loop through each file, condition, and channel.
nFiles = numel(params.files{2});
for i_fid = 1:nFiles;
    
    nConds = size(params.tdict{i_fid}.conds, 1);
    for i_cond = 1:nConds;
        
        % find the pulse on indicies, and identify an analysis window that
        % is appropriate for this condition
        t_anlyWindowStart = 0.002; % in sec
        t_anlyWindowEnd = 0.010;
        t_anlyStartIdx = floor(t_anlyWindowStart .* params.ax{i_fid}.head.sampRate);
        t_anlyEndIdx = ceil(t_anlyWindowEnd .* params.ax{i_fid}.head.sampRate);
        
        nCh = size(params.avg.trace_pA{i_fid},2);
        for i_ch = 1:nCh
            
            % figure out what the driving force is, and make sure that the
            % Vhold is close enough to the desired value specified by the
            % "isolatedCurrents" field in the physology_notes.m
            vhold = params.avg.vhold{i_fid}{i_cond, i_ch};
            drivingForce = [];
            fieldName = [];
            match = false;
            i_rev = 1;
            while ~match && i_rev<=size(params.isolatedCurrents,1);
                Erev = params.isolatedCurrents{i_rev, 3};
                if numel(Erev)>1
                    Erev = Erev(i_ch);
                end
                
                match = abs(vhold-Erev)<1;
                if match
                    drivingForce = Erev - params.isolatedCurrents{i_rev, 4};
                    drivingForce = abs(drivingForce); % in mV
                    fieldName = params.isolatedCurrents{i_rev, 1};
                end
                
                i_rev = i_rev+1;
            end
            
            if ~match
                continue
            end
            
            % pull out the max conductance for each pulse
            pOnIdx = params.tdict{i_fid}.pOnIdx{i_cond, i_ch};
            nPulses = numel(pOnIdx);
            for i_pulse = 1:nPulses
                
                % define an analysis window for each pulse
                assert(t_anlyEndIdx < (pOnIdx(2)-pOnIdx(1)), 'ERROR: analysis window includes following pulse interval');
                anlyWindow = pOnIdx(i_pulse)+t_anlyStartIdx : pOnIdx(i_pulse)+t_anlyEndIdx;
                snippet = params.avg.trace_pA{i_fid}{i_cond, i_ch}(anlyWindow);
                
                % now determine the maximal deviation from zero (positive
                % or negative). 
                [~, maxidx] = max(abs(snippet));
                maxval_pA = snippet(maxidx);
                
                % convert to conductance
                maxval_amps = maxval_pA ./ 1e12;
                drivingForce_volts = drivingForce ./ 1e3;
                maxval_siemans = maxval_amps ./ drivingForce_volts;
                maxval_nS = maxval_siemans .* 1e9;
                
                % store the maxval conductance
                
            end
            
        end
    end
end

