function [troughidx, peakidx]  = anlyMod_getWFepochs(snippet, tt, condition, pWidth, photoDelay, direction)

% initialize the outputs
[troughidx, peakidx] = deal(nan);


if any(strcmpi(condition,  {'nbqx_apv_cd2', 'nbqx_apv'})) % should be a negativity followed by a positivity
        
        trough_window = (tt >= pWidth+photoDelay) & (tt <= 0.0065);
        troughval = min(snippet(trough_window)); % only look after the pulse has come on
        troughidx = find(snippet == troughval);
        assert(numel(troughidx)==1, 'ERROR: too many trough vals')
        
    elseif any(regexpi(condition, 'ttx')) % could be negative or positive depending on distance to stim site
        
        trough_window = (tt >= pWidth+photoDelay) & (tt <= 0.0065);
        if strcmpi(direction, 'outward')
            troughval = max(snippet(trough_window));
        elseif strcmpi(direction, 'inward')
            troughval = min(snippet(trough_window)); % only look after the pulse has come on
        end
        troughidx = find(snippet == troughval);
        assert(numel(troughidx)==1, 'ERROR: too many trough vals')
        
    elseif any(strcmpi(condition,   {'FV_Na', 'FV_Na_Ca2_mGluR'}))
        
        trough_window = (tt >= pWidth+photoDelay) & (tt <= 0.004);
        troughval = min(snippet(trough_window)); % only look after the pulse has come on
        troughidx = find(snippet == troughval);
        troughtime = tt(troughidx);
        assert(numel(troughidx)==1, 'ERROR: too many trough vals')

        peak_window = (tt > troughtime) & (tt <= 0.0065);
        peakval = max(snippet(peak_window)); % only look after the pulse has come on
        peakidx = find(snippet == peakval);
        assert(numel(peakidx)==1, 'ERROR: too many peak vals')
        assert(peakidx > troughidx, 'ERROR: negativity does not lead the positivity')
        
    elseif any(strcmpi(condition,  'synapticTransmission'))
        
        trough_window = (tt >= pWidth+photoDelay) & (tt <= 0.0065);
        if strcmpi(direction, 'outward') % direction of the fEPSP could depend on distance to stim site.
            troughval = max(snippet(trough_window));
        elseif strcmpi(direction, 'inward')
            troughval = min(snippet(trough_window)); % only look after the pulse has come on
        end
        troughidx = find(snippet == troughval);
        assert(numel(troughidx)==1, 'ERROR: too many trough vals')
        
end
