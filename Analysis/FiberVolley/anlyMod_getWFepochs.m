function [troughidx, peakidx]  = anlyMod_getWFepochs(snippet, tt, condition, pWidth, photoDelay)

% initialize the outputs
[troughidx, peakidx] = deal(nan);

switch condition
    case {'nbqx_apv_cd2_ttx', 'nbqx_apv_cd2', 'nbqx_apv'}
        
        trough_window = (tt >= pWidth+photoDelay) & (tt <= 0.0065);
        troughval = min(snippet(trough_window)); % only look after the pulse has come on
        troughidx = find(snippet == troughval);
        if numel(troughidx) > 1
            troughidx = troughidx(1);
            warning('found too many trough indicies')
            keyboard
        end
            
        %assert(numel(troughidx)==1, 'ERROR: too many trough vals')
        
    case  'FV_Na'
        
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
        
    case 'synapticTransmission'
        
        trough_window = (tt >= pWidth+photoDelay) & (tt <= 0.0065);
        troughval = min(snippet(trough_window)); % only look after the pulse has come on
        troughidx = find(snippet == troughval);
        assert(numel(troughidx)==1, 'ERROR: too many trough vals')
        
end
