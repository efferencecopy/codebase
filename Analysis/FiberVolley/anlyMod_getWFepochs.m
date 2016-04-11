function [troughidx, peakidx]  = anlyMod_getWFepochs(snippet, tt, condition, pWidth, photoDelay, direction, opsin)

% initialize the outputs
[troughidx, peakidx] = deal(nan);


if any(strcmpi(condition,  {'nbqx_apv_cd2', 'nbqx_apv'})) % should be a negativity followed by a positivity
    
    trough_window = (tt >= pWidth+photoDelay) & (tt <= 0.0065);
    troughval = min(snippet(trough_window)); % only look after the pulse has come on
    troughidx = (snippet == troughval) & trough_window; % make sure matches are in the window.
    troughidx = find(troughidx, 1, 'first');
    
elseif any(regexpi(condition, 'ttx')) && ~strcmpi(opsin, 'estim'); % could be negative or positive depending on distance to stim site
    
    trough_window = (tt >= pWidth+photoDelay) & (tt <= 0.0065);
    
    if strcmpi(direction, 'outward')
        troughval = max(snippet(trough_window));
    elseif strcmpi(direction, 'inward')
        troughval = min(snippet(trough_window)); % only look after the pulse has come on
    end
    troughidx = (snippet == troughval) & trough_window; % make sure matches are in the window.
    troughidx = find(troughidx, 1, 'first');
    
    
elseif any(regexpi(condition, 'ttx')) && strcmpi(opsin, 'estim'); % special case for estim
    
    trough_window = (tt >= 0) & (tt <= pWidth+photoDelay);
    troughval = min(snippet(trough_window)); % only look after the pulse has come on
    troughidx = (snippet == troughval) & trough_window; % make sure matches are in the window.
    troughidx = find(troughidx, 1, 'first');
    
    peakval = max(snippet(trough_window)); % only look after the pulse has come on
    peakidx = (snippet == peakval) & trough_window; % make sure matches are in the window.
    peakidx = find(peakidx, 1, 'first');
    
elseif any(strcmpi(condition,   {'FV_Na', 'FV_Na_Ca2_mGluR'}))
    
    trough_window = (tt >= pWidth+photoDelay) & (tt <= 0.004);
    troughval = min(snippet(trough_window)); % only look after the pulse has come on
    troughidx = (snippet == troughval) & trough_window; % make sure matches are in the window.
    troughidx = find(troughidx, 1, 'first');
    troughtime = tt(troughidx);
    
    peak_window = (tt > troughtime) & (tt <= 0.0065);
    peakval = max(snippet(peak_window)); % only look after the pulse has come on
    peakidx = (snippet == peakval) & peak_window; % make sure matches are in the window.
    peakidx = find(peakidx, 1, 'first');
    assert(peakidx > troughidx, 'ERROR: negativity does not lead the positivity')
    
elseif any(strcmpi(condition,  'synapticTransmission'))
    
    trough_window = (tt >= pWidth+photoDelay) & (tt <= 0.0065);
    if strcmpi(direction, 'outward') % direction of the fEPSP could depend on distance to stim site.
        troughval = max(snippet(trough_window));
    elseif strcmpi(direction, 'inward')
        troughval = min(snippet(trough_window)); % only look after the pulse has come on
    end
    troughidx = (snippet == troughval) & trough_window; % make sure matches are in the window.
    troughidx = find(troughidx, 1, 'first');
    
end
