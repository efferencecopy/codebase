function out = subtractEPSC(epsc, raw)

end

% NOTES
%
% From Cruz-Martin & Schweizer 2008
%
%  The scaled EPSC (at -45) can be found by taking the driving force for
%  excitation at -45 and dividing by the driving force at Erev_IPSC. This
%  is the fraction of the driving force for excitation at -45 relative to
%  the driving force at a holding potential for pure excitation. Than
%  multiply this by the current observed in the case where you only
%  measured excitation (Erev_IPSC):
%
% Iscale = (Vm - Erev_EPSC) / (Erev_IPSC -Erev_EPSC) * I_pure_epsc


% METHOD TWO
%
% Start with the combined E/I plot. Find where the inward current stops
% growing linearly (maybe by dI/dt = constant, or the second derivitive or
% something...) Then use the linear portion as the analysis domain and
% scale the pure_EPSC to best match the E/I curve over the analysis
% window...