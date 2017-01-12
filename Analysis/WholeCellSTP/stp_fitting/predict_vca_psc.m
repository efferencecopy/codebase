function pred = predict_vca_psc(pOnTimes, d, tau_d, f, tau_f, A0)
%PREDICT EPSC FROM MODEL OF SYNAPTIC PLASTICITY
% 
%     pred = predictPSCfromTau(pOnTimes, d, tau_d, f, tau_f, A0)
% 
%     Predict EPSC amplitudes from dynamical model of short term synaptic
%     plasticity. The model assumes that the EPSC amplitude is equal to:
% 
%     EPSCn = EPSCo * D1 * D2 ... * Dn * F1 * F2 ... * Fn
% 
%     where D terms are depressive (betwen zero and 1) and the F terms are
%     facilitating (between 1 and Inf). The D and F terms relax back to 1 and
%     Zero, respectively, according to the associated Tau (either tau_d or
%     tau_f)
% 
%     Inputs: d, f, tau_d, and tau_f can be scalars (only one D/F term) or
%     vectors (multiple D or F terms). The dimensionality of d and tau_d must
%     be the same. The dimensionality of f and tau_f must be the same. But the
%     dimensionality of d and f can be different (i.e., there can be different
%     numbers of Depressive and Facilitating terms.
% 
%     Inputs: pOnTimes should be a vector of spike times in the presynaptic
%     fibers (or times of laser pulses for optical stimulation).
% 
%   Charlie Hass 2016

pred = nan(numel(pOnTimes),1); % by convention, pulse numbers go DOWN columns, so pred needs to be a column vector

% solve for the state variables after the first spike
pred(1) = A0;
D = 1 .* d;
F = 1 + f;
ipi = [0, diff(pOnTimes)];

for i_pulse = 2:numel(pOnTimes);
    
    % let the system recover according to the time constants and the
    % asympototic values of D and F (all = 1)
    D = 1 - ((1-D) .* exp(-ipi(i_pulse)./tau_d));
    F = 1 + ((F-1) .* exp(-ipi(i_pulse)./tau_f));
    
    assert(all(D>=0) && all(D<=1), 'ERROR: depressive terms oob')
    assert(all(F>=1), 'ERROR: facilitating terms oob')
    
    % make a prediction regarding the current spike
    pred(i_pulse) = prod(cat(1, A0(:), D(:), F(:))); % equivalent to (A0 * D1 * D2 ... * Dn * F1 * F2 ... * Fn)
    
    % now add the per-spike plasticity (d, f).
    D = D .* d;
    F = F + f;
end


