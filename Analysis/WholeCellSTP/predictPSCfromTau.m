function pred = predictPSCfromTau(pOnTimes, d, dTau, f, fTau, A0)

% define the state variables
d1 = d(1);
d2 = d(2);
tau_d1 = dTau(1);
tau_d2 = dTau(2);
f1 = f(1);
tau_f1 = fTau(1);

pred = nan(numel(pOnTimes),1); % by convention, pulse numbers go DOWN columns, so pred needs to be a column vector

% solve for the variables after the first spike
pred(1) = A0;
D1 = 1 .* d1;
D2 = 1 .* d2;
F1 = 1 + f1;
ipi = [0, diff(pOnTimes)];

for i_pulse = 2:numel(pOnTimes);
    
    % let the system recover according to the time constants and the
    % asympototic values of D and F (all = 1)
    D1 = 1 - ((1-D1) .* exp(-ipi(i_pulse)./tau_d1));
    D2 = 1 - ((1-D2) .* exp(-ipi(i_pulse)./tau_d2));
    F1 = 1 + ((F1-1) .* exp(-ipi(i_pulse)./tau_f1));
    
    assert(D1>=0 && D1<=1)
    assert(D2>=0 && D2<=1)
    assert(F1>=1)
    
    % make a prediction regarding the current spike
    pred(i_pulse) = A0 .* D1 .* D2 .* F1;
    
    % now add the per-spike plasticity (d, f).
    D1 = D1 .* d1;
    D2 = D2 .* d2;
    F1 = F1 + f1;
end