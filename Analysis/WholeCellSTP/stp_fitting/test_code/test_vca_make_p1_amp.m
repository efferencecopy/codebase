function [p1_amps, noiseless_p1_amps] = test_vca_make_p1_amp(n_sweeps, A0_init, A0_sigma, nsparams)

% default p1 amps are just A0_init
p1_amps = ones(n_sweeps, 1) .* A0_init;
noiseless_p1_amps = p1_amps; 

% add noise if needed
if exist('A0_sigma', 'var') && ~isempty(A0_sigma);
    p1_amps = p1_amps + normrnd(0, A0_sigma, [n_sweeps, 1]);
    p1_amps = max(p1_amps, zeros(size(p1_amps))); % don't let the EPSC be negative;
end

% add a nonstationarity if needed
if exist('nsparams', 'var') && ~isempty(nsparams);
    switch nsparams.type
        case 'linear'
            line = linspace(1, 1-nsparams.prcnt_rundown, n_sweeps);
            noiseless_p1_amps .* p1_amps;
            p1_amps = p1_amps .* line(:);
        case 'none'
            % do nothing
    end
end
