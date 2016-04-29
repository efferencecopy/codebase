function [d, f, dTau, fTau] = fitTau2STP(in, psctype, channel, method)

% inputs:
%
%  in     -> a structure of corresponding to a single whole cell recording.
%            It should contain one or more 'expt' fields which has the data
%            for a specific pulse train condition
%
%  in.epxt.(pType).pOnTimes        -> the time of pulses for this pulse-train stimulus
%  in.expt.(pType).stats.EPSCamp   -> the amplitudes of the EPSC in response to eacn pulse
%
%  psctype  -> either 'EPSCamp' or 'IPSCamp'
%
%  channel  -> which headstage contains the data?


% aggregate data across pulse train types
pOnTimes = {};
raw = {};
ttypes = fieldnames(in.expt);
for i_type = 1:numel(ttypes)
    pOnTimes{i_type} = in.expt.(ttypes{i_type}).pOnTimes;
    raw{i_type} = in.expt.(ttypes{i_type}).stats.(psctype){channel};
    raw{i_type} = mean(raw{i_type},3); % mean across sweeps
end


%
% now go about finding the best fitting STP parameters.
%
%%%%%%%%%%%%%%%%%%%%
guesses = [0.6 0.8 .5 5 1.5 2]; 
upperbound = [1 1 5 15 5 15]; % [d1, d2, tau_d1, tau_d2, f1, tau_f1]
lowbound = [0 0 0 0 0 0]; % [d1, d2, tau_d1, tau_d2, f1, tau_f1]
opts = optimoptions(@fmincon, 'Algorithm', 'interior-point');
problem = createOptimProblem('fmincon', 'objective', @fittau_rms,...
                            'x0', guesses,...
                            'lb', lowbound,...
                            'ub', upperbound,...
                            'options', opts);

switch method
    case 'multistart'
        ms = MultiStart;
        ms.UseParallel = 'always';
        [fitparams, ~, exitflag] = run(ms, problem, 50);
    case 'global'
        gs = GlobalSearch;
        [fitparams, ~, exitflag] = run(gs, problem);
end


% assign the outputs
if exitflag <= 0
    fitparams = nan(size(fitparams));
end
d = fitparams(1:2);
dTau = fitparams(3:4);
f = fitparams(5);
fTau = fitparams(6);



    function err = fittau_rms(params)

        k_d = params(1:2);
        tau_d = params(3:4);
        k_f = params(5);
        tau_f = params(6);

        pred = {};
        for i_cond = 1:numel(raw)
            A0 = raw{i_cond}(1);
            pred{i_cond} = predictPSCfromTau(pOnTimes{i_cond}, k_d, tau_d, k_f, tau_f, A0);
        end

        % pool the errors across pulse train types. Ingore the first pulse
        % (becuse the error is artifactually = zero). Do some error
        % checking along the way.
        err = cellfun(@(x,y) x(2:end)./y(2:end), raw, pred, 'uniformoutput', false);
        sizeMatch = cellfun(@(x,y) numel(x)==(numel(y)-1), err, raw); % subtracting off 1 to account for the fact that I'm ignoring the first pulse
        assert(all(sizeMatch), 'ERROR: unexpected dimensions after step 1')
        
        err = cellfun(@(x) x(:), err, 'uniformoutput', false);
        
        err = cat(1, err{:});
        sizeMatch = sum(cellfun(@(y) (numel(y)-size(y,3)), raw)) == numel(err);
        assert(sizeMatch, 'ERROR: unexpected dimensions after step 2')
        
        % compute the log of the fractional error
        err = sum(abs(log10(err))); % big negative powers of ten are good

    end

end %main function







%
%
%  TESTING ROUTINES
%
%




function test_recovery() %#ok<*DEFNU>

    % this simulation just verifies that the dynamical variables (D1, D2,
    % F1) recover according to 1st order kinetics

    D1o = 0.5;
    D2o = 0.8;
    tau_d1 = 0.5;
    tau_d2 = 1;
    F1o = 1.8;
    tau_f1 = 0.8;
    
    % The simulation assumes that a presynaptic action potential occurs at
    % time zero, and then keeps track of D1, D2, and F1 following the
    % spike.
    dt = 0.001;
    tt = 0:dt:5;
    [D1, D2, F1] = deal(nan(1, numel(tt)));
    D1(1) = D1o;
    D2(1) = D2o;
    F1(1) = F1o;
    
    
    for i_t = 2:numel(tt);
       
        % let the system recover according to the time constants and the
        % asympototic values of D and F (all = 0)
        D1(i_t) = 1 - ((1-D1o) .* exp(-tt(i_t)./tau_d1));
        D2(i_t) = 1 - ((1-D2o) .* exp(-tt(i_t)./tau_d2));
        F1(i_t) = 1 + ((F1o-1) .* exp(-tt(i_t)./tau_f1));
    end

    
    figure, hold on,
    plot(tt, [D1(:), D2(:), F1(:)])
   
end



