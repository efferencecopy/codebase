function [d, f, dTau, fTau] = fitTau2STP(raw, pOnTimes, ~, method)

% inputs:
%
%  raw: The PSC amplitudes {cell array}
%  pOnTimes: Stimulus times in seconds,  {cell array, one for each stim type}
%  method:  'MultiStart', or 'GlobalSearch'

if isempty(raw) || isempty(pOnTimes)
    error('No data were provided')
end

%
% find the best fitting STP parameters.
%
%%%%%%%%%%%%%%%%%%%%
guesses = [0.6 0.8 .5 5 1.5 2]; 
upperbound = [1 1 3 3 5 3]; % [d1, d2, tau_d1, tau_d2, f1, tau_f1]
lowbound = [0 0 0 0 0 0]; % [d1, d2, tau_d1, tau_d2, f1, tau_f1]
opts = optimoptions(@fmincon, 'Algorithm', 'interior-point',...  
                               'TolFun', 1e-8,...
                               'TolCon', 1e-8); 
problem = createOptimProblem('fmincon', 'objective', @fittau_rms,...
                            'x0', guesses,...
                            'lb', lowbound,...
                            'ub', upperbound,...
                            'options', opts);

switch method
    case 'multistart'
        ms = MultiStart;
        ms.UseParallel = 'always';
        ms.StartPointsToRun = 'all';
        ms.TolFun = 0;
        [fitparams, ~, exitflag] = run(ms, problem, 100);
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

        pred = cellfun(@(x) nan(size(x)), raw, 'uniformoutput', false);
        for i_cond = 1:numel(raw)
            for i_trl = 1:size(raw{i_cond}, 3)
                A0 = raw{i_cond}(1,1,i_trl);
                pred{i_cond}(:,1,i_trl) = predictPSCfromTau(pOnTimes{i_cond}, k_d, tau_d, k_f, tau_f, A0);
            end
        end

        % pool the errors across pulse train types. Ingore the first pulse
        % (becuse the error is artifactually = zero). Do some error
        % checking along the way.
        err = cellfun(@(x,y) abs(log(x./y)), raw, pred, 'uniformoutput', false); % log of ratios
        %err = cellfun(@(x,y) (x-y).^2,  raw, pred, 'uniformoutput', false); % SSE
        sizeMatch = cellfun(@(x,y) numel(x)==numel(y), err, raw); % subtracting off 1 to account for the fact that I'm ignoring the first pulse
        assert(all(sizeMatch), 'ERROR: unexpected dimensions after step 1')
        
        
        % compute the log of the fractional error. When the data match the
        % prediction, the ratio will be 1 and the log(ratio) will be zero.
        % This means we want the sum of the log(ratios) to be as close as
        % possible to zero. Take the ABS(log(ratio)) to force the code to
        % minimize the difference from zero b/c deviations on the low and
        % high side are both bad.
        err = cellfun(@(x) x(:), err, 'uniformoutput', false);
        err = cat(1, err{:});
        err = sum(err); 

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

    D1o = 0.5; % i.e., 1 * d1
    D2o = 0.8; % i.e., 1 * d2
    tau_d1 = 0.5;
    tau_d2 = 1;
    F1o = 1.8;  % i.e., 1 + f1
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



