function [fitparams, fval] = fit_vca_model(raw, pOnTimes, model)

% inputs:
%
%  raw: The PSC amplitudes {cell array}
%  pOnTimes: Stimulus times in seconds,  {cell array, one for each stim type}
%  model: a string indicating the number of facilitating and depressing
%         terms. Examples: 'DF', 'DDF', 'DDFF'

% assert an ERR_TYPE
ERR_TYPE = 'RMS'; % can be 'RMS' or 'LnQ'

% check input args
assert(~isempty(raw), 'ERROR: no data were supplied')
assert(~isempty(pOnTimes), 'ERROR: no pulse times were supplied')
check_vca_model(model); % will throw and error if there are problems

%
% find the best fitting STP parameters.
%
%%%%%%%%%%%%%%%%%%%%
% make a function called parse_model_coeffs that takes a vector of all
% coefficients (like the one returned from generate_guesses_and_bounds or
% from FMINUNC and returns vectors for D, F, tauD, tauF
% 
[guesses, upper_bounds, lower_bounds] = generate_vca_guesses_and_bounds(model);
N_pts_per_dim = 4;
N_cust_start_points = N_pts_per_dim .^ numel(model);
custpts = get_vca_startpoints(lower_bounds, upper_bounds, N_pts_per_dim);
N_start_points = max([N_cust_start_points, 25]); % a minimum of 100 start points;
N_start_points = min([N_start_points, 500]); % a maximum of 5000 start points

opts = optimoptions(@fmincon, 'Algorithm', 'interior-point',...  
                               'TolFun', 1e-8,...
                               'TolX',   1e-8);
problem = createOptimProblem('fmincon', 'objective', @fit_vca_err,...
                            'x0', guesses,...
                            'lb', lower_bounds,...
                            'ub', upper_bounds,...
                            'options', opts);


ms = MultiStart;
ms.UseParallel = 'always';
[fitparams, fval, exitflag, output, manymins] = run(ms, problem, N_start_points);


% fix the fitparams if the output was junk
if exitflag <= 0
    fitparams = nan(size(fitparams));
end




    function err = fit_vca_err(params)

        [d, tau_d, f, tau_f] = parse_vca_model_coeffs(params, model);

        pred = cellfun(@(x) nan(size(x)), raw, 'uniformoutput', false);
        for i_cond = 1:numel(raw)
            for i_trl = 1:size(raw{i_cond}, 3)
                A0 = raw{i_cond}(1,1,i_trl);
                pred{i_cond}(:,1,i_trl) = predict_vca_psc(pOnTimes{i_cond}, d, tau_d, f, tau_f, A0);
            end
        end

        % pool the errors across pulse train types. Ingore the first pulse
        % (becuse the error is artifactually = zero). Do some error
        % checking along the way.
        switch ERR_TYPE
            case 'LnQ'
                err = cellfun(@(x,y) abs(log(x./y)), raw, pred, 'uniformoutput', false); % log of ratios
            case 'RMS'
                err = cellfun(@(x,y) (x-y).^2,  raw, pred, 'uniformoutput', false); % squared error
        end
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
        switch ERR_TYPE
            case 'LnQ'
                err = sum(err); % for log of ratios
            case 'RMS'
                err = sqrt(mean(err)); % for RMS error
        end

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



