function [fpar, fval] = fitDetectionSurface(colors, alphas, initparams, npsoiters)


    % unpack the color directions used. turn them into unit vectors.
    if iscell(colors)
        colors = cat(1, colors{:});
    end
    colordirs = bsxfun(@rdivide, colors, sqrt(sum(colors.^2, 2))); % make sure they're unit vecs
    colordirs = [colordirs; -colordirs];


    % unpack the threshold estimates.
    if iscell(alphas)
        alphas = cat(1, alphas{:}); %in CC b/w 0 and 100
    end
    alphas = [alphas; alphas];

    % find some initial guesses. these could be user supplied, or could be
    % estimated using a variety of methods.
    if isvector(initparams) && isnumeric(initparams) %actual initial guesses have been supplied
        if numel(initparams)~=10
            error('incorrect number of inital guesses')
        end
        
    elseif strcmpi(initparams, 'ellipsoid')
        
        cordinates = bsxfun(@times, colordirs, alphas);
        
        D = [cordinates.^2,...
            2.*cordinates(:,1) .* cordinates(:,2),...
            2.*cordinates(:,1) .* cordinates(:,3),...
            2.*cordinates(:,2) .* cordinates(:,3)];
        lssoln = (D' * D) \(D' * ones(size(cordinates,1),1));
        
        
        A = [lssoln(1) lssoln(4) lssoln(5);...
            lssoln(4) lssoln(2) lssoln(6);...
            lssoln(5) lssoln(6) lssoln(3)];
        [evecs, evals] = eig(A);
        betaGuess = 2;
        initparams = [betaGuess; reshape(evecs*sqrt(evals),9,1)];
        
    elseif strcmpi(initparams, 'pso')
        
        problem = setupPSOoptions();        
        problem.fitnessfcn = @coleFitErr_local;
        
        fpars_pso = nan(npsoiters, 10);
        fvals_pso = nan(npsoiters,1);
        for a = 1:npsoiters
            %[fpars_pso(a,:), fvals_pso(a)] = pso(problem);
            options = optimoptions('pso', 'TolFun', 1e-10, 'MaxIter', 10e5, 'SelfRecognitionCoefficient', 0.5, 'SocialRecognitionCoefficient', 0.3, 'NumParticles', 200);
            [fpars_pso(a,:), fvals_pso(a)] = pso(@coleFitErr_local, 10, [1; repmat(-100, 9, 1)], [3; repmat(100, 9, 1)], options);
            fvals_pso(a)
        end
        
        [~, idx_pso] = min(fvals_pso);
        initparams = fpars_pso(idx_pso,:);
    end
    
    
    % now run the fminsearch using the inital guess supplied (or estimated)
    options = optimset('MaxFunEvals',50000,'MaxIter',50000,'TolFun',10^-6,'Display','none');
    [fpar, fval] = fminsearch(@(x) coleFitErr_local(x), initparams, options);
    
    % check the output of the fminsearch against the best value of the pso.
    % The fminsearch should be better.
    if strcmpi(initparams, 'pso')
        if fvals_pso(idx_pso) < fval
            fval = fvals_pso(idx_pso);
            fpar = fpars_pso(idx_pso,:);
        end
    end
    
    
    %
    % END OF MAIN CODE (SUBFUNCTIONS BELOW)
    %
        
    
    
    % this has to be nested b/c i don't know how to use non-nested fxn with
    % pso.m
    function err = coleFitErr_local(params)
        
        mechanisms = reshape(params(2:end),3,3);
        beta = params(1);
        wholesum = sum(abs(colordirs * mechanisms).^beta, 2); % The "response" of the system assuming a particular set of mechanisms
        
        predr = (1./wholesum).^(1/beta);
        residual = log(alphas)-log(predr);
        err = sum(residual.^2);
        
    end
end




function problem = setupPSOoptions
    % set things up for the partical swarm optimization
    problem.Aineq = [];
    problem.bineq = [];
    problem.Aeq = [];
    problem.beq = [];
    problem.LB = [1; repmat(-100, 9, 1)];
    problem.UB = [3; repmat(100, 9, 1)];
    problem.nonlcon = [];
    problem.nvars = 10;
    problem.HybridFcn = {@fminsearch};
    problem.options.TolCon = 1e-15; % it seems like these make more of a difference than the others..
    problem.options.TolFun = 1e-15; % it seems like these make more of a difference than the others..
    problem.options.CognitiveAttraction = 0.5; % it seems like these make more of a difference than the others.. typical default value is 0.5
    problem.options.SocialAttraction = 0.3; %default is 1.25
    problem.options.Generations = 5e3;
    problem.options.PopulationSize = 150;
    problem.options.ConstrBoundary = 'penalize';
    problem.options.PopInitRange = [1, 3; repmat([-100, 100], 9, 1)];
end
