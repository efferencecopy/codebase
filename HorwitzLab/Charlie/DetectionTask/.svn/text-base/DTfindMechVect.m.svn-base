function [mechVects, exitVals, exitFlag] = DTfindMechVect(LMS, counts, fitFxn)

%
%  DTfindMechVect
%    EXAMPLE: [mechVect, exitFlag] = DTfindMechVect(LMS, counts, fitFxn);
%
% This function is designed to take in CRF data from DT MoCS experiments
% and estimate the vector in color space that best predicts the joint
% Cardinal and Intermediate CRFs. Analytically there should be two color
% directions that predict the data equally well, so mechVects is a 2x3
% matrix (one row for each color direction). fitFxn can be either 'sse', or
% 'deviance'

% CAH 08/11
%

%define some parameters that are called by subfunctions:
lvmUnit = [1 -1 0]./norm([1 -1 0]);
sUnit = [0 0 1];
options.MaxIter = 10000;
options.MaxFunEvals = 10000;


% %clean up the inputs and remove the zero contrast conditions...
% l_zero = sum(abs(LMS),2)==0;
% counts(l_zero) = [];
% LMS(l_zero,:) = [];


%determine the two lowest points on the error surface and use them as
%initial guesses to the fminsearch routine.
nThetas = 500;
thetas = mat2cell(linspace(0, pi, nThetas)', ones(nThetas,1));

%find the minima of the error function
errSurf = cellfun(@minerr, thetas);
dx = sign(diff(errSurf)); %looking for changes in the sign of the first derivitive...
idx = diff(dx) == 2;
idx = [false; false; idx]; %prepending 2 zeros b/c I called diff twice


% minima at the edges won't get counted by default, so make sure to include
% them:
if (sign(dx(2)) - sign(dx(end))) == 2;
    idx(1) = true;
end

%assign the guesses. 
thetas = cell2mat(thetas);
initGuess = thetas(idx);
if numel(initGuess) > 2
    errs = errSurf(idx);
    [~, inds] = sort(errs);
    initGuess = initGuess(inds); %order them from lowest to highest deviance values;
    initGuess = initGuess(1:2); %just take the lowest two.
end

%run the minimization
for a = 1:length(initGuess)
    [mechTheta(a), exitVals(a), exitFlag(a)] = fminsearch(@minerr, initGuess(a), options);
end

%assign the outputs
mechVects = (cos(mechTheta)'*lvmUnit) + (sin(mechTheta)'*sUnit);
mechVects = bsxfun(@rdivide, mechVects, sqrt(sum(mechVects.^2,2))); %make sure they are unit vectors

if size(mechVects,1) <2;
    keyboard
end

% NESTED SUBFUNCTIONS %


    function err = minerr(theta)
        clrVec = (cos(theta)*lvmUnit(:)) + (sin(theta)*sUnit(:));
        clrVec = clrVec./norm(clrVec); %make it a unit vector
        proj = abs(LMS*clrVec(:));
        
        %mechinism direction case
        [beta, dev] = glmfit(proj, counts, 'poisson', 'link', 'log');
        
        switch lower(fitFxn)
            case {'deviance'}
                err = dev;
            case 'sse'
                mod = glmval(beta, proj, 'log');
                err = sum((mod-counts).^2);
        end
    end


end