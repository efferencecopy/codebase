%% Testing DTfindMechVect

fin

%define a bunch of things:
lvmUnit = [1/sqrt(2), -1/sqrt(2), 0];
sUnit = [0 0 1];
swmUnit = [0.0636, -0.0636, -0.45];
swmUnit = swmUnit./norm(swmUnit);
swlUnit = swmUnit .* [1 1 -1];

%open the figure and hold onto the handles
f = figure;
s(1) = subplot(1,2,1);
s(2) = subplot(1,2,2);

% specify the actual mechanism vector and two color directions to probe
mechTheta = linspace(0,pi,15);
for t = 1:length(mechTheta)
    mechUnit = (cos(mechTheta(t))*lvmUnit) + (sin(mechTheta(t))*sUnit);
    cardUnit = lvmUnit;
    intUnit = swmUnit;
    beta = [-.2 35]; % just an offset and contrast term
    
    
    % make a CRF where the log(counts) is a linear function of contrast
    cardNorms = logspace(log10(0.04), log10(0.11), 7);
    cardLMS = cardNorms(:) * cardUnit;
    cardProj = abs(cardLMS * mechUnit(:));
    
    intNorms = logspace(log10(.04), log10(0.11), 7);
    intLMS = intNorms(:) * intUnit;
    intProj = abs(intLMS * mechUnit(:));
    
    %turn the projections into firing rates using an exponential and some of
    %the experimentally derived betas
    cardCounts = glmval(beta', cardProj(:), 'log');
    intCounts = glmval(beta', intProj(:), 'log');
    counts = ([cardCounts;intCounts]);
    
    
    %run the analysis
    [guesses, exitVals, flags] = DTfindMechVect([cardLMS;intLMS], counts, 'deviance');
    
    exitVals
    [~, idx] = min(exitVals);
    mechGuess = guesses(idx,:);
    [mechUnit ; mechGuess]
    
    
    subplot(s(1)), cla, hold on
    plot(abs(cardLMS*cardUnit(:)), cardCounts, 'k');
    plot(abs(intLMS*cardUnit(:)), intCounts, 'b');
    plot(abs(intLMS*mechGuess(:)), intCounts, 'bo:');
    plot(abs(cardLMS*mechGuess(:)), cardCounts, 'ko:')
    set(gca, 'xscale', 'log')
    hold off
    
    
    %now compute the deviance at a bunch of potential mechanism directions
    colorDirs = linspace(0,pi,1000);
    dev = nan(1, numel(colorDirs));
    for a = 1:numel(colorDirs);
        clrVec = (cos(colorDirs(a))*lvmUnit(:)) + (sin(colorDirs(a))*sUnit(:));
        clrVec = clrVec./norm(clrVec); %make it a unit vector
        proj = abs([cardLMS;intLMS]*clrVec(:));
        
        %mechinism direction case
        [mechBeta, dev(a)] = glmfit(proj, counts, 'poisson', 'link', 'log');
    end
    
    
    subplot(s(2)), cla, hold on,
    cardTheta = mod(atan(cardUnit(3)./(cardUnit*lvmUnit(:))), pi);
    intTheta = mod(atan(intUnit(3)./(intUnit*lvmUnit(:))), pi);
    mechGuessTheta = mod(atan(mechGuess(3)./(mechGuess*lvmUnit(:))), pi);
    plot(colorDirs*(180/pi), dev);
    plot(cardTheta*(180/pi), 0, 'ko', 'markersize', 15)
    plot(intTheta*(180/pi), 0, 'bo', 'markersize', 15)
    plot(mechTheta(t)*(180/pi), 0, 'm^', 'markersize', 15)
    plot(mechGuessTheta*(180/pi),0,'rV', 'markersize', 15);
    hold off
    
    pause(2);
end


