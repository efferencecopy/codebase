% Fitting a gabor to an STA.  This may be useful for a "foward
% correlation" style experiment for testing the idea that luminance edges
% facilitate the responses of color-opponent V1 neurons and maybe also be
% useful for setting up stimuli in the detection task.
%
% GDLH 8/1/08

% Evaluate this stuff to load in an STA to work on
%stro=nex2stro('C:\NexFiles\Kali\K060508001.nex');

cd('C:\NexFiles\Kali') %go to the directory with all the nex files
fnames = fnamesFromTxt('Alltest.txt'); %load a matrix of names to the expts
numExpts = size(fnames, 1);
ExitFlagVect = nan*ones(numExpts,1);
for a = 1:numExpts;
    stro = nex2stro([fnames(a, :),'.nex']);
    disp(stro.sum.fileName)
    nstixperside = stro.sum.exptParams.nstixperside;
    noisetypeidx = find(strcmp(stro.sum.trialFields(1,:),'noise_type'));
    spikeidx = find(strcmp(stro.sum.rasterCells(1,:),'sig001a'));

    npixperstix = 3;  % This is given by REX
    nstixperside = 10; % This is given by REX

    % Getting the STA
    maxT = 9;
    tmpstro = stro;
    Lgunnoise = stro.trial(:,noisetypeidx) == 1;
    Lconenoise = stro.trial(:,noisetypeidx) == 2;
    tmpstro.ras(Lconenoise,:) = [];
    tmpstro.trial(Lconenoise,:) = [];
    out = getWhtnsStats(tmpstro,maxT,'STCOVmex', {nstixperside^2, 3, maxT}, spikeidx);
    tmpstro = [];
    STAs = out{1};
    STCs = out{2};
    nspikes = out{3};
    peakframe = find(sum(STAs.^2)==max(sum(STAs.^2)));
    im = reshape(STAs(:,peakframe),[size(STAs,1)/3, 3]);
    [u,s,v] = svd(im');
    if (max(abs(v(:,1))) ~= max(v(:,1)))
        v = -v;
        u = -u;
    end
    gunweights = u(:,1);
    im = reshape(v(:,1),[sqrt(size(STAs,1)/3),sqrt(size(STAs,1)/3)]);
    
    out = gaborfit(im); 
    im = im./max(abs(im(:)));

    % Getting initguess parameters
    imsq = im.^2;
    interval = [1:nstixperside];
    interval = interval-ceil(median(interval));
    initguess.xoffset = (interval*sum(imsq)')./sum(imsq(:));
    initguess.yoffset = -(interval*sum(imsq')')./sum(imsq(:));

    % Looking at the power spectrum
    powerspectra = abs(fftshift(fft2(im))).^2;
    powerspectra = powerspectra .* (powerspectra > .1*max(powerspectra(:)));
    powerspectra = powerspectra./sum(powerspectra(:));
    phases = angle(fftshift(fft2(im)));

    [i,j] = meshgrid([1:10]-6,[1:10]-6);
    angles = atan2(i,j)-pi/2;
    amp = sqrt(i.^2+j.^2);
    [peaksfy,peaksfx] = ind2sub(size(powerspectra),find(powerspectra == max(powerspectra(:)),1));
    initguess.lambda = nstixperside./sum(sum(amp.*powerspectra));
    cosines = powerspectra.*cos(angles);
    sines = powerspectra.*sin(angles);
    if (peaksfy<=6)
        meansin = mean(mean(sines(1:6,:)));
    else
        meansin = mean(mean(sines(6:10,:)));
    end
    if (peaksfx<=6)
        meancos = mean(mean(cosines(:,1:6)));
    else
        meancos = mean(mean(cosines(:,6:10)));
    end
    initguess.theta = atan2(meansin,meancos);
    initguess.theta = mod(initguess.theta + pi/2,2*pi);

    if (initguess.lambda>2*nstixperside)
        initguess.lambda = 5*nstixperside;
        initguess.sigma = sqrt(sum(abs(im(:)) > 0.5));
        initguess.gamma = 1;
        initguess.phi = 0;
    else
        initguess.sigma = initguess.lambda/3.5;
        initguess.gamma = 1.4;
        proj = [initguess.xoffset initguess.yoffset]*[cos(pi/2-initguess.theta); -sin(pi/2-initguess.theta)];
        [i,j] = ind2sub(size(powerspectra),find(max(powerspectra(:))==powerspectra,1));
        initguess.phi = phases(i,j);
        err1 = gaborfiterr([initguess.theta,initguess.lambda,initguess.phi,initguess.sigma,initguess.gamma,initguess.xoffset,initguess.yoffset],im);
        err2 = gaborfiterr([initguess.theta,initguess.lambda,initguess.phi+pi,initguess.sigma,initguess.gamma,initguess.xoffset,initguess.yoffset],im);
        if (err2<err1)
            initguess.phi = initguess.phi+pi;
        end
    end

    % Setting up a lattice for sampling the gabor
    interval = [1:nstixperside*2*npixperstix];
    interval = interval-ceil(median(interval));
    [X, Y] = meshgrid(interval(1:2:end),interval);

    % Sampling the initial guess gabor
    X = X-initguess.xoffset*2*npixperstix;
    Y = Y+initguess.yoffset*2*npixperstix;
    xprime = X.*cos(-initguess.theta)+Y.*sin(-initguess.theta);
    yprime = -X.*sin(-initguess.theta)+Y.*cos(-initguess.theta);
    initguess.gabor = exp(-(xprime.^2+initguess.gamma.^2.*yprime.^2)./(2.*(initguess.sigma*2*npixperstix).^2)).* ...
        cos(2.*pi.*yprime./(initguess.lambda*2*npixperstix)-initguess.phi);

    % Preparing for fitting
    theta = initguess.theta;
    lambda = initguess.lambda;  % pixels per cycle
    phi = initguess.phi;
    sigma = initguess.sigma;
    gamma = initguess.gamma;
    xoffset = initguess.xoffset;
    yoffset = initguess.yoffset;

    % Doing the fitting
    options = optimset('MaxIter',2000,'MaxFunEvals',3000,'Display','final');
    [x,fval,exitflag] = fminsearch(@(x)gaborfiterr(x,im),[theta,lambda,phi,sigma,gamma,xoffset,yoffset],options);
    theta = x(1);
    lambda = x(2)*2*npixperstix;
    phi = x(3);
    sigma = x(4)*2*npixperstix;
    gamma = x(5);
    xoffset = x(6)*2*npixperstix;
    yoffset = x(7)*2*npixperstix;

    % Calculating the fitted gabor
    [X, Y] = meshgrid(interval(1:2:end),interval);
    X = X-xoffset; Y = Y+yoffset; % So negative numbers mean down
    xprime = X.*cos(-theta)+Y.*sin(-theta);
    yprime = -X.*sin(-theta)+Y.*cos(-theta);
    fittedgabor = exp(-(xprime.^2+gamma.^2.*yprime.^2)./(2.*sigma.^2)).*cos(2.*pi.*yprime./lambda-phi);

    % Doing the plotting
    figure; set(gcf,'DefaultAxesXTick',[],'DefaultAxesYtick',[]);
    subplot(2,3,1); imagesc(im); colormap(gray); axis square; title('Raw');
    subplot(2,3,2); imagesc(initguess.gabor); colormap(gray); axis square; title('InitGuess');
    subplot(2,3,3); imagesc(fittedgabor); colormap(gray); axis square; title('Fitted');
    % Doing it in color now
    subplot(2,3,4);
    STA = reshape(STAs(:,peakframe),[10 10 3]);
    STA = STA./(2.1*max(abs(STA(:))))+0.5;
    image(STA);
    axis square;
    nxpix = size(initguess.gabor,1);
    nypix = size(initguess.gabor,2);
    subplot(2,3,5);
    colorgabor = gunweights*reshape(initguess.gabor,1,nxpix*nypix);
    colorgabor = reshape(colorgabor,[3, nxpix, nypix]);
    colorgabor = permute(colorgabor,[2 3 1]);
    colorgabor = colorgabor./(2.1*max(abs(colorgabor(:))))+0.5;
    image(colorgabor); axis square;
    subplot(2,3,6);
    colorgabor = gunweights*reshape(fittedgabor,1,nxpix*nypix);
    colorgabor = reshape(colorgabor,[3, nxpix, nypix]);
    colorgabor = permute(colorgabor,[2 3 1]);
    colorgabor = colorgabor./(2.1*max(abs(colorgabor(:))))+0.5;
    image(colorgabor); axis square;
    title([stro.sum.fileName,': ',num2str(exitflag)]);
    ExitFlagVect(a) = exitflag;
end

