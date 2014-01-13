function p = paramsCheck(textFile)

%
% Should take a text file list of DT experiments and return the range of
% experimental parameters used in the experiments. The code will return:
%   p.rf       => the [x,y] of the stimuli
%   p.sigma    => the sd of the gabors
%   p.sf       => the spatial frequencies tested
%   p.theta    => the orientations tested
%   p.gamma    => the ratio of SDx and SDy
%   p.size     => the size of the gabor in number of standard deviations
%   p.speed    => speed in cyc/sec
%
% These entries include ONLY THE UNIQUE values. Thus if the RF pos was
% consistent across expts, than there will only be one entry
%
% CAH 11/9

fNames = fnamesFromTxt2(textFile);
nExpts = size(fNames,1);
p = struct('rf', {nan(nExpts,2)},...
           'colorDirs', {nan(nExpts, 9)},...
           'bkgndrgb', {nan(nExpts, 3)},...
           'sigma', {nan(nExpts,1)},...
           'sf', {repmat({[]}, nExpts,1)},...
           'theta', {nan(nExpts,1)},...
           'gamma', {nan(nExpts,1)},...
           'size', {nan(nExpts,1)},...
           'speed', {nan(nExpts,1)},...
           'length', {nan(nExpts,1)},...
           'frames', {nan(nExpts,1)});


for a = 1:nExpts
    disp(fNames{a}{1})
    fpath = findfile(fNames{a}{1}, '/Volumes/NO BACKUP/NexFiles/Greg/Apollo/')
    stro = nex2stro(fpath);
    DTindicies

    % find the parameters:
    p.rf(a,:) = [stro.sum.exptParams.rf_x, stro.sum.exptParams.rf_y];
    p.bkgndrgb(a,:) = [stro.sum.exptParams.bkgnd_r, stro.sum.exptParams.bkgnd_g, stro.sum.exptParams.bkgnd_b];
    p.colorDirs(a,:) = stro.sum.exptParams.RF_colors;
    p.sigma(a) = unique(stro.trial(:,gaborSigmaInd));
    p.theta(a) = unique(stro.trial(:, gaborThetaInd));
	p.gamma(a) = unique(stro.trial(:, gaborGammaInd));
    p.size(a) = stro.sum.exptParams.flash_size;
    p.speed(a) = unique(stro.trial(:,driftRateInd));
    p.length(a) = stro.sum.exptParams.flash_length;
    
    %gabor sf
    lambdas = unique(stro.trial(:, gaborLambdaInd));
    sfs = (1./lambdas) .* stro.sum.exptParams.pixperdeg;
	p.sf{a} = sfs;
    
    %wire frames
    frames = stro.trial(:,framesPresentInd);
    frames(isnan(frames)) = 1; %if the entry is a nan, that means there were frames.
    p.frames(a) = any(frames);
    
end
