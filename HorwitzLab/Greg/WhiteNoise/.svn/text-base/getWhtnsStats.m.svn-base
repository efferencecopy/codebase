function out = getWhtnsStats(stro,nframesback,module, initargs, spikename)
% out = getWhtnsStats(stro,maxT,module,initargs, <spikename>)
%
%	 Based on getWhtnsStats3, this function will iteratively call a "module"
% function with vectors of RBGs and spikes (once per trial).  This is a
% work in progress.
%
%    Initargs should be a cell array containing the initialization
%    arguments appropriate for the module called.
%
% GDLH 3/5/08

nstixperside = stro.sum.exptParams.nstixperside; %get number of stixels per side
gammaTable = stro.sum.exptParams.gamma_table; %get gamma_table
gammaTable = reshape(gammaTable,length(gammaTable)/3,3); %reshapse gamma_table into three columns
invgammaTable = InvertGamma(gammaTable,1); %invert gamma_table
ngammasteps = size(invgammaTable,1); %get number of rows of gamma_table (65536)
seedidx = strcmp(stro.sum.trialFields(1,:),'seed'); %get seed index from trialFields
noisetypeidx = strcmp(stro.sum.trialFields(1,:),'noise_type'); %get noise_type index from trialFields
nframesidx = strcmp(stro.sum.trialFields(1,:),'num_frames'); %get nframes index from trialFields
stimonidx = strcmp(stro.sum.trialFields(1,:),'stim_on'); %get stimon index from trialFields
muidxs = [find(strcmp(stro.sum.trialFields(1,:),'mu1')),... %get mu indices into vector from trialFields
         find(strcmp(stro.sum.trialFields(1,:),'mu2')),...
         find(strcmp(stro.sum.trialFields(1,:),'mu3'))];
sigmaidxs = [find(strcmp(stro.sum.trialFields(1,:),'sigma1')),... %get sigma indices into vector from trialFields
            find(strcmp(stro.sum.trialFields(1,:),'sigma2')),...
            find(strcmp(stro.sum.trialFields(1,:),'sigma3'))];
% bkgndRGBidxs = [find(strcmp(stro.sum.trialFields(1,:),'bkgnd_r')),...
%                find(strcmp(stro.sum.trialFields(1,:),'bkgnd_g')),...
%                find(strcmp(stro.sum.trialFields(1,:),'bkgnd_b'))];
if (nargin < 5)
    spikeidx = find(strcmp(stro.sum.rasterCells(1,:),'sig001a')); %get spike index from rasterCells
else
    spikeidx = find(strcmp(stro.sum.rasterCells(1,:),spikename)); %get spike index from rasterCells
end
msperframe = 1000/stro.sum.exptParams.framerate; %calculate msec per frame
ntrials = size(stro.trial,1); %get number of trials

% compute once for gun noise
xx = linspace(stro.sum.exptParams.gauss_locut/1000, stro.sum.exptParams.gauss_hicut/1000, ...
    ngammasteps); %dividing gauss by 1000 and making equal intervals so that there are 65536 values
yy = norminv(xx');

feval(module,'init',initargs);

hwait = waitbar(0,'Finding stimuli...'); %display progress bar
for i = 1:ntrials %get values and insert into given column
    nframes = stro.trial(i,nframesidx);
    if nframes == 0, continue; end
    seed = stro.trial(i,seedidx);
    noisetype = stro.trial(i,noisetypeidx);
    mu = stro.trial(i,muidxs)/1000;
    sigma = stro.trial(i,sigmaidxs)/1000;
    
    if (noisetype == 1)  % Gaussian gun
        invnormcdf = bsxfun(@plus, bsxfun(@times, yy, sigma), mu);
        randnums = getEJrandnums(3*nstixperside^2*nframes, seed);
        randnums = reshape(randnums, [nstixperside^2*3, nframes]);
        for j = 1:3
            idxs = (1:nstixperside^2)+nstixperside^2*(j-1);
            randnums(idxs,:) = reshape(invnormcdf(randnums(idxs,:)+1,j),[length(idxs),nframes]);
        end
    elseif (noisetype == 2)  % Binary cone
        colordirlms = sign(fullfact([2,2,2])-1.5);
        randnums = getEJrandnums(nstixperside^2*nframes, seed);
        randnums = reshape(randnums, [nstixperside^2, nframes]);
        randnums = mod(randnums, size(colordirlms,1))+1;
        tmp = colordirlms(randnums,:);
        tmp = reshape(tmp, [nstixperside^2 nframes 3]);
        tmp = permute(tmp,[1 3 2]);
        randnums = reshape(tmp,[nstixperside^2*3 nframes]);
    else
        continue;
    end
    t_stimon = stro.trial(i, stimonidx);
    spiketimes = (stro.ras{i,spikeidx}-t_stimon)*1000;  % converting to ms
    frametimes = linspace(0, nframes*msperframe, nframes)+(msperframe/2)';
    spiketimes(spiketimes < nframesback*msperframe) = [];
    spiketimes(spiketimes > frametimes(end)) = [];
    n = hist(spiketimes, frametimes);
    feval(module,randnums(:),n);
    if (ishandle(hwait))
        waitbar(i/ntrials,hwait);
    else
        break;
    end
end
if (ishandle(hwait))
   close(hwait);
   out = feval(module,'return');
end
eval(['clear ',module]);
