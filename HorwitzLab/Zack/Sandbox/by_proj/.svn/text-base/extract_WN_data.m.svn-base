function extract_WN_data(stro, spike_channel, output_dir)
if ~nargin, return; end
if nargin < 2 || isempty(spike_channel)
    spike_channel = 1;
end
if nargin < 3 || isempty(output_dir)
    output_dir = '';
else
    warning('OFF', 'MATLAB:MKDIR:DirectoryExists');
    mkdir(fullfile(output_dir));
    warning('ON', 'MATLAB:MKDIR:DirectoryExists');
end

nstixperside = stro.sum.exptParams.nstixperside;
ngammasteps = 2^16;
trial_fields = stro.sum.trialFields(1,:);
seedidx = strcmp(trial_fields, 'seed');
noisetypeidx = strcmp(trial_fields, 'noise_type');
nframesidx = strcmp(trial_fields, 'num_frames');
stimonidx = strcmp(trial_fields, 'stim_on');

sigmaidxs = strncmp(trial_fields, 'sigma', 5);
muidxs = strncmp(trial_fields, 'mu', 2);

spike_channels = {'sig001a' 'sig001b'};
spikeidx = strcmp(stro.sum.rasterCells(1,:), spike_channels{spike_channel});

gauss_locut = stro.sum.exptParams.gauss_locut;
gauss_hicut = stro.sum.exptParams.gauss_hicut;
msperframe = 1000/stro.sum.exptParams.framerate;

xx = linspace(gauss_locut/1000, gauss_hicut/1000, ngammasteps);
yy = norminv(xx');

L = stro.trial(:,noisetypeidx) == 1 & stro.trial(:,nframesidx) > 0;
stro.ras(~L,:) = [];
stro.trial(~L,:) = [];

[~,filename] = fileparts(stro.sum.fileName);
filepath = fullfile(output_dir, filename);
fraw = fopen([filepath '.raw'], 'wb');
fspiketrain = fopen([filepath '.isk'], 'wt');

c = onCleanup(@() close_files(fraw,fspiketrain));

ntrials = size(stro.trial,1);
for ii = 1:ntrials
    nframes = stro.trial(ii,nframesidx);
    seed = stro.trial(ii,seedidx);
    mu = stro.trial(ii,muidxs)/1000;
    sigma = stro.trial(ii,sigmaidxs)/1000;
    
    invnormcdf = bsxfun(@plus, bsxfun(@times, yy, sigma), mu);
    randnums = getEJrandnums(3*nstixperside^2*nframes, seed);
    randnums = reshape(randnums, [nstixperside^2*3 nframes]);
    for j = 1:3
        idxs = (1:nstixperside^2)+nstixperside^2*(j-1);
        randnums(idxs,:) = reshape(invnormcdf(randnums(idxs,:)+1,j), [nstixperside^2 nframes]);
    end
    
    t_stimon = stro.trial(ii, stimonidx);
    spiketimes = (stro.ras{ii,spikeidx}-t_stimon)*1000;
    spiketimes(spiketimes < 0) = [];
    frametimes = linspace(0, nframes*msperframe, nframes)+(msperframe/2);
    spiketimes(spiketimes > frametimes(end)) = [];
    spike_train = hist(spiketimes, frametimes);
    
    fprintf(fspiketrain, '%d\n', spike_train(:));
    fwrite(fraw, randnums(:), 'double');
end

function close_files(varargin)
for ii = 1:nargin
    try
        fclose(varargin{ii});
    catch
        continue
    end
end
