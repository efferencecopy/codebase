% npixperstix nstixperside framerate gauss_locut gauss_hicut first_seed sum_nframes
% For loading one or two files.
spikenum = 1;
% fnames = {{'K010809002' 'K010809004' 'K010809008'}}; % Blue-yellow

% For loading a whole mess of files
[fnames, spikenum] = fnamesFromTxt2('C:\Documents and Settings\zacklb\Desktop\MatlabCode\Zack\Sandbox\by_proj\byopp.txt');

noisetype = 'GUN';
spike_channel_names = {'sig001a' 'sig001b'};
for ii = 1:length(fnames)
    stro = cellfun(@(x) {nex2stro(findfile(x))}, fnames{ii});
    stro = strocat(stro{:});
    
    nstixperside = stro.sum.exptParams.nstixperside;
    noisetypeidx = find(strcmp(stro.sum.trialFields(1,:), 'noise_type'));
    
    maxT = 9;
    Lgunnoise = stro.trial(:,noisetypeidx) == 1;
    Lconenoise = stro.trial(:,noisetypeidx) == 2;
    if (strcmp(noisetype, 'CONE'))
        stro.ras(Lgunnoise,:) = [];
        stro.trial(Lgunnoise,:) = [];
    elseif (strcmp(noisetype,'GUN'))
        stro.ras(Lconenoise,:) = [];
        stro.trial(Lconenoise,:) = [];
    end
    
    out = getWhtnsStats(stro, maxT, 'STCOVfull', {nstixperside^2, 3, 1, maxT}, spike_channel_names{spikenum(ii)});
    [~,filename] = fileparts(stro.sum.fileName);
    STCGUI(out{1}, out{2}, nstixperside, maxT, filename)
    
    if length(fnames) > 1, pause; end
end
