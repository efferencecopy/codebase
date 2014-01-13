load wn_params
data = nan(size(WNfiles,1), 7);
seedidx = 5; nframesidx = 12; noiseidx = 16;
for ii = 1:length(WNfiles)
    [~,filename] = fileparts(WNfiles{ii});
    if filename(1) == 'A' || filename(1) == 'F', continue; end
    try
        stro = nex2stro(WNfiles{ii});
    catch
        delete(findall(0,'type','figure','tag','TMWWaitbar'));
        continue
    end
    Lgun = stro.trial(:,noiseidx) == 1;
    stro.trial(~Lgun,:) = [];
    if ~isempty(stro.trial) && length(stro.sum.rasterCells) < 10 && stro.trial(end,seedidx) ~= 1
        data(ii,:) = [length(stro.trial) stro.sum.exptParams.npixperstix stro.sum.exptParams.nstixperside stro.sum.exptParams.framerate ...
            stro.trial(1,seedidx) stro.trial(end,seedidx) sum(stro.trial(:,nframesidx))];
    end
end
