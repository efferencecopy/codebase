function spikename = getSpikenum(stro,opt)
%
% If the optional second argument is set to 'first'
% we return the first spike channel.  Otherwise we present a 
% dialog box.

if (nargin < 2)
    opt = 0;
end
choices = {};
for i = 1:length(stro.sum.rasterCells)
    if (strncmp(stro.sum.rasterCells(i),'sig0',4))
        choices{length(choices)+1} = stro.sum.rasterCells{i};
    end
end

%weed out the lfp channels (if present)
lfpChannels = strfind(choices, '_wf');
lfpChannels = cellfun(@any, lfpChannels);
choices(lfpChannels) = [];

%present a dialog box.
if (strcmp(opt,'first') || length(choices) == 1)
    spikename = char(choices{1});
else
    a = listdlg('PromptString','Which spike?','SelectionMode','Single','ListSize',[150 50],'ListString',choices);
    if (isempty(a))
        a = 1;
    end
    spikename = char(choices{a});
end