function params = anlyMod_EI_IO(params)


% Game plan:
%
% 1) Several files could be in the queue for analysis. They should already
% be abf2obj'd, but they will likely have an unknown number of interleaved
% stimulus types. I need some way of account for all the different stimulus
% types.
%
%



% iterate over the files.
%
% outerleave to retrieve the number of conditions
%
% go through each condition CH by CH. the list of valid channels should be
% given by the outerleave tdict AND by the params.exclude field.
%
% calculate the average of each group, store it in the
% params.isolatedCurrents field along with the holding potential.
%
% plot the raw values (baseline subtracted) and the average.
%
% calculate conductance on the fly (also send this to the isolatedData
% fields
%
% calculate peaks on a peak by peak basis


for i_ax = 1:numel(params.ax)
   
   ledIdx = params.ax{i_ax}.idx.LED_470;
   tdict = outerleave(params.ax{i_ax}, ledIdx);
    
end



