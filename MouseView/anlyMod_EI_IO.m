function params = anlyMod_EI_IO(params)


% Game plan:
%
% 1) Several files could be in the queue for analysis. They should already
% be abf2obj'd, but they will likely have an unknown number of interleaved
% stimulus types. I need some way of account for all the different stimulus
% types. 
%
%