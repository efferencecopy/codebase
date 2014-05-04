% NOTES
%
% ### TO DO
%
% Export the PSCs for each stimulus location (ordered in the same fashion
% as the file names) and the location matrix itself.
%
% Save the params stuff in the MDB. Do some error checking to make certian
% that the stuff in the params field lands in the correct entry on the mdb.
% Also, I shouldn't save the abfobj's becuase they are redundant with the
% raw data, so I should flag this and omit it from the mdb. Otherwise the
% mdb will get huge.
%
% Add the analysis tags to the searchable text fields? (this would be
% overwritten when new files are added or when updates are folded into the
% mdb.
%
% Make a text file. It should document each HVA and neurons that are
% located there. Maybe make a different text file for L2/3 PY cells, L2/3
% FS cells.
%
% For now, the text file should include
% * mouse name
% * cell number
% * stimulation location
% * HVA location
%
%
%
% ### ANALYSIS
%
% Look at norm PSC amp train for each TF. All cells plotted on the same set
% of axes (one axis for each TF, one figure for each HVA)
%
% Look at PPratio P1:P2 as function of TF. All cells get plotted on the
% same set of axes (one axis for each HVA).
%
% Head to head comparison for HVAs: plot the average PPratio vs TF for each
% HVA on the same set of axes.
%
% Some how compare the effect of different stimulus location so that we can
% figure out where to stimulate (instead of stimulating at the soma)
%
% Plot the PN:P1 ratio for each TF, fit an exponential to the data sets
% where there is only monotonic depression. Compare tau's across HVAs












