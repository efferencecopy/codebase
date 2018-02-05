%% TODO LIST
 
% make sure that the imaging parameters are consistent within a mouse

% make sure that the um/pix are the same for ALL images and ALL mice. This
% is important b/c my analysis is based on pixels. If the imaging params
% are different from mouse to mouse or HVA to HVA than I can't compare
% fairly across conditions

% make sure I'm aggregating fluorescence correctly. Currently(? asof 1/2018) I'm
% computing the mean

% figure out how to take the correct number of z-planes

% eliminate AK090314_A, AK_090314_C due to bi-modality in the z-plane.
% probably due to inappropriate confocal calibration. MAYBE: they can be
% used for relative depths of expression, but they maybe shouldn't be used
% for absolute comparisons between HVAs? Although perhaps the error is
% systematic enough that it doesn't matter...

% subtract off the auto-fluorescence. Since this value will change (if we
% modify how fluorescence is aggregated across x/y) than perhaps I should
% store the profile coordinates that map onto the area where
% auto-fluorescence can be calculated. OR: should I subtract off a linear
% trendline from the average profiles? Should I make this subtraction
% common across all HVAs within a mouse?

%% LOAD IN ALL THE DATA

fin
dat = load_all_data();
save \\crash.dhe.duke.edu\charlie\pre_processed_data\axon_density_preprocessed.mat dat -v7.3
fprintf('All done importing the data\n')


%% LOAD AN EXISTING (PRE-PROCESSED) .MAT FILE

% the load command defines 'dat', which is a structure array
fin
load('\\crash.dhe.duke.edu\charlie\pre_processed_data\axon_density_preprocessed.mat')


%% QC PLOT: Z-PROFILES

close all; clc
all_hvas = {'lm', 'pm', 'al', 'am'};
NORMALIZE = false;
plot_z_profiles(dat, all_hvas, NORMALIZE);


%% PLOT THE RAW IMAGES WITH LAYER BOUNDARIES

close all; clc
all_hvas = {'lm', 'pm', 'al', 'am'};
plot_raw_images_with_layer_boundaries(dat, all_hvas);


%% PROFILE OF FLUORESCENCE INTENSITY (ONE PER IMAGE)

close all; clc
all_hvas = {'lm', 'pm', 'al', 'am'};
plot_density_profiles_from_all_slices(dat, all_hvas);


%% AVERAGE FLUORESCENCE INTENSITY (WITHIN MOUSE)

close all; clc

NORMALIZE = false;
all_hvas = {'lm', 'al', 'pm', 'am'};
plot_avg_density_profile_within_mice(dat, all_hvas, NORMALIZE)



