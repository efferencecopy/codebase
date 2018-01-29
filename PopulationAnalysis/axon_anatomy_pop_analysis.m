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
plot_z_profiles(dat, all_hvas);


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

NORMALIZE = true;
all_hvas = {'lm', 'al', 'pm', 'am'};
plot_avg_density_profile_within_mice(dat, all_hvas, NORMALIZE)



