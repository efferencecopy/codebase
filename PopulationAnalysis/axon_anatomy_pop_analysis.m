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

% fin
% dat = load_all_data();
save \\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\charlie\pre_processed_data\axon_density_preprocessed_13z_tlx_emx_only_baselined.mat dat -v7.3
fprintf('All done importing the data\n')


%% LOAD AN EXISTING (PRE-PROCESSED) .MAT FILE

% the load command defines 'dat', which is a structure array
fin
load('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\charlie\pre_processed_data\axon_density_preprocessed_13z_tlx_emx_only_baselined.mat')


%% QC PLOT: Z-PROFILES

close all; clc
all_hvas = {'lm', 'pm', 'al', 'am'};
NORMALIZE = true;
plot_z_profiles(dat, all_hvas, NORMALIZE);

% now plot the qc profiles individually
plot_z_profiles_for_each_mouse(dat, all_hvas, NORMALIZE);


%% PLOT THE RAW IMAGES WITH LAYER BOUNDARIES

close all; clc
all_hvas = {'lm', 'pm', 'al', 'am'};
plot_raw_images_with_layer_boundaries(dat, all_hvas);


%% PROFILE OF FLUORESCENCE INTENSITY (ONE PER IMAGE)

close all; clc
all_hvas = {'lm', 'pm', 'al', 'am'};
plot_density_profiles_from_all_slices(dat, all_hvas);



%% SUBTRACT OFF AUTOFLUROESCENCE

close all; clc
all_hvas = {'lm', 'al', 'pm', 'am'};
dat = subtract_autofluorescence(dat, all_hvas);


%% AVERAGE FLUORESCENCE INTENSITY (WITHIN MOUSE)

close all; clc

NORMALIZE = 'max';
PROFILE_TYPE = 'baselined';
all_hvas = {'lm', 'al', 'pm', 'am'};
plot_avg_density_profile_within_mice(dat, all_hvas, NORMALIZE, PROFILE_TYPE);


%% AVERAGE FLUORESCENCE INTENSITY (ACROSS MICE)

close all; clc

NORMALIZE = 'max';
PROFILE_TYPE = 'baselined';
all_hvas = {'lm', 'al', 'pm', 'am'};

hva_aggs = aggregate_data_across_mice(dat, all_hvas, NORMALIZE, PROFILE_TYPE);

plot_avg_density_profile_across_mice(hva_aggs, all_hvas);


%% HEAT MAP OF RELATIVE STRENGTH OF AXON INPUTS

close all; clc

all_hvas = {'lm', 'al', 'pm', 'am'};
dat = get_average_layer_boundaries(dat, all_hvas);

% do not normalize the data to anything (will happen in plot function)
hva_aggs = aggregate_data_across_mice(dat, all_hvas, 'none', 'baselined');

% make the heat maps.
mline_data = plot_absolute_axon_density_across_mice(hva_aggs, all_hvas, 'self');
mline_data = plot_absolute_axon_density_across_mice(hva_aggs, all_hvas, 'lm');


%% STATS: EMX layer specificity Oneway Anova with tukey

close all; clc

mline_data = plot_absolute_axon_density_across_mice(hva_aggs, all_hvas, 'self');

emx_data = mline_data.emx; % N_layers x N_hvas x N_mice
emx_data = permute(emx_data, [3,2,1]); % N_mice x N_hvas x N_layers

n_expected_rows = size(emx_data, 1) * size(emx_data, 2);
n_layers = 5;
emx_data = reshape(emx_data, n_expected_rows, n_layers);

[~,~,STATS] = anova1(log10(emx_data), {'L1', 'L2/3', 'L4', 'L5', 'L6'}, 'on');


comparison = multcompare(STATS);

%% STATS: EMX HVA specificity.

close all; clc

mline_data = plot_absolute_axon_density_across_mice(hva_aggs, all_hvas, 'lm');

emx_data = mline_data.emx; % N_layers x N_hvas x N_mice
emx_data = permute(emx_data, [3,1,2]); % N_mice x N_layers x N_hvas

n_expected_rows = size(emx_data, 1) * size(emx_data, 2);
N_hvas = length(all_hvas);
emx_data = reshape(emx_data, n_expected_rows, N_hvas);

[~,~,STATS] = anova1(log10(emx_data), all_hvas, 'on');
multcompare(STATS);


%% STATS:EMX two way anova

close all; clc

all_hvas = {'lm', 'al', 'pm', 'am'};
dat = get_average_layer_boundaries(dat, all_hvas);
hva_aggs = aggregate_data_across_mice(dat, all_hvas, 'none', 'baselined');
mline_data = plot_absolute_axon_density_across_mice(hva_aggs, all_hvas, 'lm');


N_hvas = length(all_hvas);
N_layers = 5;
layer_mtx = repmat({'l1'; 'l23'; 'l4'; 'l5'; 'l6'}, 1, N_hvas);
hva_mtx = repmat(all_hvas, N_layers, 1);

% deal with the EMX data
emx_data = mline_data.emx; % N_layers x N_hvas x N_mice
emx_group_layers = repmat(layer_mtx, 1, 1, size(emx_data, 3));
emx_group_hvas = repmat(hva_mtx, 1, 1, size(emx_data, 3));
emx_group_genotype = repmat({'emx'}, size(emx_data));


% create 2 cols (lat, med)
Y = emx_data(:);
group_layers = emx_group_layers(:);
group_hvas = emx_group_hvas(:);

group = {group_hvas, group_layers};


[P, T, STATS, TERMS] = anovan(Y, group, 'varnames', {'hva', 'layers'});
figure
multcompare(STATS, 'dimension', [1,2])

%% STATS:TLX vs. EMX

close all; clc

all_hvas = {'lm', 'al', 'pm', 'am'};
dat = get_average_layer_boundaries(dat, all_hvas);
hva_aggs = aggregate_data_across_mice(dat, all_hvas, 'none', 'baselined');
mline_data = plot_absolute_axon_density_across_mice(hva_aggs, all_hvas, 'lm');


N_hvas = length(all_hvas);
N_layers = 5;
layer_mtx = repmat({'l1'; 'l23'; 'l4'; 'l5'; 'l6'}, 1, N_hvas);
hva_mtx = repmat(all_hvas, N_layers, 1);

% deal with the EMX data
emx_data = mline_data.emx; % N_layers x N_hvas x N_mice
emx_group_layers = repmat(layer_mtx, 1, 1, size(emx_data, 3));
emx_group_hvas = repmat(hva_mtx, 1, 1, size(emx_data, 3));
emx_group_genotype = repmat({'emx'}, size(emx_data));

% deal with the TLX data
tlx_data = mline_data.tlx3; % N_layers x N_hvas x N_mice
tlx_group_layers = repmat(layer_mtx, 1, 1, size(tlx_data, 3));
tlx_group_hvas = repmat(hva_mtx, 1, 1, size(tlx_data, 3));
tlx_group_genotype = repmat({'tlx'}, size(tlx_data));


% create 2 cols (lat, med)
Y = cat(1, emx_data(:), tlx_data(:));
group_layers = cat(1, emx_group_layers(:), tlx_group_layers(:));
group_hvas = cat(1, emx_group_hvas(:), tlx_group_hvas(:));
group_genotype = cat(1, emx_group_genotype(:), tlx_group_genotype(:));

group = {group_genotype, group_hvas, group_layers};


[P, T, STATS, TERMS] = anovan(Y, group, 'varnames', {'genotype', 'hva', 'layers'});
figure
multcompare(STATS)


%% STATS: TLX HVA specificity.

close all; clc
all_hvas = {'lm', 'al', 'pm', 'am'};
mline_data = plot_absolute_axon_density_across_mice(hva_aggs, all_hvas, 'lm');

tlx3_data = mline_data.tlx3; % N_layers x N_hvas x N_mice
tlx3_data = permute(tlx3_data, [3,1,2]); % N_mice x N_layers x N_hvas

n_expected_rows = size(tlx3_data, 1) * size(tlx3_data, 2);
N_hvas = length(all_hvas);
tlx3_data = reshape(tlx3_data, n_expected_rows, N_hvas);

[~,~,STATS] = anova1(tlx3_data, all_hvas, 'on');
multcompare(STATS);




%% STATS: EMX layer by layer

close all; clc

layer_to_analyze = 3;

all_hvas = {'lm', 'al', 'pm', 'am'};
dat = get_average_layer_boundaries(dat, all_hvas);
hva_aggs = aggregate_data_across_mice(dat, all_hvas, 'none', 'baselined');
mline_data = plot_absolute_axon_density_across_mice(hva_aggs, all_hvas, 'lm');

N_layers = 1;
hva_mtx = repmat(all_hvas, N_layers, 1);

% deal with the EMX data
emx_data = mline_data.emx(layer_to_analyze, :, :); % 1_layer x N_hvas x N_mice
emx_group_hvas = repmat(hva_mtx, 1, 1, size(emx_data, 3));

% create 2 cols (lat, med)
Y = emx_data(:);
group_hvas = emx_group_hvas(:);

group = {group_hvas};


[P, T, STATS, TERMS] = anovan(Y, group, 'varnames', {'hva'});
figure
multcompare(STATS)


%% STATS: TLX vs. EMX main effect of layers

close all; clc

all_hvas = {'lm', 'al', 'pm', 'am'};
dat = get_average_layer_boundaries(dat, all_hvas);
hva_aggs = aggregate_data_across_mice(dat, all_hvas, 'none', 'baselined');
mline_data = plot_absolute_axon_density_across_mice(hva_aggs, all_hvas, 'lm');


N_hvas = length(all_hvas);
N_layers = 5;
layer_mtx = repmat({'l1'; 'l23'; 'l4'; 'l5'; 'l6'}, 1, N_hvas);
hva_mtx = repmat(all_hvas, N_layers, 1);

% deal with the EMX data
emx_data = mline_data.emx; % N_layers x N_hvas x N_mice
emx_group_layers = repmat(layer_mtx, 1, 1, size(emx_data, 3));
emx_group_hvas = repmat(hva_mtx, 1, 1, size(emx_data, 3));
emx_group_genotype = repmat({'emx'}, size(emx_data));

% deal with the TLX data
tlx_data = mline_data.tlx3; % N_layers x N_hvas x N_mice
tlx_group_layers = repmat(layer_mtx, 1, 1, size(tlx_data, 3));
tlx_group_hvas = repmat(hva_mtx, 1, 1, size(tlx_data, 3));
tlx_group_genotype = repmat({'tlx'}, size(tlx_data));


% create 2 cols (lat, med)
Y = cat(1, emx_data(:), tlx_data(:));
group_layers = cat(1, emx_group_layers(:), tlx_group_layers(:));
group_hvas = cat(1, emx_group_hvas(:), tlx_group_hvas(:));
group_genotype = cat(1, emx_group_genotype(:), tlx_group_genotype(:));

group = {group_genotype, group_hvas, group_layers};


[P, T, STATS, TERMS] = anovan(Y, group, 'varnames', {'genotype', 'hva', 'layers'});
figure
multcompare(STATS, 'dimension', [2,3])









