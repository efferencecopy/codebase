% Load the data sets (this takes a while)

fin

DATASET = 'interneurons';

file_path = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\charlie\pre_processed_data\';
switch DATASET
    case 'wcstp'
        load(strcat(file_path, 'wcstp_180219_with_ddff_fits_and_grand_fit.mat'));
        dat_wcstp = dat;
        clear dat
    case 'passive'
        load(strcat(file_path, 'passive_props_180219.mat'));
        dat_passive = dat;
        clear dat
    case 'interneurons'
        load(strcat(file_path, 'interneurons_180219.mat'));
        dat_interneurons = dat;
        clear dat
    case 'multipower'
        load(strcat(file_path, 'interneurons_multipower_180219.mat'));
        dat_multipower = dat;
        clear dat
end

% CD to figure dir for easier saving
cd('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\charlie\wcstp_manuscript_stuff\wcstp_figures_thrid_pass')


% define a set of attributes for each analysis group
% {CellType, Layer,  BrainArea,  OpsinType}
% Brain Area can be: 'AL', 'PM', 'AM', 'LM', 'any', 'med', 'lat'. CASE SENSITIVE
plotgroups_passive_all = {
    'PY', 'L23', 'PM', 'any';...
    'PY', 'L23', 'AM', 'any';...
    'PY', 'L23', 'LM', 'any';...
    'PY', 'L23', 'AL', 'any';...
    };

plotgroups_passive_ml = {
    'PY', 'L23', 'med', 'any';...
    'PY', 'L23', 'lat', 'any';...
    };

plotgroups_wcstp_ml = {
    'PY', 'L23', 'med', 'chief';...
    'PY', 'L23', 'lat', 'chief';...
    };

plotgroups_wcstp_all = {
    'PY', 'L23', 'PM', 'chief';...
    'PY', 'L23', 'LM', 'chief';...
    'PY', 'L23', 'AM', 'chief';...
    'PY', 'L23', 'AL', 'chief';...
    };


%% Figure 2: DC injections example cell

close all; clc;


% % example cell for PY: CH_180118_B site2 ch2 AM
example_idx = 221;
example_ch = 2;

% example cell for SOM
% example_idx = 39;
% example_ch = 1;

% example Cell for PV 
% example_idx = 54;
% example_ch = 1;

fig_dc_injections(dat_passive, example_idx, example_ch)


%% Figure 2: Input resistance (pyramidal cells)

close all; clc

fig_Rin(dat_passive, plotgroups_passive_all);


%% Figure 2: Membrane time constant

close all; clc

fig_tau(dat_passive, plotgroups_passive_all);


%% Figure 2: Input - output functions

close all; clc

fig_input_output_curves(dat_passive, plotgroups_passive_all)

%% Figure3: Example PY cell EPSCs to Trains

close all; clc

% example PY cells for STP:
ex_cells = {
'CH_170220_B', 1, 2; ...
'CH_170829_A', 2, 1; ...
'CH_170829_A', 2, 2; ...
'CH_170829_C', 3, 2; ...  % keep
'CH_180104_B', 1, 1; ...
'CH_180118_A', 3, 1; ...
'CH_180118_C', 1, 2; ...  % keeps
};

% % example cells for PTP from lateral HVAs
% ex_cells = {
% 'CH_180104_C', 1, 2; ...
% 'CH_170112_C', 4, 2; ...
% 'CH_161022_B', 2, 2; ...
% 'CH_160929_B', 1, 1; ...
% 'CH_160915_B', 2, 2; ...
% };

% % example cells for PTP from medial HVAs
% ex_cells = {
% 'CH_160915_B', 1, 2; ...
% 'CH_161026_A', 2, 1; ...
% };

for i_cell = 1:size(ex_cells,1)
    fig_example_stp_single_cell(dat_wcstp, ex_cells{i_cell, 1}, ex_cells{i_cell, 2}, ex_cells{i_cell, 3});
end

%% Figure 3: Pn:P1 for all P (PY cells)

close all; clc

ppr_groups = plotgroups_wcstp_all;

options.PLOT_AVG_MANIFOLD = false;
options.FORCE_PAIRED_RECORDINGS = false;
options.LOG_SPACE = false;

[recovpop, groupdata] = fig_pnp1_ratios(dat_wcstp, ppr_groups, pprpop, options);
quant_test_pprs(recovpop, groupdata, ppr_groups)


%% Figure 6: Input resistance (interneurons)

close all; clc

plotgroups_passive_ins = {
    'all_som', 'any', 'any', 'any';...
    'all_pv', 'any', 'any', 'any';...
    };

fig_Rin(dat_interneurons, plotgroups_passive_ins);



%% Figure 7?: Pn:P1 for all P (Inter neurons)

close all; clc

plotgroups_wcstp_ins = {
    'all_pv', 'L23', 'med', 'any';...
    'all_pv', 'L23', 'lat', 'any';...
};

options.PLOT_AVG_MANIFOLD = false;
options.FORCE_PAIRED_RECORDINGS = false;
options.LOG_SPACE = true;

[recovpop, groupdata] = fig_pnp1_ratios(dat_wcstp, plotgroups_wcstp_ins, pprpop, options);
quant_test_pprs(recovpop, groupdata, plotgroups_wcstp_ins)


%% Figure 3: P10 vs P3 scatter plot

close all; clc

fig_ppr_scatter(dat_wcstp, plotgroups_wcstp_ml)




