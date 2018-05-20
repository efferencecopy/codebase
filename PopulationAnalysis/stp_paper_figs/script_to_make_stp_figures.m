% Load the data sets (this takes a while)

fin

DATASET = 'wcstp';

switch DATASET
    case 'wcstp'
        load('\\crash.dhe.duke.edu\charlie\pre_processed_data\wcstp_180219_with_ddff_fits_and_grand_fit.mat');
        dat_wcstp = dat;
        clear dat
    case 'passive'
        load('\\crash.dhe.duke.edu\charlie\pre_processed_data\passive_props_180219.mat');
        dat_passive = dat;
        clear dat
    case 'interneurons'
        load('\\crash.dhe.duke.edu\charlie\pre_processed_data\interneurons_180219.mat');
        dat_interneurons = dat;
        clear dat
    case 'multipower'
        load('\\crash.dhe.duke.edu\charlie\pre_processed_data\interneurons_multipower_180219.mat');
        dat_multipower = dat;
        clear dat
end

% CD to figure dir for easier saving
cd('\\crash.dhe.duke.edu\charlie\wcstp_manuscript_stuff\wcstp_figures_second_pass')


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

%% Figure 2: DC injections

close all; clc;

% example cell: CH_180118_B site2 ch2 AM
example_idx = 221;
example_ch = 2;

fig_dc_injections(dat_passive, example_idx, example_ch)


%% Figure 2: Input resistance

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

% example cells:
ex_cells = {
'CH_170220_B', 1, 2; ...
'CH_170829_A', 2, 1; ...
'CH_170829_A', 2, 2; ...
'CH_170829_C', 3, 2; ...  % keep
'CH_180104_B', 1, 1; ...
'CH_180118_A', 3, 1; ...
'CH_180118_C', 1, 2; ...  % keeps
};


for i_cell = 1:size(ex_cells,1)
    fig_example_stp_single_cell(dat_wcstp, ex_cells{i_cell, 1}, ex_cells{i_cell, 2}, ex_cells{i_cell, 3});
end


%% Figure 3: Pn:P1 for all P (PY cells)

close all; clc

ppr_groups = plotgroups_wcstp_all;

options.PLOT_AVG_MANIFOLD = true;
options.FORCE_PAIRED_RECORDINGS = false;
options.LOG_SPACE = false;

[recovpop, groupdata] = fig_pnp1_ratios(dat_wcstp, ppr_groups, pprpop, options);
quant_test_pprs(recovpop, groupdata, ppr_groups)



%% Figure 7?: Pn:P1 for all P (Inter neurons)

close all; clc

plotgroups_wcstp_ins = {
    'all_som', 'L23', 'med', 'any';...
    'all_som', 'L23', 'lat', 'any';...
    'all_pv', 'L23', 'med', 'any';...
    'all_pv', 'L23', 'lat', 'any';...
    'PY', 'L23', 'med', 'any';...
    'PY', 'L23', 'lat', 'any';...
};

options.PLOT_AVG_MANIFOLD = true;
options.FORCE_PAIRED_RECORDINGS = true;
options.LOG_SPACE = true;

[recovpop, groupdata] = fig_pnp1_ratios(dat_wcstp, plotgroups_wcstp_ins, pprpop, options);
quant_test_pprs(recovpop, groupdata, plotgroups_wcstp_ins)


%% Figure 3: P10 vs P3 scatter plot

close all; clc

fig_ppr_scatter(dat_wcstp, plotgroups_wcstp_ml)



