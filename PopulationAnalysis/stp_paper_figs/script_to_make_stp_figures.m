% Load the data sets (this takes a while)

fin

DATASET = 'passive';

switch DATASET
    case 'wcstp'
        load('\\crash.dhe.duke.edu\charlie\pre_processed_data\wcstp_180219.mat');
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

% define a set of attributes for each analysis group
% {CellType, Layer,  BrainArea,  OpsinType}
% Brain Area can be: 'AL', 'PM', 'AM', 'LM', 'any', 'med', 'lat'. CASE SENSITIVE
plotgroups_all = {
    'PY', 'L23', 'PM', 'any';...
    'PY', 'L23', 'LM', 'any';...
    'PY', 'L23', 'AM', 'any';...
    'PY', 'L23', 'AL', 'any';...
    };

plotgroups_ml = {
    'PY', 'L23', 'med', 'any';...
    'PY', 'L23', 'lat', 'any';...
    };


%% Figure 2: DC injections

close all; clc;

% example cell: CH_180118_B site2 ch2 AM
example_idx = 221;
example_ch = 2;

fig_dc_injections(dat_passive, example_idx, example_ch)


%% Figure 2: Input resistance

close all; clc

fig_Rin(dat_passive, plotgroups_ml);


%% Figure 2: Membrane time constant

close all; clc

fig_tau(dat_passive, plotgroups_all);


%% Figure 2: Input - output functions

close all; clc

fig_input_output_curves(dat_passive, plotgroups_all)







