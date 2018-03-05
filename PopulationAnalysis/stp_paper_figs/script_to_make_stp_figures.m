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


%% Figure 1: DC injections

close all; clc;

% example cell: CH_180118_B site2 ch2 AM
example_idx = 221;
example_ch = 2;

fig_1_dc_injections(dat_passive, example_idx, example_ch)


%% Figure 1: Input resistance

close all; clc

% define a set of attributes for each analysis group
% {CellType, Layer,  BrainArea,  OpsinType}
% Brain Area can be: 'AL', 'PM', 'AM', 'LM', 'any', 'med', 'lat'. CASE SENSITIVE
plotgroups = {
    'PY', 'L23', 'PM', 'any';...
    'PY', 'L23', 'LM', 'any';...
    'PY', 'L23', 'AM', 'any';...
    'PY', 'L23', 'AL', 'any';...
    };

fig_1_Rin(dat_passive, plotgroups);


%% Figure 1: Membrane time constant

close all; clc

% define a set of attributes for each analysis group
% {CellType, Layer,  BrainArea,  OpsinType}
% Brain Area can be: 'AL', 'PM', 'AM', 'LM', 'any', 'med', 'lat'. CASE SENSITIVE
plotgroups = {
    'PY', 'L23', 'PM', 'any';...
    'PY', 'L23', 'LM', 'any';...
    'PY', 'L23', 'AM', 'any';...
    'PY', 'L23', 'AL', 'any';...
    };

fig_1_tau(dat_passive, plotgroups);








