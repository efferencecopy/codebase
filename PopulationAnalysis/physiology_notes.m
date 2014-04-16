%% TEMPLATE VERSION (mouse name and cell num)

%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = '';      % The mouse's name
params.cellNum = nan;    % The neuron number that day
params.DCsteps = '';    % DC steps for Rin and cell identification
params.photo = '';      % To assess where the light stimulus was, and the HOA that contains each cell
params.validCh = '';    % 'HS2_' or 'HS1_'
params.files = {};      % File names of the raw data. <Nx1> cell array
params.skipSweeps = {}; % In case I need to ignore certain sweeps
params.vHold = nan;     % The holding potential for vClamp experiments. One value for each expt. <Nx1>
params.stimLoc = [];    % The (x,y) coordinates of the objective at the locations stimulated with the LED. <Nx2> matrix
params.legTxt = {};     % Text that will appear in figures to annotate each data file. <Nx1> cell array
params.tags = {};


%
% ANALYZE OR ADD TO MDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('SUPPRESS_ANALYSIS', 'var')
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    invitroAnalysisOverview(params)
end

if exist('ADD_TO_MDB', 'var')
    [~, idx] = mdb.search(params.mouse);
    mdb.mice{idx}.popAnly{params.cellNum} = params;
end



%% CH_020314_B Cell 1


%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% I think that this cell is located in V1. The injection site was very
% anterior, and the center of the injection is at about the same A/P
% positin as the thalamic expression. So, this cell is either in V1, or in
% a HOA just anterior to V1. Too close to call

% Pop Analysis:
% Not added to any population analysis.


%% CH_020314_B Cell 2
fin

%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I think that this cell is located in LM. The slice anterior to this has a
% HOA that is displaced medially (AL?) and this slice has only weak
% thalamic signal. Series resistance is pretty good, but jumps between the
% 3rd and 4th files.
%
% Some facilitation at 40 and 20 Hz, depression at the other freqs. This is
% a Layer 2/3 PY cell. Funky doublets with big current injections.


%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'CH_020314_B';      % The mouse's name
params.cellNum = 2;    % The neuron number that day
params.DCsteps = '2014_02_20_0021';    % DC steps for Rin and cell identification
params.photo = 'cell2_5x_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.validCh = 'HS2_';    % 'HS2_' or 'HS1_'
params.files = {'2014_02_20_0025';...
                '2014_02_20_0026';...
                '2014_02_20_0027';...
                '2014_02_20_0032';...
                '2014_02_20_0033'};      % File names of the raw data. <Nx1> cell array
params.skipSweeps = {}; % In case I need to ignore certain sweeps
params.vHold = -85.* ones(numel(params.files), 1);     % The holding potential for vClamp experiments. One value for each expt. <Nx1>
params.stimLoc = [];    % The (x,y) coordinates of the objective at the locations stimulated with the LED. <Nx2> matrix
params.legTxt = {'5', '40', '60', '2', '20'};     % Text that will appear in figures to annotate each data file. <Nx1> cell array


%
% ANALYZE OR ADD TO MDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('SUPPRESS_ANALYSIS', 'var')
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    invitroAnalysisOverview(params)
end

if exist('ADD_TO_MDB', 'var')
    [~, idx] = mdb.search(params.mouse);
    mdb.mice{idx}.popAnly{params.cellNum} = params;
end


%% CH_020314_B Cell 3
fin

%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I think that this cell is located in LM. The slice anterior to this has a
% HOA that is displaced medially (AL?) and this slice has only weak
% thalamic signal. Series resistance is pretty good, but jumps between the
% 3rd and 4th files.


%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'CH_020314_B';      % The mouse's name
params.cellNum = 3;    % The neuron number that day
params.DCsteps = '2014_02_20_0040';    % DC steps for Rin and cell identification
params.photo = 'cell3_5x_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.validCh = 'HS2_';    % 'HS2_' or 'HS1_'
params.files = {'2014_02_20_0049';...
                '2014_02_20_0052';...
                '2014_02_20_0057'};      % File names of the raw data. <Nx1> cell array
params.skipSweeps = {}; % In case I need to ignore certain sweeps
params.vHold = -85.* ones(numel(params.files), 1);     % The holding potential for vClamp experiments. One value for each expt. <Nx1>
params.stimLoc = [];    % The (x,y) coordinates of the objective at the locations stimulated with the LED. <Nx2> matrix
params.legTxt = {'20', '40', '5'};     % Text that will appear in figures to annotate each data file. <Nx1> cell array


%
% ANALYZE OR ADD TO MDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('SUPPRESS_ANALYSIS', 'var')
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    invitroAnalysisOverview(params)
end

if exist('ADD_TO_MDB', 'var')
    [~, idx] = mdb.search(params.mouse);
    mdb.mice{idx}.popAnly{params.cellNum} = params;
end



%% CH_020314_C Cells 1-4

%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This mouse had an injection that was quite anterior, but it's possible
% that LM and AL were well separated. Nonetheless, these cells were
% recorded in unknown areas (some quite anterior). The data files are a
% mixture of holding potentials (some at 0mV) and frequencies. All
% stimulation at the soma. No evoked responses at the axon bundle.

% Pop Analysis: none of the cells added to pop anly.



%% CH_020314_D Cell 2

%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I think that this cell is located in either posterior AL or anterior LM.
% It seems like there is a slight transition between the lateral location
% of the axon fields for this cell and the field in the adjcent slice.
% Also, the slice posterior to this one has no thalamic expression, but
% this slice does have some signal in the thalamus. 

% Type: L2/3 PY cell. There's some facilitation of 20 and 40 Hz stimuli,
% but depression for 5 Hz stimuli.

%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'CH_020314_D';      % The mouse's name
params.cellNum = 2;    % The neuron number that day
params.DCsteps = '2014_02_22_0015';    % DC steps for Rin and cell identification
params.photo = 'cell2_5x_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.validCh = 'HS2_';    % 'HS2_' or 'HS1_'
params.files = {'2014_02_22_0022';...
                '2014_02_22_0024';...
                '2014_02_22_0027'};      % File names of the raw data. <Nx1> cell array
params.skipSweeps = {}; % In case I need to ignore certain sweeps
params.vHold = -85 .* ones(numel(params.files), 1);     % The holding potential for vClamp experiments. One value for each expt. <Nx1>
params.stimLoc = [];    % The (x,y) coordinates of the objective at the locations stimulated with the LED. <Nx2> matrix
params.legTxt = {'20', '5', '40'};     % Text that will appear in figures to annotate each data file. <Nx1> cell array


%
% ANALYZE OR ADD TO MDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('SUPPRESS_ANALYSIS', 'var')
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    invitroAnalysisOverview(params)
end

if exist('ADD_TO_MDB', 'var')
    [~, idx] = mdb.search(params.mouse);
    mdb.mice{idx}.popAnly{params.cellNum} = params;
end


%% CH_032414_A Cell 3

fin

%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I didn't save any images from the physiology rig, and there wasn't any
% histology, but this cell has a small facilitation for 20/40 Hz and nearly
% unity P1:P2 ratio. 



%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'CH_032414_A';      % The mouse's name
params.cellNum = 3;    % The neuron number that day
params.DCsteps = '2014_04_06_0009';    % DC steps for Rin and cell identification
params.photo = '';      % To assess where the light stimulus was, and the HOA that contains each cell
params.validCh = 'HS2_';    % 'HS2_' or 'HS1_'
params.files = {'2014_04_06_0011', '2014_04_06_0012', '2014_04_06_0013'};      % File names of the raw data. <Nx1> cell array
params.skipSweeps = {}; % In case I need to ignore certain sweeps
params.vHold = [-85 -85 -85];     % The holding potential for vClamp experiments. One value for each expt. <Nx1>
params.stimLoc = [];    % The (x,y) coordinates of the objective at the locations stimulated with the LED. <Nx2> matrix
params.legTxt = {'20', '5', '40'};     % Text that will appear in figures to annotate each data file. <Nx1> cell array
params.tags = {};


%
% ANALYZE OR ADD TO MDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('SUPPRESS_ANALYSIS', 'var')
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    invitroAnalysisOverview(params)
end

if exist('ADD_TO_MDB', 'var')
    [~, idx] = mdb.search(params.mouse);
    mdb.mice{idx}.popAnly{params.cellNum} = params;
end



%% CH_032414_B Cell 1

%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Weak evoked responses from a neuron in an unknown HOA. No picture saved,
% no histology. All stimuli delivered to the soma.


%% CH_032414_B Cell 2

%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% possibly from area AL, but no picutures or histology were saved. All
% stimuli delivered to the soma. 



%% CH_032414_B Cell 3

fin

%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I don't know where this cell is from, but possibly LM/AL. Read the notes
% in the docubase from the previous cell (Cell 2). This is an interneuron.
% There is facilitation when the light targets the soma, but depression
% when the light targets an interior location. TTX was washed on at the
% end. No responses in TTX.


%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'CH_032414_B';      % The mouse's name
params.cellNum = 3;    % The neuron number that day
params.DCsteps = '2014_04_07_0020';    % DC steps for Rin and cell identification
params.photo = '';      % To assess where the light stimulus was, and the HOA that contains each cell
params.validCh = 'HS2_';    % 'HS2_' or 'HS1_'
params.files = {'2014_04_07_0024';...
                '2014_04_07_0025';...
                '2014_04_07_0026';...
                '2014_04_07_0028';...
                '2014_04_07_0029';...
                '2014_04_07_0030'};      % File names of the raw data. <Nx1> cell array
params.skipSweeps = {}; % In case I need to ignore certain sweeps
params.vHold = -85.*ones(numel(params.files),1);     % The holding potential for vClamp experiments. One value for each expt. <Nx1>
params.stimLoc = [0, 0;...
                  0, 0;...
                  0, 0;...
                  -151 178;...
                  -151 178;...
                  -151 178];    % The (x,y) coordinates of the objective at the locations stimulated with the LED. <Nx2> matrix
params.legTxt = {'40 soma', '20 soma', '5 soma', '40 dist', '20 dist', '5 dist'};     % Text that will appear in figures to annotate each data file. <Nx1> cell array
params.tags = {'interneuron TTX'};


%
% ANALYZE OR ADD TO MDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('SUPPRESS_ANALYSIS', 'var')
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    invitroAnalysisOverview(params)
end

if exist('ADD_TO_MDB', 'var')
    [~, idx] = mdb.search(params.mouse);
    mdb.mice{idx}.popAnly{params.cellNum} = params;
end



%% CH_032414_C Cell 3

%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A layer 2/3 PY cell. Possibly from area LM. Facilitation at the soma for
% 20 Hz stimuli but depression when stimulating axons further away... Check
% that the pulse widths are the same across stimulating locations.
%
% Added this cell to a pop anly about TF dependence of PPR for somatic and
% axonal stimulation.


%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'CH_032414_C';      % The mouse's name
params.cellNum = 3;    % The neuron number that day
params.DCsteps = '2014_04_09_0015';    % DC steps for Rin and cell identification
params.photo = 'cell_3_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.validCh = 'HS2_';    % 'HS2_' or 'HS1_'
params.files = {'2014_04_09_0018';...
                '2014_04_09_0025';...
                '2014_04_09_0026';...
                '2014_04_09_0027'};  % File names of the raw data. <Nx1> cell array
params.skipSweeps = {[], [], [], [23:25]}; % In case I need to ignore certain sweeps
params.vHold = -85*ones(size(params.files));     % The holding potential for vClamp experiments. One value for each expt. <Nx1>
params.stimLoc = [0 0;...
                  -12 234;...
                  -12 234;...
                  -12 234];    % The (x,y) coordinates of the objective at the locations stimulated with the LED. <Nx2> matrix
params.legTxt = {'soma 20', 'dist 20', 'dist 5', 'dist 40'};     % Text that will appear in figures to annotate each data file. <Nx1> cell array
params.tags = {};


%
% ANALYZE OR ADD TO MDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('SUPPRESS_ANALYSIS', 'var')
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    invitroAnalysisOverview(params)
end

if exist('ADD_TO_MDB', 'var')
    [~, idx] = mdb.search(params.mouse);
    mdb.mice{idx}.popAnly{params.cellNum} = params;
end




%% CH_032414_C Cell 4

%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = '';      % The mouse's name
params.cellNum = nan;    % The neuron number that day
params.DCsteps = '';    % DC steps for Rin and cell identification
params.photo = '';      % To assess where the light stimulus was, and the HOA that contains each cell
params.validCh = '';    % 'HS2_' or 'HS1_'
params.files = {};      % File names of the raw data. <Nx1> cell array
params.skipSweeps = {}; % In case I need to ignore certain sweeps
params.vHold = nan;     % The holding potential for vClamp experiments. One value for each expt. <Nx1>
params.stimLoc = [];    % The (x,y) coordinates of the objective at the locations stimulated with the LED. <Nx2> matrix
params.legTxt = {};     % Text that will appear in figures to annotate each data file. <Nx1> cell array
params.tags = {};


%
% ANALYZE OR ADD TO MDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('SUPPRESS_ANALYSIS', 'var')
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    invitroAnalysisOverview(params)
end

if exist('ADD_TO_MDB', 'var')
    [~, idx] = mdb.search(params.mouse);
    mdb.mice{idx}.popAnly{params.cellNum} = params;
end




%% EB_031014_D Cell 2

fin

%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'EB_031014_D';      % The mouse's name
params.cellNum = 2;    % The neuron number that day
params.DCsteps = '2014_04_01_0012';    % DC steps for Rin and cell identification
params.photo = 'cell_2_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.validCh = 'HS2_';    % 'HS2_' or 'HS1_'
params.files = {'2014_04_01_0017';...
               '2014_04_01_0020';...
               '2014_04_01_0029';...
               '2014_04_01_0030';...
               '2014_04_01_0032';...
               '2014_04_01_0035';...
               '2014_04_01_0039'};      % File names of the raw data. <Nx1> cell array
params.skipSweeps = {}; % In case I need to ignore certain sweeps
params.vHold = -85 .* ones(numel(params.files), 1);     % The holding potential for vClamp experiments. One value for each expt. <Nx1>
params.stimLoc = [0 0;...
                 -81 290;...
                 -183 451;...
                 -183 451;...
                 -183 451;...
                 -473 437;...
                 0 0];    % The (x,y) coordinates of the objective at the locations stimulated with the LED. <Nx2> matrix
params.legTxt = {'Soma FS Open';...
                 'Cortex 1';...
                 'Cortex 2 4.5 volts';...
                 'Cortex 2 3 volts';...
                 'Cortex 2 10 volts';...
                 'Cortex 3 ';...
                 'Soma repeat FS open'};     % Text that will appear in figures to annotate each data file. <Nx1> cell array


%
% ANALYZE OR ADD TO MDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('SUPPRESS_ANALYSIS', 'var')
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    invitroAnalysisOverview(params)
end

if exist('ADD_TO_MDB', 'var')
    [~, idx] = mdb.search(params.mouse);
    mdb.mice{idx}.popAnly{params.cellNum} = params;
end



%% EB_031014_B Cell 2
 
fin

%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'EB_031014_B';      % The mouse's name
params.cellNum = 2;    % The neuron number that day
params.DCsteps = '2014_03_28_0023';    % DC steps for Rin and cell identification
params.photo = 'cell_2_tdTomato_5x_470nm';      % To assess where the light stimulus was, and the HOA that contains each cell
params.validCh = 'HS1_';    % 'HS2_' or 'HS1_'
params.files = {'2014_03_28_0026';...
                '2014_03_28_0027';...
                '2014_03_28_0029';...
                '2014_03_28_0030';...
                '2014_03_28_0032';...
                '2014_03_28_0036';...
                '2014_03_28_0037';...
                '2014_03_28_0039'};      % File names of the raw data. <Nx1> cell array
params.skipSweeps = {}; % In case I need to ignore certain sweeps
params.vHold =  -85 .* ones(numel(params.files), 1);     % The holding potential for vClamp experiments. One value for each expt. <Nx1>
params.stimLoc = [0 0;...
                  0 0;...
                  215 -51;...
                  176 117;...
                  299 409;...
                  -70 291;...
                 -70 291;...
                 -70 291];    % The (x,y) coordinates of the objective at the locations stimulated with the LED. <Nx2> matrix
params.legTxt = { 'soma FS open';...
                  'soma FS closed';...
                  'Lateral site 1';...
                  'Lateral site 2';...
                  'Lateral site 3';...
                  'Lateral site 4';...
                  'HalfMoon Incorrect';...
                  'HalfMoon Correct'};     % Text that will appear in figures to annotate each data file. <Nx1> cell array


%
% ANALYZE OR ADD TO MDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('SUPPRESS_ANALYSIS', 'var')
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    invitroAnalysisOverview(params)
end

if exist('ADD_TO_MDB', 'var')
    [~, idx] = mdb.search(params.mouse);
    mdb.mice{idx}.popAnly{params.cellNum} = params;
end


%% EB_031014_B Cell 1

%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%
% Params
%%%%%%%%%%%%%%%%%%%%%%%%%
params.mouse = 'EB_031014_B';      % The mouse's name
params.cellNum = 1; 
params.DCsteps = '2014_03_28_0000';    % DC steps for Rin and cell identification
params.photo = 'cell_1_tdTomato_5x';      % To assess where the light stimulus was, and the HOA that contains each cell
params.validCh = 'HS2_';    % 'HS2_' or 'HS1_'
params.files = {'2014_03_28_0009';...
                '2014_03_28_0011',;...
                '2014_03_28_0013';...
                '2014_03_28_0016';...
                '2014_03_28_0019';...
                '2014_03_28_0020';...
                '2014_03_28_0021';...
                '2014_03_28_0022'};      % File names of the raw data. <Nx1> cell array
params.skipSweeps = {}; % In case I need to ignore certain sweeps
params.vHold = -85 .* ones(numel(params.files), 1);     % The holding potential for vClamp experiments. One value for each expt. <Nx1>
params.stimLoc = [-179 -48;...
                  -179 -48;...
                  -288 -118;...
                  0 0;...
                  -178 -169;...
                  -178 -169;...
                  -153 -147;...
                  -113 -111];    % The (x,y) coordinates of the objective at the locations stimulated with the LED. <Nx2> matrix
params.legTxt = {'distal pos 1 FS closed';...
                 'distal pos 1 FS open';...
                 'distal pos 2';...
                 'Soma FS closed';...
                 'distal pos 3';...
                 'distal pos 3 half moon';...
                 'distal pos 4 half moon';...
                 'distal pos 5 half moon'};     % Text that will appear in figures to annotate each data file. <Nx1> cell array

%
% ANALYZE OR ADD TO MDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('SUPPRESS_ANALYSIS', 'var')
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    invitroAnalysisOverview(params)
end

if exist('ADD_TO_MDB', 'var')
    [~, idx] = mdb.search(params.mouse);
    mdb.mice{idx}.popAnly{params.cellNum} = params;
end

