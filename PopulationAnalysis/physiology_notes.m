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


%
% ANALYZE OR ADD TO MDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('SUPPRESS_ANALYSIS', 'var')
    params.fxns = {@pulseTrains};
    invitroAnalysisOverview(params)
end

if exist('ADD_TO_MDB', 'var')
    [~, idx] = mdb.search(params.mouse);
    mdb.mice{idx}.popAnly{params.cellNum} = params;
end



%% EB_031014_B Cell 2

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
params.vHold = -85;     % The holding potential for vClamp experiments. One value for each expt. <Nx1>
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
    params.fxns = {@pulseTrains};
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
    params.fxns = {@pulseTrains};
    invitroAnalysisOverview(params)
end

if exist('ADD_TO_MDB', 'var')
    [~, idx] = mdb.search(params.mouse);
    mdb.mice{idx}.popAnly{params.cellNum} = params;
end

