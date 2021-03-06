%% TEMPLATE VERSION (mouse name and cell num)

fin

%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = '';      % The mouse's name
params.cellNum = nan;   % The neuron number that day
params.DCsteps = '';    % DC steps for Rin and cell identification
params.photo = '';      % To assess where the light stimulus was, and the HOA that contains each cell
params.validCh = 'HS2_';% 'HS2_' or 'HS1_'
params.files = {};      % File names of the raw data. <Nx1> cell array
params.skipSweeps = {}; % In case I need to ignore certain sweeps
params.vHold = nan;     % The holding potential for vClamp experiments. One value for each expt. <Nx1>
params.stimLoc = [];    % The (x,y) coordinates of the objective at the locations stimulated with the LED. <Nx2> matrix
params.legTxt = {};     % Text that will appear in figures to annotate each data file. <Nx1> cell array
params.tags = {};
params.subtractexp = false;
params.filter = 4000;

% 
% %
% % ANALYZE OR ADD TO PARAMSDB
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
%     params.fxns = {@anlyMod_pulseTrains_stimLoc};
%     params = invitroAnalysisOverview(params);
% end
% 
% if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
%     addPopAnlyParamsToMDB(params);
% end



%% AK_072814_A pair 1

fin
%
% Notes
%
%%%%%%%%%%%%%%%%%
% Ok data from CH2, less than great data from CH1. I'm not certain that
% cell 2 is an IN neuron, so this should be verified.
%
% brain area: possibly LM
% popAnalysis: E/I and A/N ratios
%

%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'AK_072814_A';      % The mouse's name
params.cellNum = 1;    % The neuron number that day
params.photo = 'AK_072814_A_pair1_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.files = {'2014_08_12_', [2:10]};  % <file name prefix, suffix>
params.groups = {'control', [2:5];...
                 'nbqxGabazine', [6:10];...
                 'NMDAR', [6:10]};
params.excludeHS1 = {{'_0003', [1:7]}, {'_0002', [1:13, 15:24, 26]}};
params.excludeHS2 = {{'_0003', [2,3,5:8,11,13]}, {'_0002', [1:10, 13:17, 21:26]}};

% stuff for E/I and AMPA/NMDA ratios
% key for isolatedCurrents = {<current><group><Vhold><Erev>}
% Erev is to calculate driving force for conversion from pA to pS
params.isolatedCurrents = {'excit', 'control', -75, 20;...
                           'inhib', 'control', 15, -75;...
                           'ampa', 'control', -75, 20;...
                           'nmda', 'nbqxGabazine', 50, 20};
params.tags = {};
params.filter = 8e3;
params.celldepth = [norm([-102+35 237+8]), norm([102 229])];


%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_optoIV, @anlyMod_EIbalance, @anlyMod_NMDAR};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end



%% AK_072814_B pair 1

fin
%
% Notes
%
%%%%%%%%%%%%%%%%%
% Nice data from 2 PY cells that are likely in PM. Both have good access
% resisitance. 
%
% brain area: likely PM
% popAnalysis: adding to E/I A/N
%

%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'AK_072814_B';      % The mouse's name
params.cellNum = 1;    % The neuron number that day
params.photo = 'AK_072814_B_pair1_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.files = {'2014_08_15_', [2:6]};  % <file name prefix, suffix>
params.groups = {'control', [2:3];...
                 'nbqxGabazine', [4:6];...
                 'NMDAR', [4:6]};
params.excludeHS1 = {{'_0002', [1,9,10,14]}};
params.excludeHS2 = {{'_0002', [1,9,10,14]}};

% stuff for E/I and AMPA/NMDA ratios
% key for isolatedCurrents = {<current><group><Vhold><Erev>}
% Erev is to calculate driving force for conversion from pA to pS
params.isolatedCurrents = {'excit', 'control', -75, 15;...
                           'inhib', 'control', 15, -75;...
                           'ampa', 'control', -75, 15;...
                           'nmda', 'nbqxGabazine', 50, 15};
params.tags = {};
params.filter = 2e3;
HS1loc = [10 28];
HS2loc = [-4 8];
Pialoc = [44 229];
params.celldepth = [norm(HS1loc-Pialoc), norm([HS2loc-Pialoc])];


%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_optoIV, @anlyMod_EIbalance, @anlyMod_NMDAR};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end



%% AK_072814_B cell 3

fin
%
% Notes
%
%%%%%%%%%%%%%%%%%
% Fairly clean data from only a single channel. There is signigicant
% holding current required at +50 and -75 mV after the drugs have washed
% on. The resulting Vclamp error causes the NMDA conductance to be
% underestamated... I should institue a check of this...
%
% brain area: likely LM
% popAnalysis: adding to E/I and A/N. A/N should be filtered out due to
% Vclamp errors.
%

%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'AK_072814_B';      % The mouse's name
params.cellNum = 3;    % The neuron number that day
params.photo = 'AK_072814_B_cell3_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.files = {'2014_08_15_', [12:16]};  % <file name prefix, suffix>
params.groups = {'control', [12:13];...
                 'nbqxGabazine', [14:16]};
params.excludeHS1 = {};
params.excludeHS2 = {};

% stuff for E/I and AMPA/NMDA ratios
% key for isolatedCurrents = {<current><group><Vhold><Erev>}
% Erev is to calculate driving force for conversion from pA to pS
params.isolatedCurrents = {'excit', 'control', -75, 15;...
                           'inhib', 'control', 15, -75;...
                           'ampa', 'control', -75, 15;...
                           'nmda', 'nbqxGabazine', 50, 15};
params.tags = {};
params.filter = 2e3;

HS1loc = [70 2.8];
HS2loc = [nan];
Pialoc = [-225 261];
params.celldepth = [norm(HS1loc-Pialoc), norm([HS2loc-Pialoc])];


%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_optoIV, @anlyMod_EIbalance};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end



%% AK_072814_C pair 1
fin

%
% Notes
%
%%%%%%%%%%%%%%%%%
% Nice data from CH1, junky data from CH2 (the access resistance goes
% downhill fast).
%
% brain area: likely pm, but need to check once histology is done.
% popAnalysis: Adding Cell 1 to E/I and A/N pop list.
%

%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'AK_072814_C';      % The mouse's name
params.cellNum = 1;    % The neuron number that day
params.photo = 'AK_072814_C_pair1_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.files = {'2014_08_13_', [1:9]};  % <file name prefix, suffix>
params.groups = {'control', [1:4];...
                 'nbqxGabazine', [5:9];...
                 'NMDAR', [5:9]};
params.excludeHS1 = {{'_0001', [7,8,10,11,13,14,20,25]},{'_0002', [2,5,6,10:12,19,20,25:28]},{'_0003', [15:25]}, {'_0009', [24]}};
params.excludeHS2 = {{'_0001', [7,8,10,11,13,14,20,25]},{'_0002', [2,5,6,10:12,19,20,25:28]},{'_0003', [15:25]}};

% stuff for E/I and AMPA/NMDA ratios
% key for isolatedCurrents = {<current><group><Vhold><Erev>}
% Erev is to calculate driving force for conversion from pA to pS
params.isolatedCurrents = {'excit', 'control', -75, 15;...
                           'inhib', 'control', 15, -75;...
                           'ampa', 'control', -75, 15;...
                           'nmda', 'nbqxGabazine', 50, 15};
params.tags = {};
params.filter = 2e3;

HS1loc = [-1.5 -23.2];
HS2loc = [15 -15];
Pialoc = [64 167];
params.celldepth = [norm(HS1loc-Pialoc), norm([HS2loc-Pialoc])];


%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_optoIV, @anlyMod_EIbalance, @anlyMod_NMDAR};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end


%% AK_072814_C cell 2

fin
%
% Notes
%
%%%%%%%%%%%%%%%%%
% Limited data from a single PY cell? No synaptic blockers, so only useful
% for E/I analysis. For this neuron I stimulated two differet spatial
% locations, but the data are not interpretable due to the fact that Vclamp
% errors are high on the last (distal) location due to massive amount of
% required holding current. The apparent E/I ratio goes from 6.7 when
% stimulating at the soma to 24 at the axons. But I think this could
% partially be due to Vclamp errors or rundown. 
%
%
% brain area: likely LM
% popAnalysis: E/I only for cell 2
%

%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'AK_072814_C';      % The mouse's name
params.cellNum = 2;    % The neuron number that day
params.photo = 'AK_072814_C_cell2_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell

% at the soma
params.files = {'2014_08_13_', [11,13]};  % <file name prefix, suffix>
params.groups = {'control', [11,13]};

% % over axons
% params.files = {'2014_08_13_', [16,17]};  % <file name prefix, suffix>
% params.groups = {'control', [16,17]};

params.excludeHS1 = {};
params.excludeHS2 = {{'_0011', [7,10,16]}, {'_0016', [4:6,9,16,17,20,22,23,25,26,28]}};

% stuff for E/I and AMPA/NMDA ratios
% key for isolatedCurrents = {<current><group><Vhold><Erev>}
% Erev is to calculate driving force for conversion from pA to pS
params.isolatedCurrents = {'excit', 'control', -75, 15;...
                           'inhib', 'control', 15, -75};
params.tags = {};
params.filter = 2e3;


HS1loc = [nan];
HS2loc = [0 55];
Pialoc = [-299 251];
params.celldepth = [norm(HS1loc-Pialoc), norm([HS2loc-Pialoc])];

%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_optoIV, @anlyMod_EIbalance};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end

%% AK_081814_C Cell 2
fin

%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% These data are from a SOM+ cell in the injection site of a ChIEF
% injection. These recordings show very nice TF dependent facilitation.


%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'AK_081814_C';      % The mouse's name
params.cellNum = 2;    % The neuron number that day
params.DCsteps = '2014_09_04_0001';    % DC steps for Rin and cell identification
params.photo = '';      % To assess where the light stimulus was, and the HOA that contains each cell
params.validCh = 'HS2_';    % 'HS2_' or 'HS1_'
params.files = {'2014_09_04_0002';...
                '2014_09_04_0003';...
                '2014_09_04_0004'};      % File names of the raw data. <Nx1> cell array
params.skipSweeps = {}; % In case I need to ignore certain sweeps
params.vHold = -85.* ones(numel(params.files), 1);     % The holding potential for vClamp experiments. One value for each expt. <Nx1>
params.stimLoc = zeros(5,2);    % The (x,y) coordinates of the objective at the locations stimulated with the LED. <Nx2> matrix
params.legTxt = {'10', '20', '40'};     % Text that will appear in figures to annotate each data file. <Nx1> cell array
params.subtractexp = false;
params.filter = 2000;

%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end

%% AK_082514_C Cell 1
fin

%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% These data are from a SOM+ cell in the injection site of a ChR2
% injection. There's basically no facilitation, and probably some
% depression at both TFs.


%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'AK_082514_C';      % The mouse's name
params.cellNum = 1;    % The neuron number that day
params.DCsteps = '';    % DC steps for Rin and cell identification
params.photo = '';      % To assess where the light stimulus was, and the HOA that contains each cell
params.validCh = 'HS1_';    % 'HS2_' or 'HS1_'
params.files = {'2014_09_17_0000';...
                '2014_09_17_0001'};      % File names of the raw data. <Nx1> cell array
params.skipSweeps = {}; % In case I need to ignore certain sweeps
params.vHold = -85.* ones(numel(params.files), 1);     % The holding potential for vClamp experiments. One value for each expt. <Nx1>
params.stimLoc = zeros(5,2);    % The (x,y) coordinates of the objective at the locations stimulated with the LED. <Nx2> matrix
params.legTxt = {'20', '40'};     % Text that will appear in figures to annotate each data file. <Nx1> cell array
params.subtractexp = false;
params.filter = 2000;

%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end


%% AK_082514_C Cell 2
fin

%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% These data are from a SOM+ cell in the injection site of a ChR2
% injection. There's basically no facilitation, and probably some
% depression at both TFs.


%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'AK_082514_C';      % The mouse's name
params.cellNum = 2;    % The neuron number that day
params.DCsteps = '';    % DC steps for Rin and cell identification
params.photo = '';      % To assess where the light stimulus was, and the HOA that contains each cell
params.validCh = 'HS1_';    % 'HS2_' or 'HS1_'
params.files = {'2014_09_17_0002';...
                '2014_09_17_0003'};      % File names of the raw data. <Nx1> cell array
params.skipSweeps = {}; % In case I need to ignore certain sweeps
params.vHold = -85.* ones(numel(params.files), 1);     % The holding potential for vClamp experiments. One value for each expt. <Nx1>
params.stimLoc = zeros(5,2);    % The (x,y) coordinates of the objective at the locations stimulated with the LED. <Nx2> matrix
params.legTxt = {'20', '40'};     % Text that will appear in figures to annotate each data file. <Nx1> cell array
params.subtractexp = false;
params.filter = 2000;

%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end


%% AK_092914_A cell 1

fin

% NOTES
%%%%%%%%%%%%%%%%%%%%%%%%
%
% brain area: AL
% popAnalysis: EIAN, NMDAR
%

%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'AK_092914_A';      % The mouse's name
params.cellNum = 1;    % The neuron number that day
params.photo = 'AK_092914_A_cell1_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.files = {'2014_10_16_', [0:8]};  % <file name prefix, suffix>
params.groups = {'control', [0,1];...
                 'nbqxGabazine', [2,3];...
                 'NMDAR', [2:8]};
params.excludeHS1 = {};
params.excludeHS2 = {{'_0005', [4,8]}, {'_0008', [1:11, 23:27]}};
params.tags = {};
params.filter = 800;

HS1loc = [nan];
HS2loc = [0 0];
Pialoc = [-203 283];
params.celldepth = [norm(HS1loc-Pialoc), norm([HS2loc-Pialoc])];


% stuff for E/I and AMPA/NMDA ratios
% key for isolatedCurrents = {<current><group><Vhold><Erev>}
% Erev is to calculate driving force for conversion from pA to pS
params.isolatedCurrents = {'excit', 'control', -72, 15;...
                           'inhib', 'control', 15, -72;...
                           'ampa', 'control', -72, 15;...
                           'nmda', 'nbqxGabazine', 50, 15};


%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_optoIV,  @anlyMod_EIbalance, @anlyMod_NMDAR};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end



%% AK_092914_A cell 2

fin

% NOTES
%%%%%%%%%%%%%%%%%%%%%%%%
%
% brain area: AL
% popAnalysis: EIAN, NMDAR
%

%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'AK_092914_A';      % The mouse's name
params.cellNum = 2;    % The neuron number that day
params.photo = 'AK_092914_A_cell2_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.files = {'2014_10_16_', [9:12]};  % <file name prefix, suffix>
params.groups = {'control', [9,10];...
                 'nbqxGabazine', [11,12]};
params.excludeHS1 = {};
params.excludeHS2 = {{'_0012', [5:11]}};
params.tags = {};
params.filter = 800;


% stuff for E/I and AMPA/NMDA ratios
% key for isolatedCurrents = {<current><group><Vhold><Erev>}
% Erev is to calculate driving force for conversion from pA to pS
params.isolatedCurrents = {'excit', 'control', -72, 15;...
                           'inhib', 'control', 15, -72;...
                           'ampa', 'control', -72, 15;...
                           'nmda', 'nbqxGabazine', 50, 15};

HS1loc = [nan];
HS2loc = [0 0];
Pialoc = [-181 -139];
params.celldepth = [norm(HS1loc-Pialoc), norm([HS2loc-Pialoc])];

                       
                       
%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_optoIV,  @anlyMod_EIbalance, @anlyMod_NMDAR};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end



%% AK_092914_B Cell 1

fin

% NOTES
%%%%%%%%%%%%%%%%%%%%%%%%
%
% brain area: PM
% popAnalysis: EAIN, NMDAR
%


%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'AK_092914_B';      % The mouse's name
params.cellNum = 1;    % The neuron number that day
params.photo = 'AK_092914_B_cell1_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.files = {'2014_10_17_', [0:8]};  % <file name prefix, suffix>
params.groups = {'control', [0,1];...
                 'nbqxGabazine', [2,3];...
                 'NMDAR', [2:8]};
params.excludeHS1 = {};
params.excludeHS2 = {{'_0000', 1}, {'_0005', 3:14}, {'_0008', 10}};
params.tags = {};
params.filter = 2e3;


% stuff for E/I and AMPA/NMDA ratios
% key for isolatedCurrents = {<current><group><Vhold><Erev>}
% Erev is to calculate driving force for conversion from pA to pS
params.isolatedCurrents = {'excit', 'control', -72, 15;...
                           'inhib', 'control', 15, -72;...
                           'ampa', 'control', -72, 15;...
                           'nmda', 'nbqxGabazine', 50, 15};


HS1loc = [nan];
HS2loc = [0 0];
Pialoc = [-108 160];
params.celldepth = [norm(HS1loc-Pialoc), norm([HS2loc-Pialoc])];

                       
                       
                       
%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_optoIV,  @anlyMod_EIbalance, @anlyMod_NMDAR};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end



%% AK_092914_B Cell 2

fin

% NOTES
%%%%%%%%%%%%%%%%%%%%%%%%
% Pretty good data, but marginal Rs. Only one cell. Seems like a PY cell
% based off the amount of ihnibition at +17 mV. 
%
% brain area: Likely AL, but I need to verify with the histology
% popAnalysis: adding to EIAN and NMDAR
%


%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'AK_092914_B';      % The mouse's name
params.cellNum = 2;    % The neuron number that day
params.photo = 'AK_092914_B_cell2_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.files = {'2014_10_17_', [10:17]};  % <file name prefix, suffix>
params.groups = {'control', [10,11];...
                 'nbqxGabazine', [12,13];...
                 'NMDAR', [12:17]};
params.excludeHS1 = {};
params.excludeHS2 = {};
params.tags = {};
params.filter = 800;



HS1loc = [nan];
HS2loc = [0 0];
Pialoc = [-68 135];
params.celldepth = [norm(HS1loc-Pialoc), norm([HS2loc-Pialoc])];


% stuff for E/I and AMPA/NMDA ratios
% key for isolatedCurrents = {<current><group><Vhold><Erev>}
% Erev is to calculate driving force for conversion from pA to pS
params.isolatedCurrents = {'excit', 'control', -72, 17;...
                           'inhib', 'control', 17, -72;...
                           'ampa', 'control', -72, 17;...
                           'nmda', 'nbqxGabazine', 50, 17};


%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_optoIV,  @anlyMod_EIbalance, @anlyMod_NMDAR};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end



%% AK_092914_C Cell 1

fin

% NOTES
%%%%%%%%%%%%%%%%%%%%%%%%
% Seems like a unknown cell type in L4. Not much else to say. Only a little
% data.
%
% brain area: AL
% popAnalysis: EIAN (EI only)
%


%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'AK_092914_C';      % The mouse's name
params.cellNum = 1;    % The neuron number that day
params.photo = 'AK_092914_C_cell1_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.files = {'2014_10_18_', [0,2]};  % <file name prefix, suffix>
params.groups = {'control', [0,2]};
params.excludeHS1 = {};
params.excludeHS2 = {};
params.tags = {};
params.filter = 800;



HS1loc = [nan];
HS2loc = [0 0];
Pialoc = [-25 -440];
params.celldepth = [norm(HS1loc-Pialoc), norm([HS2loc-Pialoc])];


% stuff for E/I and AMPA/NMDA ratios
% key for isolatedCurrents = {<current><group><Vhold><Erev>}
% Erev is to calculate driving force for conversion from pA to pS
params.isolatedCurrents = {'excit', 'control', -72, 15;...
                           'inhib', 'control', 15, -72};


%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_optoIV,  @anlyMod_EIbalance, @anlyMod_NMDAR};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end



%% AK_092914_C cell 2

fin

% NOTES
%%%%%%%%%%%%%%%%%%%%%%%%
% Somewhat chaotic data. There weren't the correct Vholds and Rs comp
% settings to measure E/I or A/N ratios, but there is consistent enough
% data to include the NMDAR data. NB-488 fill seems to suggest that this
% cell is a PY cell.
%
% brain area: AL
% popAnalysis:  EIAN, NMDAR
%


%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'AK_092914_C';      % The mouse's name
params.cellNum = 2;    % The neuron number that day
params.photo = 'AK_092914_C_cell2';      % To assess where the light stimulus was, and the HOA that contains each cell
params.files = {'2014_10_18_', [5,9,10,12,18:25]};  % <file name prefix, suffix>
params.groups = {'control', [5,9];...
                 'nbqxGabazine', [10,12];...
                 'NMDAR', [18:25]};
params.excludeHS1 = {};
params.excludeHS2 = {{'_0025', [6]}};
params.tags = {};
params.filter = 800;


HS1loc = [nan];
HS2loc = [0 0];
Pialoc = [-43 -225];
params.celldepth = [norm(HS1loc-Pialoc), norm([HS2loc-Pialoc])];


% stuff for E/I and AMPA/NMDA ratios
% key for isolatedCurrents = {<current><group><Vhold><Erev>}
% Erev is to calculate driving force for conversion from pA to pS
params.isolatedCurrents = {'excit', 'control', -72, 15;...
                           'inhib', 'control', 15, -72;...
                           'ampa', 'control', -72, 15;...
                           'nmda', 'nbqxGabazine', 50, 15};


%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_optoIV,  @anlyMod_EIbalance, @anlyMod_NMDAR};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end



%% AK_092914_C Cell 3

fin

% NOTES
%%%%%%%%%%%%%%%%%%%%%%%%
% Pretty clean data. This cell is deep (upper L5 or lower L4). Likely a PY
% cell, although there isn't much spontaneous activity on to this cell at
% any Vhold.
%
%
% brain area: AL? Need to confirm with histology
% popAnalysis: Adding to EIAN and NMDAR. 
%


%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'AK_092914_C';      % The mouse's name
params.cellNum = 3;    % The neuron number that day
params.photo = 'AK_092914_C_cell3_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.files = {'2014_10_18_', [26:36]};  % <file name prefix, suffix>
params.groups = {'control', [28,29];...
                 'nbqxGabazine', [30,31];...
                 'NMDAR', [30:36]};
params.excludeHS1 = {};
params.excludeHS2 = {{'_0035', [13:19]}};
params.tags = {};
params.filter = 800;


HS1loc = [nan];
HS2loc = [0 0];
Pialoc = [328 480];
params.celldepth = [norm(HS1loc-Pialoc), norm([HS2loc-Pialoc])];


% stuff for E/I and AMPA/NMDA ratios
% key for isolatedCurrents = {<current><group><Vhold><Erev>}
% Erev is to calculate driving force for conversion from pA to pS
params.isolatedCurrents = {'excit', 'control', -72, 17;...
                           'inhib', 'control', 17, -72;...
                           'ampa', 'control', -72, 17;...
                           'nmda', 'nbqxGabazine', 50, 17};


%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_optoIV,  @anlyMod_EIbalance, @anlyMod_NMDAR};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end







%% AK_092914_D Cell 4

fin

% NOTES
%%%%%%%%%%%%%%%%%%%%%%%%
% this cell is barely in L2/3 and arguably in L4. I should sanity check
% this with someone. Otherwise clean data, although there may be some Ca2+
% current at Vhold -20mV in the presence of synaptic blockers.
%
% brain area: Possibly AL, but I'll need to check
% popAnalysis: Adding to EIAN and NMDAR
%


%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'AK_092914_D';      % The mouse's name
params.cellNum = 4;    % The neuron number that day
params.photo = 'AK_092914_D_cell4';      % To assess where the light stimulus was, and the HOA that contains each cell
params.files = {'2014_10_19_', [3:12]};  % <file name prefix, suffix>
params.groups = {'control', [3,4];...
                 'nbqxGabazine', [6,7];...
                 'NMDAR', [6:12]};
params.excludeHS1 = {};
params.excludeHS2 = {};
params.tags = {};
params.filter = 800;



HS1loc = [0 0];
HS2loc = [0 0];
Pialoc = [-240 360];
params.celldepth = [norm(HS1loc-Pialoc), norm([HS2loc-Pialoc])];


% stuff for E/I and AMPA/NMDA ratios
% key for isolatedCurrents = {<current><group><Vhold><Erev>}
% Erev is to calculate driving force for conversion from pA to pS
params.isolatedCurrents = {'excit', 'control', -72, 17;...
                           'inhib', 'control', 17, -72;...
                           'ampa', 'control', -72, 17;...
                           'nmda', 'nbqxGabazine', 50, 17};


%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_optoIV,  @anlyMod_EIbalance, @anlyMod_NMDAR};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end

%% AK_092914_D Cell 5

fin

% NOTES
%%%%%%%%%%%%%%%%%%%%%%%%
% Pretty nice data from a L2/3 PY cell in PM. 
%
% brain area: PM, PY cell.
% popAnalysis: 
%


%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'AK_092914_D';      % The mouse's name
params.cellNum = 5;    % The neuron number that day
params.photo = 'AK_092914_D_cell5';      % To assess where the light stimulus was, and the HOA that contains each cell
params.files = {'2014_10_19_', [14:24]};  % <file name prefix, suffix>
params.groups = {'control', [14,15];...
                 'nbqxGabazine', [17,18];...
                 'NMDAR', [17:24]};
params.excludeHS1 = {};
params.excludeHS2 = {};
params.tags = {};
params.filter = 2e3;



HS1loc = [0 0];
HS2loc = [0 0];
Pialoc = [54 -184];
params.celldepth = [norm(HS1loc-Pialoc), norm([HS2loc-Pialoc])];


% stuff for E/I and AMPA/NMDA ratios
% key for isolatedCurrents = {<current><group><Vhold><Erev>}
% Erev is to calculate driving force for conversion from pA to pS
params.isolatedCurrents = {'excit', 'control', -72, 17;...
                           'inhib', 'control', 17, -72;...
                           'ampa', 'control', -72, 17;...
                           'nmda', 'nbqxGabazine', 50, 17};


%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_optoIV,  @anlyMod_EIbalance, @anlyMod_NMDAR};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end



%% AK_101314_A Cell 1

fin

% NOTES
%%%%%%%%%%%%%%%%%%%%%%%%
% Okay recording quality, but I forgot to aquire data at -72mV prior to
% drugs, so I can't recover A/N or E/I ratios. The cell fill looks nice
% though... Also, I think that this is and interneuron, but this is based
% on little data. I should consult the cell fill.
%
% brain area: PM
%
% popAnalysis: NMDAR
%
% Layer: 2/3 
%
% cell type: Not exactly sure, but could be a PY cell. There's little
% physiology evidence to back this up (not much spontaneous activity at any
% Vhold, but also not much data).
%


%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'AK_101314_A';      % The mouse's name
params.cellNum = 1;    % The neuron number that day
params.photo = 'AK_101314_A_cell1';      % To assess where the light stimulus was, and the HOA that contains each cell
params.files = {'2014_10_28_', [4:11]};  % <file name prefix, suffix>
params.groups = {'NMDAR', [4:11]};
params.excludeHS1 = {};
params.excludeHS2 = {};
params.tags = {};
params.filter = 2e3;
params.celldepth = [norm([2 -193]), nan];

% stuff for E/I and AMPA/NMDA ratios
% key for isolatedCurrents = {<current><group><Vhold><Erev>}
% Erev is to calculate driving force for conversion from pA to pS
params.isolatedCurrents = {};



HS1loc = [0];
HS2loc = [0];
Pialoc = [2 -193];
params.celldepth = [norm(HS1loc-Pialoc), norm([HS2loc-Pialoc])];


%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_optoIV, @anlyMod_NMDAR};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end


%% AK_101314_A Cell 2

fin

% NOTES
%%%%%%%%%%%%%%%%%%%%%%%%
% This cell was recorded in the same slice as cell1 from AK_101314_A but
% after it had already been exposed to synaptic blockers. I kept the
% blockers on, and recorded this cell. It will only contribute to NMDAR
% analysis, and I don't know what type of cell this is... 
%
% brain area: PM
% popAnalysis: NMDAR
% cell type: not sure, need to check the histology
%


%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'AK_101314_A';      % The mouse's name
params.cellNum = 2;    % The neuron number that day
params.photo = 'AK_101314_A_cell2';      % To assess where the light stimulus was, and the HOA that contains each cell
params.files = {'2014_10_28_', [13:20]};  % <file name prefix, suffix>
params.groups = {'NMDAR', [13:20]};
params.excludeHS1 = {};
params.excludeHS2 = {};
params.tags = {};
params.filter = 2e3;
params.celldepth = [nan, norm([26 -256])];



HS1loc = [0];
HS2loc = [0];
Pialoc = [26 -256];
params.celldepth = [norm(HS1loc-Pialoc), norm([HS2loc-Pialoc])];


% stuff for E/I and AMPA/NMDA ratios
% key for isolatedCurrents = {<current><group><Vhold><Erev>}
% Erev is to calculate driving force for conversion from pA to pS
params.isolatedCurrents = {};


%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_optoIV, @anlyMod_NMDAR};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end



%% AK_101314_C Pair 1

fin

% NOTES
%%%%%%%%%%%%%%%%%%%%%%%%
% Nice paired data from two cells in AL. Need to determine what type of
% cell they are...
%
% brain area: AL
% popAnalysis: EIAN, NMDAR
%


%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'AK_101314_C';      % The mouse's name
params.cellNum = 1;    % The neuron number that day
params.photo = 'AK_101314_C_pair1';      % To assess where the light stimulus was, and the HOA that contains each cell
params.files = {'2014_10_29_', [0:10]};  % <file name prefix, suffix>
params.groups = {'control', [0,1];...
                 'nbqxGabazine', [2,3];...
                 'NMDAR', [4:10]};
params.excludeHS1 = {};
params.excludeHS2 = {};
params.tags = {};
params.filter = 2e3;
params.celldepth = [norm([-263, -403]), norm([-262+56, -389+125])];



HS1loc = [1 14];
HS2loc = [-56 -125];
Pialoc = [-262 -389];
params.celldepth = [norm(HS1loc-Pialoc), norm([HS2loc-Pialoc])];


% stuff for E/I and AMPA/NMDA ratios
% key for isolatedCurrents = {<current><group><Vhold><Erev>}
% Erev is to calculate driving force for conversion from pA to pS
params.isolatedCurrents = {'excit', 'control', -72, 17;...
                           'inhib', 'control', 17, -72;...
                           'ampa', 'control', -72, 17;...
                           'nmda', 'nbqxGabazine', 50, 17};


%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_optoIV,  @anlyMod_EIbalance, @anlyMod_NMDAR};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end


%% AK_101314_C Pair 1 (spatial dependence)

fin

% NOTES
%%%%%%%%%%%%%%%%%%%%%%%%
% These recordings were taken after NBQX and Gabazine had washed out, but
% in normal ACSF. The currents are likely NMDAR mediated though (i think).
% I'm just trying to understand if stimulation location affects the ratio
% of currents measured in different cell types.
%
% brain area: AL
% popAnalysis: none.
%


%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'AK_101314_C';      % The mouse's name
params.cellNum = 100; % NOTE: I'm making this number big so that it doesn't interfere with the params info from the previous analysis cell (same neurons)
params.photo = 'AK_101314_C_pair1';      % To assess where the light stimulus was, and the HOA that contains each cell
params.files = {'2014_10_29_', [11:15]};  % <file name prefix, suffix>
params.groups = {'NMDAR', [11:15]};
params.excludeHS1 = {};
params.excludeHS2 = {};
params.tags = {};
params.filter = 2e3;

params.stimLoc = [-67, -67;...
                  -132 -178;...
                  -209 -284;...
                  -3, 6;...
                  -76 -67];


% stuff for E/I and AMPA/NMDA ratios
% key for isolatedCurrents = {<current><group><Vhold><Erev>}
% Erev is to calculate driving force for conversion from pA to pS
params.isolatedCurrents = {};


%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_optoIV};
    params = invitroAnalysisOverview(params);
end


% % extra analysis for this pair of cells
% peak_pA_HS1 = cellfun(@mean, params.ivdat.NMDAR.peakBySweep_pA{1});
% peak_pA_HS2 = cellfun(@mean, params.ivdat.NMDAR.peakBySweep_pA{2});
% figure, hold on,
% plot(peak_pA_HS1, '-k.')
% plot(peak_pA_HS2, '-b.')



%% BOOKMARK FOR AK MICE

% make sure that dat_DB script does not create duplicate entries, and that
% one run of dat_DB for one script does not intrued (or becomes affected
% by) prior dat_DB runs. To do this, make the mdb from scratch. this will
% add processing time...

% Make sure that the kinetics are the same for each Vhold. It's possible
% that the NMDAR currents are contaminated by Ca2+ current for Vholds close
% to -20 mV.

% make sure that all the physiology_notes have "note" sections and "cell
% type" sections and "brain area" sections

% institute a check to make sure that all the parameters are the same
% acorss a set of .abf files.

% Need to add:
% AK_092914_A (EIAN) Need to look at histology, write notes, and confirm HVA location on the pop_workbook
% AK_092914_B (EIAN) Need to look at histology, write notes, and confirm HVA location on the pop_workbook
% AK_092914_C (EIAN) Need to look at histology, write notes, and confirm HVA location on the pop_workbook
% AK_092914_D (EIAN) Need to look at histology, write notes, and confirm HVA location on the pop_workbook
% AK_101314_A need to verify cell type
% AK_101314_C need to verify brain region

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
params.stimLoc = zeros(5,2);    % The (x,y) coordinates of the objective at the locations stimulated with the LED. <Nx2> matrix
params.legTxt = {'5', '40', '60', '2', '20'};     % Text that will appear in figures to annotate each data file. <Nx1> cell array
params.subtractexp = false;
params.filter = 4000;

%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
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
params.stimLoc = zeros(5,2);    % The (x,y) coordinates of the objective at the locations stimulated with the LED. <Nx2> matrix
params.legTxt = {'20', '40', '5'};     % Text that will appear in figures to annotate each data file. <Nx1> cell array
params.subtractexp = false;
params.filter = 4000;

%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end



%% CH_020314_C Cells 1-4
fin
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

fin

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
params.stimLoc = zeros(3,2);    % The (x,y) coordinates of the objective at the locations stimulated with the LED. <Nx2> matrix
params.legTxt = {'20', '5', '40'};     % Text that will appear in figures to annotate each data file. <Nx1> cell array
params.subtractexp = false;
params.filter = 4000;

%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
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
params.subtractexp = false;
params.filter = 4000;

%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
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
params.subtractexp = false;
params.filter = 4000;

%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end



%% CH_032414_C Cell 3

fin

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
params.subtractexp = false;
params.filter = 4000;

%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end


%% CH_032414_C Cell 4
fin
%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Very small facilitation for 20 Hz at the soma, that turns into  unity
% P1:P2 when the light targets L4. When LED targets L5, than you see
% depression at all TRs. 
%
% This cell appears to by a L2/3 PY cell. Definitely not a FS cell, but
% little Ih. Vm rest is about -82Mv. This neuron is either in LM or AL. The
% histology is a little ambiguous, but there's a modest bias for me to call
% this a LM. 

%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'CH_032414_C';      % The mouse's name
params.cellNum = 4;    % The neuron number that day
params.DCsteps = '2014_04_09_0030';    % DC steps for Rin and cell identification
params.photo = 'cell_4_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.validCh = 'HS2_';    % 'HS2_' or 'HS1_'
params.files = {'2014_04_09_0032';...
                '2014_04_09_0035';...
                '2014_04_09_0036';...
                '2014_04_09_0037';...
                '2014_04_09_0038';...
                '2014_04_09_0039';...
                '2014_04_09_0040';...
                '2014_04_09_0041';...
                '2014_04_09_0042';...
                '2014_04_09_0043';...
                };      % File names of the raw data. <Nx1> cell array
params.skipSweeps = {[],[18]}; % In case I need to ignore certain sweeps
params.vHold = -85.*ones(size(params.files));     % The holding potential for vClamp experiments. One value for each expt. <Nx1>
params.stimLoc = [0,0;...
                  52,190;...
                  52,190;...
                  52,190;...
                  52,190;...
                  52,190;...
                  52,190;...
                  51,341;...
                  51,341;...
                  51,341];    % The (x,y) coordinates of the objective at the locations stimulated with the LED. <Nx2> matrix
params.legTxt = {'Soma', 'Pos 1, TF =20','Pos 1, TF =5','Pos 1, TF =40','Pos 1, TF = 40, 400us','Pos 1, TF = 20, 400us','Pos 1, TF = 5, 400us', 'Pos 2, TF = 5, 400us', 'Pos 2, TF = 20, 400us', 'Pos 2, TF = 40, 400us'};     % Text that will appear in figures to annotate each data file. <Nx1> cell array
params.tags = {'HOA, TF'};
params.subtractexp = false;
params.filter = 4000;

%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end


%% CH_042214_A cell 1

fin

%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% A L4 delayed onset IN. Possibly an FS cell. P2:P1 facilitation at the
% soma for 20 Hz, but depression at the distal stim site. Adding to
% population analysis.


%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'CH_042214_A';      % The mouse's name
params.cellNum = 1;    % The neuron number that day
params.DCsteps = '2014_05_12_0001';    % DC steps for Rin and cell identification
params.photo = 'cell_1_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.validCh = 'HS2_';    % 'HS2_' or 'HS1_'
params.files = {'2014_05_12_0004';...
                '2014_05_12_0006';...
                '2014_05_12_0007';...
                '2014_05_12_0008'};      % File names of the raw data. <Nx1> cell array
params.skipSweeps = {}; % In case I need to ignore certain sweeps
params.vHold = -85 .* ones(numel(params.files));     % The holding potential for vClamp experiments. One value for each expt. <Nx1>
params.stimLoc = [0 0;...
                  -130 238;...
                  -130 238;...
                  -130 238];    % The (x,y) coordinates of the objective at the locations stimulated with the LED. <Nx2> matrix
params.legTxt = {};     % Text that will appear in figures to annotate each data file. <Nx1> cell array
params.tags = {'axon stim'};
params.subtractexp = false;
params.filter = 4000;


%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end


%% CH_042214_A cell 2 part 1

fin

%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% A L2/3 PY cell in PM (According to online notes.) Access resistance is
% potentially a problem. Adding to pop list.


%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'CH_042214_A';      % The mouse's name
params.cellNum = 2;    % The neuron number that day
params.DCsteps = '2014_05_12_0010';    % DC steps for Rin and cell identification
params.photo = 'cell_2_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.validCh = 'HS2_';    % 'HS2_' or 'HS1_'
params.files = {'2014_05_12_0012';...
               '2014_05_12_0014';...
               '2014_05_12_0015';...
               '2014_05_12_0016';...
               '2014_05_12_0017'};      % File names of the raw data. <Nx1> cell array
params.skipSweeps = {}; % In case I need to ignore certain sweeps
params.vHold = -85 .* ones(numel(params.files));     % The holding potential for vClamp experiments. One value for each expt. <Nx1>
params.stimLoc = [0 0;...
                  42 238;...
                  107 394;...
                  107 394;...
                  107 394];    % The (x,y) coordinates of the objective at the locations stimulated with the LED. <Nx2> matrix
params.legTxt = {};     % Text that will appear in figures to annotate each data file. <Nx1> cell array
params.tags = {'Axon stim'};
params.subtractexp = false;
params.filter = 4000;


%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end



%% CH_042214_A cell 2 part 2
fin

%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is only for spatial light spread control. Good data for spatial
% light spread. Looks like some places evoke responses but not others.


%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'CH_042214_A';      % The mouse's name
params.cellNum = 1;    % The neuron number that day
params.DCsteps = '2014_05_12_0010';    % DC steps for Rin and cell identification
params.photo = 'cell_2_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.validCh = 'HS2_';    % 'HS2_' or 'HS1_'
params.files = {'2014_05_12_0018';...
                '2014_05_12_0019';...
                '2014_05_12_0020';...
                '2014_05_12_0021';...
                '2014_05_12_0022';...
                '2014_05_12_0023';...
                '2014_05_12_0024';...
                '2014_05_12_0025'};      % File names of the raw data. <Nx1> cell array
params.skipSweeps = {}; % In case I need to ignore certain sweeps
params.vHold = -85.*ones(size(params.files));     % The holding potential for vClamp experiments. One value for each expt. <Nx1>
params.stimLoc = [46 404;...
                  -22 406;...
                  -116 389;...
                  -185 362;...
                  -237 331;...
                  -210 350;...
                  -256 317;...
                  -46 404];    % The (x,y) coordinates of the objective at the locations stimulated with the LED. <Nx2> matrix
params.legTxt = {};     % Text that will appear in figures to annotate each data file. <Nx1> cell array
params.tags = {'Spatial light spread', 'control'};
params.subtractexp = false;
params.filter = 4000;


%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

%% CH_042214_D Cell 1

fin

%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Spatial light spread. Very nice example of light targets that are
% equidistant but only some of the locations evoke responses.


%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'CH_042214_D';      % The mouse's name
params.cellNum = 1;    % The neuron number that day
params.DCsteps = '2014_05_13_0000';    % DC steps for Rin and cell identification
params.photo = 'CH_042214_D_cell_1_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.validCh = 'HS2_';    % 'HS2_' or 'HS1_'
params.files = {'2014_05_13_0001';...
                '2014_05_13_0002';...
                '2014_05_13_0003';...
                '2014_05_13_0004';...
                '2014_05_13_0005';...
                '2014_05_13_0006';...
                '2014_05_13_0007';...
                '2014_05_13_0009';...
                '2014_05_13_0010';...
                '2014_05_13_0011';...
                '2014_05_13_0012'};      % File names of the raw data. <Nx1> cell array
params.skipSweeps = {}; % In case I need to ignore certain sweeps
params.vHold = -85 .* ones(numel(params.files));     % The holding potential for vClamp experiments. One value for each expt. <Nx1>
params.stimLoc = [50 408;...
                  -29 410;...
                  -83 401;...
                  -133 388;...
                  -185 366;...
                  -240 331;...
                  -290 290;...
                  -330 244;...
                  -363 191;...
                  -387 137;...
                  -402 80];    % The (x,y) coordinates of the objective at the locations stimulated with the LED. <Nx2> matrix
params.legTxt = {};     % Text that will appear in figures to annotate each data file. <Nx1> cell array
params.tags = {'Spatial light spread', 'control'};
params.subtractexp = false;
params.filter = 4000;


%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end


%% CH_042214_D Cell 2

fin

%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Ok data, not great. L2/3 PY cell from AL (according to notes). Ra
% increases... Facilitation at soma, depression otherwise, except at high
% TF. Some of the files have different first pulse heights despite being
% stimuluated at the same site (non-stationarities?)


%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'CH_042214_D';      % The mouse's name
params.cellNum = 2;    % The neuron number that day
params.DCsteps = '2014_05_13_0014';    % DC steps for Rin and cell identification
params.photo = 'CH_042214_D_cell_2_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.validCh = 'HS2_';    % 'HS2_' or 'HS1_'
params.files = {'2014_05_13_0017';...
                '2014_05_13_0019';...
                '2014_05_13_0020';...
                '2014_05_13_0021';...
                '2014_05_13_0022'};      % File names of the raw data. <Nx1> cell array
params.skipSweeps = {}; % In case I need to ignore certain sweeps
params.vHold = -85 .* ones(numel(params.files));     % The holding potential for vClamp experiments. One value for each expt. <Nx1>
params.stimLoc = [0 0;...
                  -166 161;...
                  -188 176;...
                  -184 193;...
                  -184 193];
params.legTxt = {};     % Text that will appear in figures to annotate each data file. <Nx1> cell array
params.tags = {};
params.subtractexp = false;
params.filter = 4000;


%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end


%% CH_042214_E Cell 1

fin

%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Neuron from an FS cell in L2/3. There is a ton of voltage clamp break
% down due to spiking. And this file has several repeats at each spatial
% location using different voltages to the LED. There are subtle STP
% differences for different voltages, but there are some differences too.
% I'm not 100% sure where this neuron is from b/c it's in a sagital slice.
%


%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'CH_042214_E';      % The mouse's name
params.cellNum = 1;    % The neuron number that day
params.DCsteps = '2014_05_16_0001';    % DC steps for Rin and cell identification
params.photo = 'CH_042214_E_cell_1_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.validCh = 'HS2_';    % 'HS2_' or 'HS1_'
params.files = {'2014_05_16_0003';...
                '2014_05_16_0004';...
                '2014_05_16_0006';...
                '2014_05_16_0007';...
                '2014_05_16_0009';...
                '2014_05_16_0010';...
                '2014_05_16_0011';...
                '2014_05_16_0012';...
                '2014_05_16_0013';...
                '2014_05_16_0014';...
                '2014_05_16_0015'};      % File names of the raw data. <Nx1> cell array
params.skipSweeps = {}; % In case I need to ignore certain sweeps
params.vHold = -85 .* ones(numel(params.files));     % The holding potential for vClamp experiments. One value for each expt. <Nx1>
params.stimLoc = [0 0;...
                  0 0;...
                  -33 250;...
                  -33 250;...
                  -40 407;...
                  -40 407;...
                  -40 407;...
                  -40 407;...
                  -40 407;...
                  -40 407;...
                  -40 407];    % The (x,y) coordinates of the objective at the locations stimulated with the LED. <Nx2> matrix
params.legTxt = {};     % Text that will appear in figures to annotate each data file. <Nx1> cell array
params.tags = {'control', 'multiple voltages'};
params.subtractexp = false;
params.filter = 4000;


%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end


%% CH_050614_B cell 1

fin

%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Potentially interesting data. The 5Hz data shows strong depression
% whereas all other TFs show facilitation (regardless of LED stim
% location). Adding to pop anly list. L2/3 PY cell in AL.


%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'CH_050614_B';      % The mouse's name
params.cellNum = 1;    % The neuron number that day
params.DCsteps = '2014_05_19_0001';    % DC steps for Rin and cell identification
params.photo = 'CH_050614_B_cell_1_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.validCh = 'HS2_';    % 'HS2_' or 'HS1_'
params.files = {'2014_05_19_0003';...
                '2014_05_19_0005';...
                '2014_05_19_0006';...
                '2014_05_19_0007';...
                '2014_05_19_0008';...
                '2014_05_19_0009';...
                '2014_05_19_0010'};      % File names of the raw data. <Nx1> cell array
params.skipSweeps = {[], [], [], [3], [], [], [21]}; % In case I need to ignore certain sweeps
params.vHold = -85 .* ones(numel(params.files));     % The holding potential for vClamp experiments. One value for each expt. <Nx1>
params.stimLoc = [0 0;...
                  -110 254;...
                  -110 254;...
                  -110 254;...
                  -211 312;...
                  -211 312;...
                  -211 312];    % The (x,y) coordinates of the objective at the locations stimulated with the LED. <Nx2> matrix
params.legTxt = {};     % Text that will appear in figures to annotate each data file. <Nx1> cell array
params.tags = {'axon stim'};
params.subtractexp = false;
params.filter = 4000;


%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end


%% CH_050614_B cell 2 part 1

fin

%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Only marginal data. Responses are quite weak, and there are some problems
% with series resistance. Otherwise, o.k. L2/3 PY cell from AL. Adding to
% pop anly list as such.

%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'CH_050614_B';      % The mouse's name
params.cellNum = 2;    % The neuron number that day
params.DCsteps = '2014_05_19_0014';    % DC steps for Rin and cell identification
params.photo = 'CH_050614_B_cell_2_tdTomato_bright';      % To assess where the light stimulus was, and the HOA that contains each cell
params.validCh = 'HS2_';    % 'HS2_' or 'HS1_'
params.files = {'2014_05_19_0016';...
                '2014_05_19_0017';...
                '2014_05_19_0018';...
                '2014_05_19_0019';...
                '2014_05_19_0020'};      % File names of the raw data. <Nx1> cell array
params.skipSweeps = {}; % In case I need to ignore certain sweeps
params.vHold = -85 .* ones(size(params.files));     % The holding potential for vClamp experiments. One value for each expt. <Nx1>
params.stimLoc = [0 0;
                  -78 218;
                  -91 330
                  -91 330
                  -91 330];    % The (x,y) coordinates of the objective at the locations stimulated with the LED. <Nx2> matrix
params.legTxt = {};     % Text that will appear in figures to annotate each data file. <Nx1> cell array
params.tags = {'axon stim'};
params.subtractexp = false;
params.filter = [1000];


%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end

%% CH_050614_B cell 2 part 2

fin

%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Control experiment for spatial light spread. Looks good. Equidistant LED
% stimulation sites.

%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'CH_050614_B';      % The mouse's name
params.cellNum = 2;    % The neuron number that day
params.DCsteps = '';    % DC steps for Rin and cell identification
params.photo = 'CH_050614_B_cell_2_tdTomato_bright';      % To assess where the light stimulus was, and the HOA that contains each cell
params.validCh = 'HS2_';    % 'HS2_' or 'HS1_'
params.files = {'2014_05_19_0021';...
                '2014_05_19_0022';...
                '2014_05_19_0023';...
                '2014_05_19_0024';...
                '2014_05_19_0025';...
                '2014_05_19_0026'};      % File names of the raw data. <Nx1> cell array
params.skipSweeps = {[], [], [], [], [],[4]}; % In case I need to ignore certain sweeps
params.vHold = -85 .* ones(size(params.files));     % The holding potential for vClamp experiments. One value for each expt. <Nx1>
params.stimLoc = [-91 330;
                  -170 267;
                  -226 221
                  -272 161
                  -298 104
                  -315 40];    % The (x,y) coordinates of the objective at the locations stimulated with the LED. <Nx2> matrix
params.legTxt = {};     % Text that will appear in figures to annotate each data file. <Nx1> cell array
params.tags = {'spatial light spread', 'control'};
params.subtractexp = false;
params.filter = 500;


%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

%% CH_063014_C Pair 1

fin

% NOTES
%%%%%%%%%%%%%%%%%%%%%%%%
% First file looking at A/N ratios and E/I ratios. There isn't enough data
% to do the full analysis, but these data are useful for debugging code.
%
% brain area: I think this is at the LM/AL border. I'm not 100% certian.
% popAnalysis: E/I & A/N ratios.
%

%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'CH_063014_C';      % The mouse's name
params.cellNum = 1;    % The neuron number that day
params.photo = 'CH_063014_C_pair1_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.files = {'2014_07_19_', [0:16]};  % <file name prefix, suffix>
params.groups = {'control', [0:7];...
                 'nbqxGabazine', [8:16];...
                 'NMDAR', [8:16]};
params.excludeHS1 = {};
params.excludeHS2 = {'_0007', '_0016'};
params.tags = {};
params.filter = 800;


HS1loc = [nan];
HS2loc = [nan];
Pialoc = [nan];
params.celldepth = [norm(HS1loc-Pialoc), norm([HS2loc-Pialoc])];


% stuff for E/I and AMPA/NMDA ratios
% key for isolatedCurrents = {<current><group><Vhold><Erev>}
% Erev is to calculate driving force for conversion from pA to pS
params.isolatedCurrents = {'excit', 'control', -80, 15;...
                           'inhib', 'control', 20, -75;...
                           'ampa', 'control', -80, 15;...
                           'nmda', 'nbqxGabazine', 40, 15};


%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_optoIV,  @anlyMod_EIbalance, @anlyMod_NMDAR};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end


%% CH_070714_A Pair 1

fin

%
% Notes
%%%%%%%%%%%%%%%%%
% This pair of cells has o.k data for A/N and E/I ratios. The cell 2 has
% crummier data and should not be included.
%
% Brain area: unknown. Very posterior, so it could be LM or even V1. It's
% a little unclear and the histology doesn't help much
%
% popAnalysis: Added to E/I and A/N for cell 1

%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'CH_070714_A';      % The mouse's name
params.cellNum = 1;    % The neuron number that day
params.photo = 'CH_070714_A_pair1_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.files = {'2014_07_23_', [0:13]};  % <file name prefix, suffix>
params.groups = {'control', [0:6];...
                 'nbqxGabazine', [7:13];...
                 'NMDAR', [7:13]};
params.excludeHS1 = {};
params.excludeHS2 = {'_0004', '_0005', '_0006'};
params.tags = {};
params.filter = 800;


HS1loc = [nan];
HS2loc = [nan];
Pialoc = [nan];
params.celldepth = [norm(HS1loc-Pialoc), norm([HS2loc-Pialoc])];


% stuff for E/I and AMPA/NMDA ratios
% key for isolatedCurrents = {<current><group><Vhold><Erev>}
% Erev is to calculate driving force for conversion from pA to pS
params.isolatedCurrents = {'excit', 'control', -75, 15;...
                           'inhib', 'control', 15, -75;...
                           'ampa', 'control', -75, 15;...
                           'nmda', 'nbqxGabazine', 40, 15};


%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_optoIV, @anlyMod_EIbalance, @anlyMod_NMDAR};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end



%% CH_070714_C Pair 1

fin
%
% Notes
%
%%%%%%%%%%%%%%%%%
% E/I and A/N ratios from two cells in L2/3 of an unknown HVA. The neuron
% on channel 2 is pretty crummy and escapes Vclamp. I omited the spikes
% from this channel, but the results after the LED pulse still contaminate
% the recording. Channel one has good data that I would like to use...
%
% brain area: unknown. The hisology is out of order, and inconclusive. I
% could go back to scrutinize...
%
% popAnalysis: E/I and A/N. 
%





%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'CH_070714_C';      % The mouse's name
params.cellNum = 1;    % The neuron number that day
params.photo = 'CH_070714_C_pair1_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.files = {'2014_07_23_', [14:29]};  % <file name prefix, suffix>
params.groups = {'control', [14:19];...
                 'nbqxGabazine', [20:29]};
params.excludeHS1 = {{'_0014', [8]}};
params.excludeHS2 = {{'_0014', [8]}, '_0016', {'_0017',[1,4,5,7,10,11,14,18,19,27,31,35,37,38]},...
                     {'_0017', [12]}, {'_0025', [2,5,6,7,8,9,14,17,22,23,25,26,28,33,36]}, '_0026', '_0027', '_0028', '_0029'};

% stuff for E/I and AMPA/NMDA ratios
% key for isolatedCurrents = {<current><group><Vhold><Erev>}
% Erev is to calculate driving force for conversion from pA to pS
params.isolatedCurrents = {'excit', 'control', -75, 15;...
                           'inhib', 'control', 15, -75;...
                           'ampa', 'control', -75, 15;...
                           'nmda', 'nbqxGabazine', 50, 15};
params.tags = {};
params.filter = 2e3;


params.celldepth = [180 180];



%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_optoIV, @anlyMod_EIbalance};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end



%% CH_071414_A Pair 1

fin
%
% Notes
%
%%%%%%%%%%%%%%%%%
% pretty decent data from 2 neurons in L2/3. The neuron on ch 2 is an
% interneuron (As assessed by the amount of spontaneous PSCs and the E/I
% balance).
%
% brain area: possibly LM b/c I think this slice is too posterior to be AL
% popAnalysis: adding to E/I and A/N ratios for both cells
%

%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'CH_071414_A';      % The mouse's name
params.cellNum = 1;    % The neuron number that day
params.photo = 'CH_071414_A_pair1_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.files = {'2014_08_01_', [0:14]};  % <file name prefix, suffix>
params.groups = {'control', [0:6];...
                 'nbqxGabazine', [7:14];...
                 'NMDAR', [7:14]};
params.excludeHS1 = {{'_0006', [3,6,10,14]}, {'_0005', [10,16,24]}, {'_0004', [7]}, '_0004'};
params.excludeHS2 = {};

% stuff for E/I and AMPA/NMDA ratios
% key for isolatedCurrents = {<current><group><Vhold><Erev>}
% Erev is to calculate driving force for conversion from pA to pS
params.isolatedCurrents = {'excit', 'control', -75, 15;...
                           'inhib', 'control', 15, -75;...
                           'ampa', 'control', -75, 15;...
                           'nmda', 'nbqxGabazine', 50, 15};
params.tags = {};
params.filter = 2e3;



HS1loc = [-99 13];
HS2loc = [-27 -27];
Pialoc = [-153 276];
params.celldepth = [norm(HS1loc-Pialoc), norm([HS2loc-Pialoc])];


%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_optoIV, @anlyMod_EIbalance, @anlyMod_NMDAR};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end


%% CH_071414_D pair 1

fin
%
% Notes
%
%%%%%%%%%%%%%%%%%
% Pretty good data from 2 L2/3 PY cells (i think). Unfortunately, the
% injection site was very medial, and I think hit PM more than V1. I'm
% guessing that the HVA I recorded from was LM. Also, variability in CH1
% may underestimate the inhibition present.
%
%
% brain area: LM?
% popAnalysis: adding, but with notes about injection issues. E/I and A/N.
%

%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'CH_071414_D';      % The mouse's name
params.cellNum = 1;    % The neuron number that day
params.photo = 'CH_071414_D_pair1_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.files = {'2014_07_30_', [0:4, 6:13]};  % <file name prefix, suffix>
params.groups = {'control', [1:4,6];...
                 'nbqxGabazine', [7:13];...
                 'NMDAR', [7:13]};
params.excludeHS1 = {'_0012', '_0013'};
params.excludeHS2 = {{'_0011' [34:40]}, '_0012', '_0013'};

% stuff for E/I and AMPA/NMDA ratios
% key for isolatedCurrents = {<current><group><Vhold><Erev>}
% Erev is to calculate driving force for conversion from pA to pS
params.isolatedCurrents = {'excit', 'control', -75, 20;...
                           'inhib', 'control', 20, -75;...
                           'ampa', 'control', -75, 20;...
                           'nmda', 'nbqxGabazine', 50, 20};
params.tags = {};
params.filter = 2e3;


HS1loc = [-6 23];
HS2loc = [37 21];
Pialoc = [10 211];
params.celldepth = [norm(HS1loc-Pialoc), norm([HS2loc-Pialoc])];



%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_optoIV, @anlyMod_EIbalance, @anlyMod_NMDAR};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end



%% CH_072214_A pair 1

fin
%
% Notes
%
%%%%%%%%%%%%%%%%%
% Decent data from CH2, but CH1 is junky. I'm assuming that cel 2 is an IN
% cell.
%
% brain area: I'm guessing this is area LM.
% popAnalysis: Added cell 2 to the list for E/I and A/N
%

%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'CH_072214_A';      % The mouse's name
params.cellNum = 1;    % The neuron number that day
params.photo = 'CH_072214_A_pair1_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.files = {'2014_08_08_', [4:6,12,13,15]};  % <file name prefix, suffix>
params.groups = {'control', [4:6];...
                 'nbqxGabazine', [12,13,15]};
params.excludeHS1 = {'_0006', '_0012', '_0015'};
params.excludeHS2 = {{'_0006', [1,2,4,6,10]}};

% stuff for E/I and AMPA/NMDA ratios
% key for isolatedCurrents = {<current><group><Vhold><Erev>}
% Erev is to calculate driving force for conversion from pA to pS
params.isolatedCurrents = {'excit', 'control', -75, 15;...
                           'inhib', 'control', 25, -75;...
                           'ampa', 'control', -75, 15;...
                           'nmda', 'nbqxGabazine', 50, 15};
params.tags = {};
params.filter = 2e3;


HS1loc = [-68 24];
HS2loc = [3 6];
Pialoc = [-28 -352];
params.celldepth = [norm(HS1loc-Pialoc), norm([HS2loc-Pialoc])];



%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_optoIV, @anlyMod_EIbalance};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end


%% CH_072214_B pair 1

fin
%
% Notes
%
%%%%%%%%%%%%%%%%%
% Good data from CH2. Junk data from CH1. CH2 has solid acess throughout
% and huge inhibition relative to excitation. 
%
% brain area: LM I think.
% popAnalysis: including CH2 to E/I and AMPA/NMDA
%

%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'CH_072214_B';      % The mouse's name
params.cellNum = 1;    % The neuron number that day
params.photo = 'CH_072214_B_pair1_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.files = {'2014_08_06_', [1:12]};  % <file name prefix, suffix>
params.groups = {'control', [1:6];...
                 'nbqxGabazine', [7:12]};
params.excludeHS1 = {{'_0005', [9:35]}, '_0006',{'_0011', [10:40]}, '_0012'};
params.excludeHS2 = {{'_0001', [4]}};

% stuff for E/I and AMPA/NMDA ratios
% key for isolatedCurrents = {<current><group><Vhold><Erev>}
% Erev is to calculate driving force for conversion from pA to pS
params.isolatedCurrents = {'excit', 'control', -75, 15;...
                           'inhib', 'control', 15, -75;...
                           'ampa', 'control', -75, 15;...
                           'nmda', 'nbqxGabazine', 50, 15};
params.tags = {};
params.filter = 500;


HS1loc = [-10 29];
HS2loc = [48 2];
Pialoc = [70 520];
params.celldepth = [norm(HS1loc-Pialoc), norm([HS2loc-Pialoc])];



%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_optoIV, @anlyMod_EIbalance};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end


%% CH_072214_C pair 1

fin
%
% Notes
%
%%%%%%%%%%%%%%%%%
% Good data from CH2. Less good data from CH1. CH2 has solid acess
% throughout. Both cells could be INs...
%
% brain area: It's possible that this is in the injection site, or PM,
% tough to say, and I don't think it should contribute to much.
%
% popAnalysis: including CH2 to E/I and AMPA/NMDA but not writing down a
% brain area.
%

%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'CH_072214_C';      % The mouse's name
params.cellNum = 1;    % The neuron number that day
params.photo = 'CH_072214_C_pair1_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.files = {'2014_08_07_', [0:11]};  % <file name prefix, suffix>
params.groups = {'control', [1:4];...
                 'nbqxGabazine', [5:11]};
params.excludeHS1 = {{'_0001', [12,13,16,21:30,32,33]}, '_0002', {'_0010', [25:40]}};
params.excludeHS2 = {{'_0004', [22:40]}, {'_0009', [2,3]}};

% stuff for E/I and AMPA/NMDA ratios
% key for isolatedCurrents = {<current><group><Vhold><Erev>}
% Erev is to calculate driving force for conversion from pA to pS
params.isolatedCurrents = {'excit', 'control', -75, 15;...
                           'inhib', 'control', 15, -75;...
                           'ampa', 'control', -75, 15;...
                           'nmda', 'nbqxGabazine', 50, 15};
params.tags = {};
params.filter = 2e3;


HS1loc = [-26 -3];
HS2loc = [2 10];
Pialoc = [-12 270];
params.celldepth = [norm(HS1loc-Pialoc), norm([HS2loc-Pialoc])];



%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_optoIV, @anlyMod_EIbalance};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end



%% CH_072214_D cell 1

fin
%
% Notes
%
%%%%%%%%%%%%%%%%%
% Less than stelar data from a single PY cell in PM. No drugs applied, so I
% only measured E/I ratio. Access resistance isn't great.
%
% brain area: PM
% popAnalysis: Adding one cell. Will only contribute to E/I analysis.
%

%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'CH_072214_D';      % The mouse's name
params.cellNum = 1;    % The neuron number that day
params.photo = 'CH_072214_D_cell_1_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.files = {'2014_08_05_', [4:7]};  % <file name prefix, suffix>
params.groups = {'control', [4:7]};
params.excludeHS1 = {};
params.excludeHS2 = {};

% stuff for E/I and AMPA/NMDA ratios
% key for isolatedCurrents = {<current><group><Vhold><Erev>}
% Erev is to calculate driving force for conversion from pA to pS
params.isolatedCurrents = {'excit', 'control', -75, 15;...
                           'inhib', 'control', 15, -75};
params.tags = {};
params.filter = 2e3;


HS1loc = [0];
HS2loc = [0];
Pialoc = [-85 238];
params.celldepth = [norm(HS1loc-Pialoc), norm([HS2loc-Pialoc])];



%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_optoIV, @anlyMod_EIbalance};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end

%% CH_081114_A pair 1

fin
%
% Notes
%
%%%%%%%%%%%%%%%%%
% Paired recording from a SOM+ cell (HS1) and a PY cell (HS2). There is
% basically no recruitment of the SOM+ cell at the LED power used in these
% files (Although there are other voltages that do recruit the SOM+ cell).
% No data for A/N ratio. I'm adding both cells to the population spread
% sheet. L2/3 cells.
% 
% brain area: AL
% popAnalysis: Adding coth cells. Will only contribute to E/I analysis.
%

%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'CH_081114_A';      % The mouse's name
params.cellNum = 1;    % The neuron number that day
params.photo = 'CH_081114_A_pair1_tdTomato_5x';      % To assess where the light stimulus was, and the HOA that contains each cell
params.files = {'2014_08_27_', [0,1]};  % <file name prefix, suffix>
params.groups = {'control', [0,1]};
params.excludeHS1 = {};
params.excludeHS2 = {};

% stuff for E/I and AMPA/NMDA ratios
% key for isolatedCurrents = {<current><group><Vhold><Erev>}
% Erev is to calculate driving force for conversion from pA to pS
params.isolatedCurrents = {'excit', 'control', -75, 20;...
                           'inhib', 'control', 20, -75};
params.tags = {};
params.filter = 2e3;


HS1loc = [15 5];
HS2loc = [58 20];
Pialoc = [-216 152];
params.celldepth = [norm(HS1loc-Pialoc), norm([HS2loc-Pialoc])];



%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_optoIV, @anlyMod_EIbalance};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end


%% CH_081114_A pair 2

fin
%
% Notes
%
%%%%%%%%%%%%%%%%%
% Paired recording from a PY cell (HS1) and a SOM+ cell (HS2). Both cells
% are in L4. There is basically no FF exictation onto the PY cell, but a
% small FF drive to the SOM+ cell with a latency of about 2 ms. At 4 ms,
% the SOM+ cell receives strong recurrent excitation. The PY cell has only
% recurrent excitation. No direct projection. Using stronger LED pulses, I
% measured the NMDA IV curve. To the extent I trust the Vclamp errors,
% there appears to be a difference in the NMDAR subtype expressed on the
% SOM+ cell.
% 
% brain area: PM
% popAnalysis: Adding both cells to E/I, A/N, and NMDAR pop files.
%

%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'CH_081114_A';      % The mouse's name
params.cellNum = 2;    % The neuron number that day
params.photo = 'CH_081114_A_pair2_tdTomato_5x';      % To assess where the light stimulus was, and the HOA that contains each cell
params.files = {'2014_08_27_', [12,13,15,16,19:25]};  % <file name prefix, suffix>
params.groups = {'control', [12,13];...
                 'nbqxGabazine', [15,16];...
                 'NMDAR', [19:25]};
params.excludeHS1 = {};
params.excludeHS2 = {{'_0013', [2]}};

% stuff for E/I and AMPA/NMDA ratios
% key for isolatedCurrents = {<current><group><Vhold><Erev>}
% Erev is to calculate driving force for conversion from pA to pS
params.isolatedCurrents = {'excit', 'control', -75, 15;...
                           'inhib', 'control', 15, -75;...
                           'ampa', 'control', -75, 15;...
                           'nmda', 'nbqxGabazine', 50, 15};
params.tags = {};
params.filter = 2e3;


HS1loc = [-16 7];
HS2loc = [-7 -7];
Pialoc = [-351 -295];
params.celldepth = [norm(HS1loc-Pialoc), norm([HS2loc-Pialoc])];



%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_optoIV, @anlyMod_EIbalance, @anlyMod_NMDAR};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end


%% CH_081114_B cell 1

fin
%
% Notes
%
%%%%%%%%%%%%%%%%%
% A single PY cell from L2/3. The access is okay, but for some reason there
% is still a substantial NMDAR current at +25 indicating that I was not at
% Erev for exciation when measureing inhibition. This could also be due to
% space clamp issues. Also, there is still a measureable NMDAR current at
% -75, which again, could be due to space clamp issues.
% 
% brain area: AL
% popAnalysis: Adding to EIAN and NMDAR analysis pending approval.
%

%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'CH_081114_B';      % The mouse's name
params.cellNum = 1;    % The neuron number that day
params.photo = 'CH_081114_B_cell1_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.files = {'2014_08_25_', [0:10]};  % <file name prefix, suffix>
params.groups = {'control', [0,1];...
                 'nbqxGabazine', [2,3,9];...
                 'NMDAR', [2:10]};
params.excludeHS1 = {{'_0010' [5,6]}};
params.excludeHS2 = {}; % all files are junk. no cell on HS2

% stuff for E/I and AMPA/NMDA ratios
% key for isolatedCurrents = {<current><group><Vhold><Erev>}
% Erev is to calculate driving force for conversion from pA to pS
params.isolatedCurrents = {'excit', 'control', -75, 25;...
                           'inhib', 'control', 25, -75;...
                           'ampa', 'control', -75, 25;...
                           'nmda', 'nbqxGabazine', 60, 25};
params.tags = {};
params.filter = 2e3;


HS1loc = [-33 38];
HS2loc = [nan];
Pialoc = [184 274];
params.celldepth = [norm(HS1loc-Pialoc), norm([HS2loc-Pialoc])];



%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_optoIV, @anlyMod_EIbalance, @anlyMod_NMDAR};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end

%% CH_081114_B cell 2

fin
%
% Notes
%
%%%%%%%%%%%%%%%%%
% A single SOM+ cell in L2/3. Decent data for E/I data, but crummy for A/N
% (no data for A/N). Also, the excitation measured is mostly recurrent
% excitation (on the basis of latency). At a latency <2ms, there is only a
% dinky little response.
%
% brain area: AL
% popAnalysis: adding to EIAN, but this cell will only contribute to E/I
%

%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'CH_081114_B';      % The mouse's name
params.cellNum = 2;    % The neuron number that day
params.photo = 'CH_081114_B_cell2_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.files = {'2014_08_25_', [12:15]};  % <file name prefix, suffix>
params.groups = {'control', [12,13]};
params.excludeHS1 = {};
params.excludeHS2 = {}; % all files are junk. no cell on HS2

% stuff for E/I and AMPA/NMDA ratios
% key for isolatedCurrents = {<current><group><Vhold><Erev>}
% Erev is to calculate driving force for conversion from pA to pS
params.isolatedCurrents = {'excit', 'control', -75, 15;...
                           'inhib', 'control', 15, -75};
params.tags = {};
params.filter = 2e3;


HS1loc = [nan];
HS2loc = [6 11];
Pialoc = [175 238];
params.celldepth = [norm(HS1loc-Pialoc), norm([HS2loc-Pialoc])];



%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_optoIV, @anlyMod_EIbalance};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end


%% CH_081114_C pair 1

fin
%
% Notes
%
%%%%%%%%%%%%%%%%%
% A PY cell and a SOM+ cell. I didn't aquire enough data to include this
% pair in a population analysis. Nothing much to note...


%% CH_081114_C pair 2

fin
%
% Notes
%
%%%%%%%%%%%%%%%%%
% This dataset should be treated with sceptisism. There aren't many sweeps,
% and the holding current on HS2 (PY) is rather large, resulting in Vclamp
% errors as high as 17 mV. Nonetheless, the SOM+ cell on HS1 seems to have
% no inhibition, but lots of excitation... 
%
% brain area: possibly AL, but could also be anterior part of LM
% popAnalysis: Adding tentatively to EIAN.
%

%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'CH_081114_C';      % The mouse's name
params.cellNum = 2;    % The neuron number that day
params.photo = 'CH_081114_C_pair2_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.files = {'2014_08_26_', [6,7]};  % <file name prefix, suffix>
params.groups = {'control', [6,7]};
params.excludeHS1 = {};
params.excludeHS2 = {}; % all files are junk. no cell on HS2

% stuff for E/I and AMPA/NMDA ratios
% key for isolatedCurrents = {<current><group><Vhold><Erev>}
% Erev is to calculate driving force for conversion from pA to pS
params.isolatedCurrents = {'excit', 'control', -75, 15;...
                           'inhib', 'control', 15, -75};
params.tags = {};
params.filter = 2e3;


HS1loc = [nan];
HS2loc = [nan];
Pialoc = [nan];
params.celldepth = [norm(HS1loc-Pialoc), norm([HS2loc-Pialoc])];


%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_optoIV, @anlyMod_EIbalance};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end


%% CH_081114_D pair 1

fin
%
% Notes
%
%%%%%%%%%%%%%%%%%
% A paired recording from a SOM cell (HS1) and a PY cell. The series
% resistance on the SOM cells is really high and I didn't estimate the Erev
% for excitation well in this cell, so I'm only going to continue to
% analyze the PY cell.%
%
% brain area: PM
% popAnalysis: PY cell added to EIAN (E/I only)
%

%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'CH_081114_D';      % The mouse's name
params.cellNum = 1;    % The neuron number that day
params.photo = 'CH_081114_D_pair1_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.files = {'2014_09_02_', [0,1]};  % <file name prefix, suffix>
params.groups = {'control', [0,1]};
params.excludeHS1 = {};
params.excludeHS2 = {}; % all files are junk. no cell on HS2

% stuff for E/I and AMPA/NMDA ratios
% key for isolatedCurrents = {<current><group><Vhold><Erev>}
% Erev is to calculate driving force for conversion from pA to pS
params.isolatedCurrents = {'excit', 'control', -70, 15;...
                           'inhib', 'control', 15, -70};
params.tags = {};
params.filter = 2e3;


HS1loc = [-21 -36];
HS2loc = [-7 -7];
Pialoc = [-82 -242];
params.celldepth = [norm(HS1loc-Pialoc), norm([HS2loc-Pialoc])];



%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_optoIV, @anlyMod_EIbalance};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end


%% CH_081114_D pair 2

fin
%
% Notes
%
%%%%%%%%%%%%%%%%%
%
% not a lot here. Paired recording from a PY cell (HS1) and a SOM cell
% (HS2). There were two different LED voltages used, but the Ra changes
% dramatically during one of them, so comparisons between the two voltages
% are unreliable. The SOM cell may get culled from the analysis due to Ra
% issues. 
%
% brain area: PM
% popAnalysis: Adding to EIAN. SOM cell might get culled due to Ra.
%

%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'CH_081114_D';      % The mouse's name
params.cellNum = 2;    % The neuron number that day
params.photo = 'CH_081114_D_pair2_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.files = {'2014_09_02_', [4:5]};  % <file name prefix, suffix>
params.groups = {'control', [4,5]};
params.excludeHS1 = {};
params.excludeHS2 = {}; % all files are junk. no cell on HS2

% stuff for E/I and AMPA/NMDA ratios
% key for isolatedCurrents = {<current><group><Vhold><Erev>}
% Erev is to calculate driving force for conversion from pA to pS
params.isolatedCurrents = {'excit', 'control', -70, 15;...
                           'inhib', 'control', 15, -70};
params.tags = {};
params.filter = 2e3;


HS1loc = [23 17];
HS2loc = [53 36];
Pialoc = [8 -198];
params.celldepth = [norm(HS1loc-Pialoc), norm([HS2loc-Pialoc])];



%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_optoIV, @anlyMod_EIbalance};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end


%% CH_081914_B all cells

%
% Notes
%
%%%%%%%%%%%%%%%%%
%
% This is a PV+ mouse that I used to estimate the Erev for Cl using my
% Cs-Gluconate internal. These data files haven't been analyzed thorougly,
% but the Erev is very close to -70mV. 
%
%



%% CH_090414_A pair 1

fin
%
% Notes
%
%%%%%%%%%%%%%%%%%
% Pretty mediocre data, but possibly useable nonetheless. Paired recording
% from SOM+ and PY cell. Both cells look like they get direct excitation at
% 2 ms, and inhibition at about 4ms.
%
% brain area: AL I think. My notes don't specify much, but based off of the
% shape of the axon cloud, and the location of the harp string scars in the
% white images, I think that it's impossible that this area is PM. 
%
% popAnalysis: adding to E/I, A/N analysis. Could be culled due to Ra.
%

%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'CH_090414_A';      % The mouse's name
params.cellNum = 1;    % The neuron number that day
params.photo = 'CH_090414_A_pair1_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.files = {'2014_09_19_', [2,4]};  % <file name prefix, suffix>
params.groups = {'control', [2,4]};
params.excludeHS1 = {};
params.excludeHS2 = {}; % all files are junk. no cell on HS2

% stuff for E/I and AMPA/NMDA ratios
% key for isolatedCurrents = {<current><group><Vhold><Erev>}
% Erev is to calculate driving force for conversion from pA to pS
params.isolatedCurrents = {'excit', 'control', -72, 15;...
                           'inhib', 'control', 15, -72};
params.tags = {};
params.filter = 2e3;


HS1loc = [nan];
HS2loc = [nan];
Pialoc = [nan];
params.celldepth = [norm(HS1loc-Pialoc), norm([HS2loc-Pialoc])];


%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_optoIV, @anlyMod_EIbalance};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end


%% CH_090414_A pair 2

fin
%
% Notes
%
%%%%%%%%%%%%%%%%%
% Paired recording from a SOM cell (HS1) and a PY cell (HS2). Decent data.
% Adding to pop analysis for E/I but HS2 cell will not work yet b/c it was
% recorded at a different Vhold than HS1 for the Iinhib case (+15 vs +17
% mV). 
%
% brain area: AL i think
% popAnalysis: adding both cells to E/I pop sheet.
%

%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'CH_090414_A';      % The mouse's name
params.cellNum = 2;    % The neuron number that day
params.photo = 'CH_090414_A_pair2_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.files = {'2014_09_19_', [5,6]};  % <file name prefix, suffix>
params.groups = {'control', [5,6]};
params.excludeHS1 = {};
params.excludeHS2 = {}; % all files are junk. no cell on HS2

% stuff for E/I and AMPA/NMDA ratios
% key for isolatedCurrents = {<current><group><Vhold><Erev>}
% Erev is to calculate driving force for conversion from pA to pS
params.isolatedCurrents = {'excit', 'control', -72,     [15 17];...
                           'inhib', 'control', [15 17], -72     };
params.tags = {};
params.filter = 2e3;


HS1loc = [-83 39];
HS2loc = [6 25];
Pialoc = [-484 246];
params.celldepth = [norm(HS1loc-Pialoc), norm([HS2loc-Pialoc])];



%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_optoIV, @anlyMod_EIbalance};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end


%% CH_090414_C pair 1

fin
%
% Notes
%
%%%%%%%%%%%%%%%%%
% Paired recording from a SOM cell (HS1) and a PY cell (HS2). Decent data.
% Adding to pop analysis for E/I and A/N, and NMDAR.
%
% brain area: This cell may be in LM or AL or at the border...
% popAnalysis: adding both cells to E/I pop sheet.
%

%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'CH_090414_C';      % The mouse's name
params.cellNum = 1;    % The neuron number that day
params.photo = 'CH_090414_C_pair1_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.files = {'2014_09_24_', [1:8, 12]};  % <file name prefix, suffix>
params.groups = {'control', [1,2];...
                 'nbqxGabazine', [3,4];...
                 'NMDAR', [3:8,12]};
params.excludeHS1 = {};
params.excludeHS2 = {'_0012'}; % all files are junk. no cell on HS2

% stuff for E/I and AMPA/NMDA ratios
% key for isolatedCurrents = {<current><group><Vhold><Erev>}
% Erev is to calculate driving force for conversion from pA to pS
params.isolatedCurrents = {'excit', 'control', [-75, -76], 15;...
                           'inhib', 'control', 15, -75;...
                           'ampa', 'control', [-75, -76], 15;...
                           'nmda', 'nbqxGabazine', 50, 15};
params.tags = {};
params.filter = 2e3;


HS1loc = [43 -67];
HS2loc = [-8 11];
Pialoc = [243 -219];
params.celldepth = [norm(HS1loc-Pialoc), norm([HS2loc-Pialoc])];



%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_optoIV, @anlyMod_EIbalance, @anlyMod_NMDAR};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end



%% CH_090414_C cell 2

fin
%
% Notes
%
%%%%%%%%%%%%%%%%%
% Recording from a single PY cell. Nice data for the most part. Access is
% poor for the later files
%
% brain area: PM.
% popAnalysis: adding to E/I A/N, NMDAR pop sheet.
%

%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'CH_090414_C';      % The mouse's name
params.cellNum = 2;    % The neuron number that day
params.photo = 'CH_090414_C_cell2_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.files = {'2014_09_24_', [19:28]};  % <file name prefix, suffix>
params.groups = {'control', [19,20];...
                 'nbqxGabazine', [21,22];...
                 'NMDAR', [21:28]};
params.excludeHS1 = {};
params.excludeHS2 = {}; % all files are junk. no cell on HS2

% stuff for E/I and AMPA/NMDA ratios
% key for isolatedCurrents = {<current><group><Vhold><Erev>}
% Erev is to calculate driving force for conversion from pA to pS
params.isolatedCurrents = {'excit', 'control', -73, 15;...
                           'inhib', 'control', 15, -73;...
                           'ampa', 'control', -73, 15;...
                           'nmda', 'nbqxGabazine', 50, 15};
params.tags = {};
params.filter = 2e3;


HS1loc = [nan];
HS2loc = [3 0];
Pialoc = [-26 262];
params.celldepth = [norm(HS1loc-Pialoc), norm([HS2loc-Pialoc])];



%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_optoIV, @anlyMod_EIbalance, @anlyMod_NMDAR};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end


%% CH_091114_F pair 1

fin
%
% Notes
%
%%%%%%%%%%%%%%%%%
% Pretty good data from 2 PY cells. The Ra and holding current are high on
% some files though, and the Vclamp error creeps up.
%
% brain area: AL
% popAnalysis: E/I, A/N, NMDAR
%

%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'CH_091114_F';      % The mouse's name
params.cellNum = 1;    % The neuron number that day
params.photo = 'CH_091114_F_pair1_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.files = {'2014_09_26_', [1:9]};  % <file name prefix, suffix>
params.groups = {'control', [1,2];...
                 'nbqxGabazine', [3:4];...
                 'NMDAR', [3:9]};
params.excludeHS1 = {{'_0002', [18]}, '_0007', '_0008', '_0009'};
params.excludeHS2 = {{'_0006', [2:13,15,17,19,21,23,25]}}; % all files are junk. no cell on HS2

% stuff for E/I and AMPA/NMDA ratios
% key for isolatedCurrents = {<current><group><Vhold><Erev>}
% Erev is to calculate driving force for conversion from pA to pS
params.isolatedCurrents = {'excit', 'control', -72, 15;...
                           'inhib', 'control', 15, -72;...
                           'ampa', 'control', -72, 15;...
                           'nmda', 'nbqxGabazine', 50, 15};
params.tags = {};
params.filter = 2e3;


HS1loc = [-46 -10];
HS2loc = [-6 9];
Pialoc = [247 -120];
params.celldepth = [norm(HS1loc-Pialoc), norm([HS2loc-Pialoc])];


%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_optoIV, @anlyMod_EIbalance, @anlyMod_NMDAR};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end



%% CH_100614_B cell 1

fin

%
% Notes
%
%%%%%%%%%%%%%%%%%
% Pretty good data, although the NMDAR data at Vhold = -20mV looks
% contaminated by Ca2+ currents. The kinetics at -20 are much slower than
% the others, and there were lots of cases where the Vclamp escaped.
%
% brain area: PM
% popAnalysis: E/I, A/N, NMDAR
%

%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'CH_100614_B';      % The mouse's name
params.cellNum = 1;    % The neuron number that day
params.photo = 'CH_100614_B_cell1_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.files = {'2014_10_22_', [0:4,6:10]};  % <file name prefix, suffix>
params.groups = {'control', [0,1];...
                 'nbqxGabazine', [2,3];...
                 'NMDAR', [2,3,4,6:10]};
params.excludeHS1 = {{'_0000', [14]}, {'_0010', [4,5,7,9,11,12,13,15,18,19,24,25,26,29,30]}};
params.excludeHS2 = {}; 
params.celldepth = [norm(-32 -232), nan];

% stuff for E/I and AMPA/NMDA ratios
% key for isolatedCurrents = {<current><group><Vhold><Erev>}
% Erev is to calculate driving force for conversion from pA to pS
params.isolatedCurrents = {'excit', 'control', -72, 17;...
                           'inhib', 'control', 17, -72;...
                           'ampa', 'control', -72, 17;...
                           'nmda', 'nbqxGabazine', 50, 17};
params.tags = {};
params.filter = 2e3;

HS1loc = [0];
HS2loc = [0];
Pialoc = [-32 -232];
params.celldepth = [norm(HS1loc-Pialoc), norm([HS2loc-Pialoc])];



%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_optoIV, @anlyMod_EIbalance, @anlyMod_NMDAR};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end



%% CH_100614_C cell 2

fin

%
% Notes
%
%%%%%%%%%%%%%%%%%
% Not great but not terrible data from a single PY cell. Using Rs
% compensation on an otherwise dismal Rs. I also think that the Rs
% compensation forced the sec_Vm recording to be a few mV different than
% what I expected.
%
% brain area: PM
% popAnalysis: EIAN, NMDAR
%

%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'CH_100614_C';      % The mouse's name
params.cellNum = 2;    % The neuron number that day
params.photo = 'CH_100614_C_cell2_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.files = {'2014_10_21_', [4:7]};  % <file name prefix, suffix>
params.groups = {'control', [4,5];...
                 'nbqxGabazine', [6,7]};
params.excludeHS1 = {};
params.excludeHS2 = {}; 
params.celldepth = [nan, norm(-120 -215)];

% stuff for E/I and AMPA/NMDA ratios
% key for isolatedCurrents = {<current><group><Vhold><Erev>}
% Erev is to calculate driving force for conversion from pA to pS
params.isolatedCurrents = {'excit', 'control', -76, 17;...
                           'inhib', 'control', 17, -76;...
                           'ampa', 'control', -76, 17;...
                           'nmda', 'nbqxGabazine', 61, 17};
params.tags = {};
params.filter = 2e3;

HS1loc = [0];
HS2loc = [0];
Pialoc = [-120 -215];
params.celldepth = [norm(HS1loc-Pialoc), norm([HS2loc-Pialoc])];



%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_optoIV, @anlyMod_EIbalance, @anlyMod_NMDAR};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end



%% CH_100614_C cell 3

fin

%
% Notes
%
%%%%%%%%%%%%%%%%%
% Pretty good data. There's a possibility to directly compare Rs
% compensation vs. no compensation. There also appears to be some Ca2+
% contamination at Vhold = -20
%
%
% brain area: AL
% popAnalysis: EIAN, NMDAR
%

%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'CH_100614_C';      % The mouse's name
params.cellNum = 3;    % The neuron number that day
params.photo = 'CH_100614_C_cell3_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.files = {'2014_10_21_', [12:22]};  % <file name prefix, suffix>
params.groups = {'control', [12,13];...
                 'nbqxGabazine', [14,15];...
                 'NMDAR', [16:22]};
params.excludeHS1 = {};
params.excludeHS2 = {}; 

% stuff for E/I and AMPA/NMDA ratios
% key for isolatedCurrents = {<current><group><Vhold><Erev>}
% Erev is to calculate driving force for conversion from pA to pS
params.isolatedCurrents = {'excit', 'control', -72, 17;...
                           'inhib', 'control', 17, -72;...
                           'ampa', 'control', -72, 17;...
                           'nmda', 'nbqxGabazine', 50, 17};
params.tags = {};
params.filter = 2e3;


HS1loc = [0];
HS2loc = [0];
Pialoc = [-361 310];
params.celldepth = [norm(HS1loc-Pialoc), norm([HS2loc-Pialoc])];


%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_optoIV, @anlyMod_EIbalance, @anlyMod_NMDAR};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end




%% CH_141020_A pair 1



fin

%
% Notes
%
%%%%%%%%%%%%%%%%%
% 
%
%
% brain area: PM
%
% cell type: Both cells get a modest amount of spontaneous inhibition. Cell
% 2 might get slightly more. Cell 1 looks more PY like based on the in
% vitro cell fill. Neither cell gets much spontaneous excitation. I'm going
% to call both of these PY cells. 
%
% layer: 2/3
% popAnalysis: EIAN, NMDAR
%

%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'CH_141020_A';      % The mouse's name
params.cellNum = 1;    % The neuron number that day
params.photo = 'CH_141020_A_pair1';      % To assess where the light stimulus was, and the HOA that contains each cell
params.files = {'2014_11_03_', [0:9]};  % <file name prefix, suffix>
params.groups = {'control', [0,1];...
                 'nbqxGabazine', [2,3];...
                 'NMDAR', [2:9]};
params.excludeHS1 = {{'_0000', [3,4,8,12,19]}};
params.excludeHS2 = {{'_0000', [3,4,8,12,19]}}; 
params.celldepth = [norm([102 178]), norm([98 189])];

% stuff for E/I and AMPA/NMDA ratios
% key for isolatedCurrents = {<current><group><Vhold><Erev>}
% Erev is to calculate driving force for conversion from pA to pS
params.isolatedCurrents = {'excit', 'control', -72, 17;...
                           'inhib', 'control', 17, -72;...
                           'ampa', 'control', -72, 17;...
                           'nmda', 'nbqxGabazine', 50, 17};
params.tags = {};
params.filter = 2e3;


%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_optoIV, @anlyMod_EIbalance, @anlyMod_NMDAR};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end



%% CH_141020_A cell 2



fin

%
% Notes
%
%%%%%%%%%%%%%%%%%
% 
% ok. data. NMDAR only
%
% brain area: likely AL, but need to consult histology
%
% cell type: und. no control data, so no spontaneous activity. Need to
% consult the confocal imaging
%
% layer: 2/3
% popAnalysis: NMDAR
%

%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'CH_141020_A';      % The mouse's name
params.cellNum = 2;    % The neuron number that day
params.photo = 'CH_141020_A_cell2';      % To assess where the light stimulus was, and the HOA that contains each cell
params.files = {'2014_11_03_', [11:17]};  % <file name prefix, suffix>
params.groups = {'NMDAR', [11:17]};
params.excludeHS1 = {};
params.excludeHS2 = {{'_0011', [1:4]}}; 
params.celldepth = [norm([124 220]), norm([124 220])];

% stuff for E/I and AMPA/NMDA ratios
% key for isolatedCurrents = {<current><group><Vhold><Erev>}
% Erev is to calculate driving force for conversion from pA to pS
params.isolatedCurrents = {};
params.tags = {};
params.filter = 2e3;


%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_optoIV, @anlyMod_NMDAR};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end



%% CH_141020_B pair 1 



fin

%
% Notes
%
%%%%%%%%%%%%%%%%%
% 
% Pretty decent data from cell2, but cell 1 is junk.
%
% brain area: pm
%
% cell type: Cell 2 is PY cell based off of in vitro cell fill. It gets
% some spontaneous excitation, but lots of spontaneous inhibition
%
% layer: 2/3
% popAnalysis:
%

%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'CH_141020_B';      % The mouse's name
params.cellNum = 1;    % The neuron number that day
params.photo = 'CH_141020_B_pair1';      % To assess where the light stimulus was, and the HOA that contains each cell
params.files = {'2014_11_04_', [2:15]};  % <file name prefix, suffix>
params.groups = {'control', [5,6];...
                 'nbqxGabazine', [7,8];...
                 'NMDAR', [9:15]};
params.excludeHS1 = {};
params.excludeHS2 = {{'_0003', [2]}, {'_0006', [2]}, {'_0010', [2]}, {'_0013', [2]}}; 
params.celldepth = [norm([nan]), norm([-15 -178])];

% stuff for E/I and AMPA/NMDA ratios
% key for isolatedCurrents = {<current><group><Vhold><Erev>}
% Erev is to calculate driving force for conversion from pA to pS
params.isolatedCurrents = {'excit', 'control', -72, 17;...
                           'inhib', 'control', 17, -72;...
                           'ampa', 'control', -72, 17;...
                           'nmda', 'nbqxGabazine', 50, 17};
params.tags = {};
params.filter = 2e3;

% params.stimLoc = [-57 75;...
%                   -49 -8];

%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_optoIV, @anlyMod_EIbalance, @anlyMod_NMDAR};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end


%% CH_141020_C  pair 1



fin

%
% Notes
%
%%%%%%%%%%%%%%%%%
% 
% ok data from channel 1, but ch2 cuts out early.
%
% brain area: likely AL, but need to check the histology
%
% cell type: both cells have spontaneous inhibition (particularly cell 1).
% I think that this is also a PY cell based off of the in vitro image of
% NB488
%
% layer: 2/3
% popAnalysis: EIAN, NMDAR
%

%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'CH_141020_C';      % The mouse's name
params.cellNum = 1;    % The neuron number that day
params.photo = 'CH_141020_C_pair1';      % To assess where the light stimulus was, and the HOA that contains each cell
params.files = {'2014_11_06_', [0:11]};  % <file name prefix, suffix>
params.groups = {'control', [0,3];...
                 'nbqxGabazine', [4,7];...
                 'NMDAR', [4,6:11]};
params.excludeHS1 = {};
params.excludeHS2 = {}; 


% stuff for E/I and AMPA/NMDA ratios
% key for isolatedCurrents = {<current><group><Vhold><Erev>}
% Erev is to calculate driving force for conversion from pA to pS
params.isolatedCurrents = {'excit', 'control', -72, 17;...
                           'inhib', 'control', 17, -72;...
                           'ampa', 'control', -72, 17;...
                           'nmda', 'nbqxGabazine', 50, 17};
params.tags = {};
params.filter = 2e3;


HS1loc = [-60 -16];
HS2loc = [-10 0];
Pialoc1 = [-193 64];
Pialoc2 = [-137 122];
params.celldepth = [norm(HS1loc-Pialoc1), norm([HS2loc-Pialoc2])];


%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_optoIV, @anlyMod_EIbalance, @anlyMod_NMDAR};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end



%% CH_141020_D pair 1 



fin

%
% Notes
%
%%%%%%%%%%%%%%%%%
% 
% pretty good data from CH2. 
%
% brain area: PM
%
% cell type: cel 2 receives a small amt of spontaneous inhibition, but lots
% of spontaneous excitation.
%
% layer: 2/3
%
% popAnalysis: EIAN, NMDAR
%

%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'CH_141020_D';      % The mouse's name
params.cellNum = 1;    % The neuron number that day
params.photo = 'CH_141020_D_pair1';      % To assess where the light stimulus was, and the HOA that contains each cell
params.files = {'2014_11_05_', [0,3:12]};  % <file name prefix, suffix>
params.groups = {'control', [0,3];...
                 'nbqxGabazine', [4,5];...
                 'NMDAR', [6:12]};
params.excludeHS1 = {};
params.excludeHS2 = {}; 


% stuff for E/I and AMPA/NMDA ratios
% key for isolatedCurrents = {<current><group><Vhold><Erev>}
% Erev is to calculate driving force for conversion from pA to pS
params.isolatedCurrents = {'excit', 'control', -72, 17;...
                           'inhib', 'control', 17, -72;...
                           'ampa', 'control', -72, 17;...
                           'nmda', 'nbqxGabazine', 50, 17};
params.tags = {};
params.filter = 2e3;

HS1loc = [73 67];
HS2loc = [20 10];
Pialoc1 = [-125 -200];
Pialoc2 = [-125 -200];
params.celldepth = [norm(HS1loc-Pialoc1), norm([HS2loc-Pialoc2])];


%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_optoIV, @anlyMod_EIbalance, @anlyMod_NMDAR};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end




%% CH_141124_C cell 1

fin

%
% Notes
%
%%%%%%%%%%%%%%%%%
% 
% A single cell in PM. This was the first neuron tested using the
% interleaved stimulus protocol generated in matlab and then imported into
% clampex
%
% Brain area: PM
% Population analysis: EI_IO_singlePulse


%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'CH_141124_C';      % The mouse's name
params.cellNum = 1;    % The neuron number that day
params.photo = 'CH_141124_C_cell1';      % To assess where the light stimulus was, and the HOA that contains each cell
params.files = {'2014_12_11_', [4:7]};  % <file name prefix, suffix>
params.groups = {'pulses', [4:7]};
params.excludeHS1 = {{'_0004', [1:3]}};
params.excludeHS2 = {{'_0005', [4,5,11,18,21,24,25,29,31]}, {'_0006', [7:11]}}; 


% stuff for E/I and AMPA/NMDA ratios
% key for isolatedCurrents = {<current><group><Vhold><Erev>}
% Erev is to calculate driving force for conversion from pA to pS
params.isolatedCurrents = {'excit', 'control', -72, 17;...
                           'inhib', 'control', 17, -72};
params.tags = {};
params.filter = 2e3;

HS1loc = [];
HS2loc = [];
Pialoc1 = [];
Pialoc2 = [];
params.celldepth = [norm(HS1loc-Pialoc1), norm([HS2loc-Pialoc2])];


%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_avgOuterleave, @anlyMod_EI_IO};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end



%% CH_141215_A cell 2

fin

%
% Notes
%
%%%%%%%%%%%%%%%%%
% 
% I think this is an RS cell, and my notes say that the NB488 fill showed a
% clear apical dendrite. The recording location was PM. The inhibition
% recorded for the single pulses is way too fast. I don't know why this is
% the case and I don't know how to account for this...
%
% Brain area: PM

%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'CH_141215_A';      % The mouse's name
params.cellNum = 2;    % The neuron number that day
params.photo = 'CH_141215_A_cell2';      % To assess where the light stimulus was, and the HOA that contains each cell
params.files = {'2015_01_10_', [6:11]};  % <file name prefix, suffix>
params.groups = {'trains', [9,11];...
                 'pulses', [6,7,8,10]};
params.excludeHS1 = {};
params.excludeHS2 = {{'_0007', [13]}, {'_0011', [41]}}; 


% stuff for E/I and AMPA/NMDA ratios
% key for isolatedCurrents = {<current><group><Vhold><Erev>}
% Erev is to calculate driving force for conversion from pA to pS
params.isolatedCurrents = {'excit', 'control', -72, 17;...
                           'inhib', 'control', 17, -72};
params.tags = {};
params.filter = 2e3;

HS1loc = [];
HS2loc = [0 0];
Pialoc1 = [];
Pialoc2 = [220 -198];
params.celldepth = [norm(HS1loc-Pialoc1), norm([HS2loc-Pialoc2])];


%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_avgOuterleave, @anlyMod_EI_IO};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end



%% CH_141215_C cell 3

fin

%
% Notes
%
%%%%%%%%%%%%%%%%%
% some interesting data recorded from a L4 (ish) SOM+ cell in a mouse
% injected with hSyn.oChIEF.citrine. there is substantial facilitation of
% the excitatory currents, and depression of the inhibitory currents. this
% is particularly promenent at 40 Hz.
% 

%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'CH_141215_C';      % The mouse's name
params.cellNum = 3;    % The neuron number that day
params.photo = 'CH_141215_C_cell3';      % To assess where the light stimulus was, and the HOA that contains each cell
params.files = {'2014_12_29_', [12,13,14,17,18,19]};  % <file name prefix, suffix>
params.groups = {'trains', [12,13,14,17,18,19]};
params.excludeHS1 = {};
params.excludeHS2 = {}; 


% stuff for E/I and AMPA/NMDA ratios
% key for isolatedCurrents = {<current><group><Vhold><Erev>}
% Erev is to calculate driving force for conversion from pA to pS
params.isolatedCurrents = {'excit', 'control', -72, 17;...
                           'inhib', 'control', 17, -72};
params.tags = {};
params.filter = 2e3;

HS1loc = [];
HS2loc = [];
Pialoc1 = [];
Pialoc2 = [];
params.celldepth = [norm(HS1loc-Pialoc1), norm([HS2loc-Pialoc2])];


%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_avgOuterleave, @anlyMod_EI_IO};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end


%% CH_150422_D site 3

fin

%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% LFP measurements. Only including here to visualize the stimulation sites
% relative to the axon cloud. to see the analysis data, look at the .m file
% in this mouse's data directory.


%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'CH_150422_D';      % The mouse's name
params.cellNum = 3;    % The neuron number that day
params.DCsteps = [];    % DC steps for Rin and cell identification
params.photo = 'CH_150422_D_site3';      % To assess where the light stimulus was, and the HOA that contains each cell
params.files = {};      % File names of the raw data. <Nx1> cell array
params.skipSweeps = {}; % In case I need to ignore certain sweeps
params.stimLoc = [0     0
                    29  -104
                   -61  -152
                   -25  -196
                   -84  -226];    % The (x,y) coordinates of the objective at the locations stimulated with the LED. <Nx2> matrix

%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {};
    params = invitroAnalysisOverview(params);
end




%% BOOKMARK FOR CH MICE
% data added to population analyses up to here

% Need to add

% CH_091114_F (ChR2, good data)
% CH_100614_A (EIAN) ## make sure this mouse has been added
% CH_100614_B (EIAN) Need to check histology
% CH_100614_C (EIAN) Need to check histology
% CH_141020_A Confirm brain region on cell 2

% CH_141201_A & _B Need to add to the sheet. Practice data file for subtracting
% the excitation from mixed currents

% CH_141215_A to _F Add to sheet, make sure that pop_anly workbooks are
% updated.




%% EB_031014_A cell 2

fin

%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shortish inter-sweep intervals. If I was interested in run down from
% sweep to sweep, this might be a good cell to look at.


%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'EB_031014_A';      % The mouse's name
params.cellNum = 2;    % The neuron number that day
params.tags = {'Rundown, APV'};



%% EB_031014_B Cell 1

fin

%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Looking at the latency to the PSC when the LED targets distal locations.
% In this experiment, I didn't do a great job of targeting ChR2+ axons, but
% there are a few stimulation places that evoke good responses at a longer
% latency than at the soma.
%
% This appears to be a L2/3 PY cell. No histology was saved, so I don't
% know which HOA this is. 
%
% 


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
params.tags = {'distal stim, latency, room temp'};
params.subtractexp = false;
params.filter = 4000;


%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

%% EB_031014_B Cell 2
 
fin

%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Looking at the latency to the PSC when the LED targets distal locations.
% There's definitely an effect of stimulus location, but the Ra also goes
% down a bit....
%
% This appears to be a L2/3 FS cell. No histology was saved, so I don't
% know which HOA this is. I think it's medial to V1.


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
params.tags = {'distal stim, latency, room temp'};
params.subtractexp = false;
params.filter = 4000;

%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end



%% EB_031014_D Cell 1

fin

%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Access is pretty bad for this cell (stable, but bad: 50 MOhms).
% Nonetheless, I tested a few different spatial locations and TFs. There
% was depression at all TFs. This neuron might be a good example of "false
% depression" at the soma. 
%
% This cell is an FS cell in L2/3 of an unknown HOA. No histology, but the
% note on Docubase suggest that it's PM. There is good thalamic signal, and
% the axon field is medial to V1.
%
% I added this cell to the distal stim/latency pop anly.



%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'EB_031014_D';      % The mouse's name
params.cellNum = 1;    % The neuron number that day
params.DCsteps = '2014_04_01_0000';    % DC steps for Rin and cell identification
params.photo = 'cell_1_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.validCh = 'HS2_';    % 'HS2_' or 'HS1_'
params.files = {'2014_04_01_0002';...
                '2014_04_01_0005';...
                '2014_04_01_0007';...
                '2014_04_01_0009';...
                '2014_04_01_0010';...
                '2014_04_01_0011'};      % File names of the raw data. <Nx1> cell array
params.skipSweeps = {}; % In case I need to ignore certain sweeps
params.vHold = -85.*(ones(size(params.files)));     % The holding potential for vClamp experiments. One value for each expt. <Nx1>
params.stimLoc = [0 0;...
                 -12 109;...
                 -178 286;...
                 -504 257;...
                 -504 257;...
                 0 0];    % The (x,y) coordinates of the objective at the locations stimulated with the LED. <Nx2> matrix
params.legTxt = {'Soma, 20 Hz';...
                 'Pos 2, , 20 Hz';...
                 'Pos 3, 20 Hz';...
                 'Pos 4, 20 Hz';...
                 'Pos 4, 5 Hz';...
                 'Soma, 5 Hz'};     % Text that will appear in figures to annotate each data file. <Nx1> cell array
params.tags = {'distal stim, latency, room temp, TF, false depression'};
params.subtractexp = false;
params.filter = 4000;

%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end



%% EB_031014_D Cell 2

fin

%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This neuron was tested at 20Hz at a variety of different locations and
% voltages to the LED. There is prominant "false facilitation", but at the
% distal LED locations, the trend (synaptic depression) is consistent. 
%
% This is a L2/3 PY cell. Likely in PM, but unknown.
%



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
params.legTxt = {'Soma';...
                 'Pos 1';...
                 'Pos 2 4.5 volts';...
                 'Pos 2 3 volts';...
                 'Pos 2 10 volts';...
                 'Pos 3 ';...
                 'Soma'};     % Text that will appear in figures to annotate each data file. <Nx1> cell array
params.tags = {'distal stim, latency, room temp, false depression'};
params.subtractexp = false;
params.filter = 4000;

%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end

%% EB_031014_E Cell 3

fin

%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Looking at the latency to response with stimulation at different cortical
% locations. This file uses different pulse widths...
%
% A L2/3 PY cell from an area lateral to V1.



%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'EB_031014_E';      % The mouse's name
params.cellNum = 3;    % The neuron number that day
params.DCsteps = '2014_03_25_0002';    % DC steps for Rin and cell identification
params.photo = 'cell_3_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.validCh = 'HS2_';    % 'HS2_' or 'HS1_'
params.files = {'2014_03_25_0007';...
         '2014_03_25_0008';...
         '2014_03_25_0009';...
         '2014_03_25_0010';...
         '2014_03_25_0011';...
         '2014_03_25_0012';...
         '2014_03_25_0014';...
         '2014_03_25_0015'};      % File names of the raw data. <Nx1> cell array
params.skipSweeps = {}; % In case I need to ignore certain sweeps
params.vHold = -85.*ones(size(params.files));     % The holding potential for vClamp experiments. One value for each expt. <Nx1>
params.stimLoc = [-257 201;...
         -357 291;...
        -678 384;...
         -833 48;...
         -302 29;...
        -302 29;...
        -302 29;...
        -302 29];    % The (x,y) coordinates of the objective at the locations stimulated with the LED. <Nx2> matrix
params.legTxt = {'Cortex 1';...
         'Cortex 2';...
         'Axons 1';...
         'Axons 2';...
         'Cortex 3 Slow';...
         'Cortex 3 Fast';...
         'Cortex 3 Half Moon';...
         'Cortex 3 No Half Moon'};     % Text that will appear in figures to annotate each data file. <Nx1> cell array
params.tags = {'distal stim, latency, room temp'};
params.subtractexp = false;
params.filter = 4000;

%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end

%% EB_031714_A Cell 1

fin

%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This cell shows strong false depression for 20 Hz LED pulses delivered to
% the soma, but not for 5 or 10 Hz. The more distal locations show less of
% the false depression, and even a very small amount of facilitation at 20
% Hz. RECORDED AT ROOM TEMP!!!
%
% This is a L2/3 PY cell. Possibly from LM, AL, or V1. I think that an
% argument could be made for this cell being in AL. There is good thalamic
% expression, and I don't think that it's in the edge of the injection
% (Although I think that V1 is labeled in this slice).



%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'EB_031714_A';      % The mouse's name
params.cellNum = 1;    % The neuron number that day
params.DCsteps = '2014_04_02_0000';    % DC steps for Rin and cell identification
params.photo = 'cell_1_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.validCh = 'HS2_';    % 'HS2_' or 'HS1_'
params.files = {'2014_04_02_0002';...
                '2014_04_02_0003';...
                '2014_04_02_0004';...
                '2014_04_02_0006';...
                '2014_04_02_0007';...
                '2014_04_02_0008';...
                '2014_04_02_0011';...
                '2014_04_02_0012';...
                '2014_04_02_0013';...
                '2014_04_02_0014'};      % File names of the raw data. <Nx1> cell array
params.skipSweeps = {}; % In case I need to ignore certain sweeps
params.vHold = -85.*ones(size(params.files));     % The holding potential for vClamp experiments. One value for each expt. <Nx1>
params.stimLoc = [0 0;...
                  0 0;...
                  0 0;...
                  -255 63;...
                  -255 63;...
                  -255 63;...
                  -517 20;...
                  -517 20;...
                  0 0;...
                  0 0];    % The (x,y) coordinates of the objective at the locations stimulated with the LED. <Nx2> matrix
params.legTxt = {'Soma TF = 5',...
                 'Soma TF = 20',...
                 'Soma TF = 10',...
                 'Pos 1, TF = 20',...
                 'Pos 1, TF = 5',...
                 'Pos 1, TF = 10',...
                 'Pos 2, TF = 20',...
                 'Pos 2, TF = 5',...
                 'Soma TF = 5',...
                 'Soma TF = 20'};     % Text that will appear in figures to annotate each data file. <Nx1> cell array
params.tags = {'False depression, TF, trains, HOA'};
params.subtractexp = false;
params.filter = 4000;

%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end


%% EB_031714_A Cell 3

fin

%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Junky data from a cell in an unknown area lateral to V1. I didn't
% stimuliate in the correct area, and the cell was at the edge of an axon
% field. I should do a better job of targeting cells in the center of the
% HOA.
%



%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'EB_031714_A';      % The mouse's name
params.cellNum = 3;    % The neuron number that day
params.DCsteps = '2014_04_02_0017';    % DC steps for Rin and cell identification
params.photo = 'cell_3_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.validCh = 'HS2_';    % 'HS2_' or 'HS1_'
params.files = {'2014_04_02_0018',...
                '2014_04_02_0022',...
                '2014_04_02_0023',...
                '2014_04_02_0024',...
                '2014_04_02_0025',...
                '2014_04_02_0026',...
                '2014_04_02_0027'};      % File names of the raw data. <Nx1> cell array
params.skipSweeps{7} = [9,10]; % In case I need to ignore certain sweeps
params.vHold = -85.*ones(size(params.files));     % The holding potential for vClamp experiments. One value for each expt. <Nx1>
params.stimLoc = [0 0;...
                  -349 117;...
                  -327 -87;...
                  -219 -229;...
                  -31 335;...
                  135 -475;...
                  404 -170];    % The (x,y) coordinates of the objective at the locations stimulated with the LED. <Nx2> matrix
params.legTxt = {};     % Text that will appear in figures to annotate each data file. <Nx1> cell array
params.tags = {};
params.subtractexp = false;
params.filter = 4000;

%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end


%% EB_031714_C Cell 1

fin

%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Not much data from this cell. Just 5 and 20 Hz trains at the soma.
% There's massive depression at 20 Hz, but I can't confirm that it's "false
% depression" b/c I didn't test distal LED stim locations. Very strong
% depression...
%
% This cell was a crummy leaky PY cell in lower L2/3. Could be in AL or LM.


%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'EB_031714_C';      % The mouse's name
params.cellNum = 1;    % The neuron number that day
params.DCsteps = '2014_04_05_0006';    % DC steps for Rin and cell identification
params.photo = 'cell1_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.validCh = 'HS2_';    % 'HS2_' or 'HS1_'
params.files = {'2014_04_05_0008', '2014_04_05_0009'};      % File names of the raw data. <Nx1> cell array
params.skipSweeps = {}; % In case I need to ignore certain sweeps
params.vHold = [-85 -85];     % The holding potential for vClamp experiments. One value for each expt. <Nx1>
params.stimLoc = [0 0; 0 0];    % The (x,y) coordinates of the objective at the locations stimulated with the LED. <Nx2> matrix
params.legTxt = {'20', '5'};     % Text that will appear in figures to annotate each data file. <Nx1> cell array
params.tags = {'False depression'};
params.subtractexp = false;
params.filter = 4000;

%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end


%% EB_031714_C Cell 2

fin

%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Only had a chance to do trains at the soma with FS open. There's
% facilitation for 20 Hz, but depression for 5 and 40 Hz.
%
% This is an FS cell from L2/3 in what is likely AL. The axon field is
% quite anterior. The callosum has fused.

fin
%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'EB_031714_C';      % The mouse's name
params.cellNum = 2;    % The neuron number that day
params.DCsteps = '2014_04_05_0014';    % DC steps for Rin and cell identification
params.photo = 'cell2_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.validCh = 'HS2_';    % 'HS2_' or 'HS1_'
params.files = {'2014_04_05_0017','2014_04_05_0018','2014_04_05_0019'};      % File names of the raw data. <Nx1> cell array
params.skipSweeps = {}; % In case I need to ignore certain sweeps
params.vHold = [-85 -85 -85];     % The holding potential for vClamp experiments. One value for each expt. <Nx1>
params.stimLoc = [0 0; 0 0; 0 0];    % The (x,y) coordinates of the objective at the locations stimulated with the LED. <Nx2> matrix
params.legTxt = {'20', '5', '40'};     % Text that will appear in figures to annotate each data file. <Nx1> cell array
params.tags = {'HOA, TF'};
params.subtractexp = false;
params.filter = 4000;

%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end

%% EB_031714_D Cell 1

fin

%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% LED light trains at a few locations in the cortex. Strong depression at
% the soma, less so at the other locaitons. I washed in TTX on this cell,
% which has the effect of suppressing evoked PSC when the light is targeted
% to any location (including the Soma)
%
% This is a funky L5A PY cell. I think it's from Area LM. The first slice
% of the series has P and LM, so I think that the slice I recorded from is
% LM. 
%


%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'EB_031714_D';      % The mouse's name
params.cellNum = 1;    % The neuron number that day
params.DCsteps = '2014_04_04_0001';    % DC steps for Rin and cell identification
params.photo = 'cell1_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.validCh = 'HS2_';    % 'HS2_' or 'HS1_'
params.files = {'2014_04_04_0004';...
                '2014_04_04_0005';...
                '2014_04_04_0007';...
                '2014_04_04_0008';...
                '2014_04_04_0010';...
                '2014_04_04_0011'};      % File names of the raw data. <Nx1> cell array
params.skipSweeps = {[], [10], [], [], [], [1]}; % In case I need to ignore certain sweeps
params.vHold = -85.*ones(size(params.files));     % The holding potential for vClamp experiments. One value for each expt. <Nx1>
params.stimLoc = [0 0;...
                  0 0;...
                  -396 57;...
                  -396 57;...
                  -232 175;...
                  -232 175];    % The (x,y) coordinates of the objective at the locations stimulated with the LED. <Nx2> matrix
params.legTxt = {'Soma 20', 'Soma 5', 'Pos 1 20 Hz', 'Pos 1: 5 Hz', 'Pos 2: 20 Hz', 'Pos 2: 5 Hz'};     % Text that will appear in figures to annotate each data file. <Nx1> cell array
params.tags = {'TTX, false depression, HOA, TF'};
params.subtractexp = false;
params.filter = 4000;

%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end

%% EB_031714_D Cell 3

fin

%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stimulated this cell at one location distal from the soma. There was a
% modest facilitation effect when the light targeted the soma for 20 and 40
% Hz stimuli. Depression for 5 hz was profound when the light targeted a
% different location.
%
% This is a L2/3 YP cell. It could be in anterior LM (weak signal in the
% thalamus). Another possibility is that it is in LI and that the area
% medial to the axon field I recorded is AL. LI might be more likly unless
% LM has a lateral overlap with AL in the A/P direction.
%



%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'EB_031714_D';      % The mouse's name
params.cellNum = 3;    % The neuron number that day
params.DCsteps = '2014_04_04_0023';    % DC steps for Rin and cell identification
params.photo = 'cell3_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.validCh = 'HS2_';    % 'HS2_' or 'HS1_'
params.files = {'2014_04_04_0025';...
                '2014_04_04_0026';...
                '2014_04_04_0027';...
                '2014_04_04_0030';...
                '2014_04_04_0031';...
                '2014_04_04_0032'};      % File names of the raw data. <Nx1> cell array
params.skipSweeps = {}; % In case I need to ignore certain sweeps
params.vHold = -85.*ones(size(params.files));     % The holding potential for vClamp experiments. One value for each expt. <Nx1>
params.stimLoc = [0 0;...
                  0 0;...
                  0 0;...
                  -16 301;...
                  -16 301;...
                  -16 301];    % The (x,y) coordinates of the objective at the locations stimulated with the LED. <Nx2> matrix
params.legTxt = {'Soma: TF = 20', 'Soma: TF = 5', 'Soma: TF = 40', 'Pos1: TF = 20', 'Pos1: TF = 5', 'Pos1: TF = 40'};     % Text that will appear in figures to annotate each data file. <Nx1> cell array
params.tags = {'TTX'};
params.subtractexp = false;
params.filter = 4000;

%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end


%% EB_040714_A cell 1

fin

%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A PY cell from PM it think. Strong depression for 5 Hz stim. More so
% than for anything else.
% adding this to the population analysis list.



%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.debug = false;
params.mouse = 'EB_040714_A';      % The mouse's name
params.cellNum = 1;    % The neuron number that day
params.DCsteps = '2014_04_21_0001';    % DC steps for Rin and cell identification
params.photo = 'cell1_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.validCh = 'HS2_';    % 'HS2_' or 'HS1_'
params.files = {'2014_04_21_0004';...
                '2014_04_21_0010';...
                '2014_04_21_0011';...
                '2014_04_21_0012'};      % File names of the raw data. <Nx1> cell array
params.skipSweeps = {[], [21]}; % In case I need to ignore certain sweeps
params.vHold = -85.*ones(size(params.files));     % The holding potential for vClamp experiments. One value for each expt. <Nx1>
params.stimLoc = [0 0;...
                  -153, 152;...
                  -153, 152;...
                  -153, 152];    % The (x,y) coordinates of the objective at the locations stimulated with the LED. <Nx2> matrix
params.legTxt = {'soma 20', 'pos 1 5 Hz', 'pos 1 20 Hz', 'pos 1 40 Hz'};     % Text that will appear in figures to annotate each data file. <Nx1> cell array
params.tags = {};
params.subtractexp = false;
params.filter = 4000;

%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end

%% EB_040714_A cell 2

fin

%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Pretty crummy data. An FS cell from AL or LM, but there was limited TF
% sampling (only 5 and 20) and the neuron required a lot of holding
% current... At the soma there is some massize false facilitation?

%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'EB_040714_A';      % The mouse's name
params.cellNum = 2;    % The neuron number that day
params.DCsteps = '2014_04_21_0016';    % DC steps for Rin and cell identification
params.photo = 'cell2_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.validCh = 'HS2_';    % 'HS2_' or 'HS1_'
params.files = {'2014_04_21_0020';...
                '2014_04_21_0023';...
                '2014_04_21_0024'};      % File names of the raw data. <Nx1> cell array
params.skipSweeps = {[4,9]}; % In case I need to ignore certain sweeps
params.vHold = -85.*ones(1,3);     % The holding potential for vClamp experiments. One value for each expt. <Nx1>
params.stimLoc = [0 0; -275 -45; -275 -45;];    % The (x,y) coordinates of the objective at the locations stimulated with the LED. <Nx2> matrix
params.legTxt = {'soma 20 Hz', '20 Hz', '5 Hz'};     % Text that will appear in figures to annotate each data file. <Nx1> cell array
params.tags = {};
params.subtractexp = false;
params.filter = 4000;

%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end
%% EB_040714_A cell 3

fin

%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% A L2/3 PY cell. Some Ih present. Pretty good data, adding to population
% anlysis text files. I think that this cell is in AL.



%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'EB_040714_A';      % The mouse's name
params.cellNum = 3;    % The neuron number that day
params.DCsteps = '2014_04_21_0027';    % DC steps for Rin and cell identification
params.photo = 'cell3_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.validCh = 'HS2_';    % 'HS2_' or 'HS1_'
params.files = {'2014_04_21_0029';...
                '2014_04_21_0031';...
                '2014_04_21_0032';...
                '2014_04_21_0033'};      % File names of the raw data. <Nx1> cell array
params.skipSweeps = {}; % In case I need to ignore certain sweeps
params.vHold = -85.*ones(1,4);     % The holding potential for vClamp experiments. One value for each expt. <Nx1>
params.stimLoc = [0 0; -300 -50; -300 -50; -300 -50;];    % The (x,y) coordinates of the objective at the locations stimulated with the LED. <Nx2> matrix
params.legTxt = {'soma 20 hz', '20', '40', '5'};     % Text that will appear in figures to annotate each data file. <Nx1> cell array
params.tags = {};
params.subtractexp = false;
params.filter = 4000;

%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end
%% EB_040714_B cell 1
fin

%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is a cool data set. There is very good agreement between the two
% distal stimulation sites, and they are both different than recordings at
% the soma. There is also very good agreement with the soma-targeted files
% even though they're at the beginning and end of the experiment.
%
% This cell is likely in PM, but the histology can't confirm this b/c the
% cortex was damaged and became detached (and was not re-attached during
% histology). The picture taken under the slice scope, and notes in the
% text suggest that the cell is from PM, so I'm tentatively putting it in
% that group.
%

%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'EB_040714_B';      % The mouse's name
params.cellNum = 1;    % The neuron number that day
params.DCsteps = '2014_04_23_0001';    % DC steps for Rin and cell identification
params.photo = 'cell_1_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.validCh = 'HS2_';    % 'HS2_' or 'HS1_'
params.files = {'2014_04_23_0006';...
                     '2014_04_23_0012';...
                     '2014_04_23_0013';...
                     '2014_04_23_0014';...
                     '2014_04_23_0017';...
                     '2014_04_23_0018';...
                     '2014_04_23_0019';...
                     '2014_04_23_0020';...
                     '2014_04_23_0021';...
                     '2014_04_23_0022'};      % File names of the raw data. <Nx1> cell array
params.skipSweeps = {}; % In case I need to ignore certain sweeps
params.vHold = -85.*ones(size(params.files));     % The holding potential for vClamp experiments. One value for each expt. <Nx1>
params.stimLoc = [0 0;...
                   -312 4;...
                   -312 4;...
                   -312 4;...
                   -430 13;...
                   -430 13;...
                   -430 13;...
                   0 0;...
                   0 0;...
                   0 0];    % The (x,y) coordinates of the objective at the locations stimulated with the LED. <Nx2> matrix
params.legTxt = {'1 soma 20', '1 40', '1 20', '1 5', '2 20', '2 40', '2 5', '2 soma 20', '2 soma 40', '2 soma 5'};     % Text that will appear in figures to annotate each data file. <Nx1> cell array
params.tags = {};
params.subtractexp = false;
params.filter = 4000;

%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end
%%  EB_040714_B cell 3

fin

%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% L2/3 PY cell. PPfacil for 20 at some distal locations, but not at all of
% them. Data are subtally different than at the soma. This cell is from PM
% probably. adding to the population analysis text file.


%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'EB_040714_B';      % The mouse's name
params.cellNum = 3;    % The neuron number that day
params.DCsteps = '2014_04_23_0029';    % DC steps for Rin and cell identification
params.photo = 'cell_3_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.validCh = 'HS2_';    % 'HS2_' or 'HS1_'
params.files = {'2014_04_23_0032';...
                '2014_04_23_0034';...
                '2014_04_23_0035';...
                '2014_04_23_0036';...
                '2014_04_23_0037';...
                '2014_04_23_0041';...
                '2014_04_23_0042';...
                '2014_04_23_0043';...
                '2014_04_23_0044';...
                '2014_04_23_0045'};      % File names of the raw data. <Nx1> cell array
params.skipSweeps = {[22:25],[],[],[],[],[],[],[6],[]}; % In case I need to ignore certain sweeps
params.vHold = -85.*ones(size(params.files));     % The holding potential for vClamp experiments. One value for each expt. <Nx1>
params.stimLoc = [0 0;...
                 -297 210;...
                 -297 210;...
                 -297 210;...
                 -297 210;...
                 -400 235;...
                 -400 235;...
                 -400 235;...
                 -400 235;...
                 0 0];    % The (x,y) coordinates of the objective at the locations stimulated with the LED. <Nx2> matrix
params.legTxt = {'soma 20 Hz';...
                'pos 1 20';...
                'pos 1 40';...
                'pos 1 5';...
                'pos 1 10';...
                'pos 2 20';...
                'pos 2 40';...
                'pos 2 5';...
                'pos 2 10';...
                'soma 20 Hz'};     % Text that will appear in figures to annotate each data file. <Nx1> cell array
params.tags = {};
params.subtractexp = false;
params.filter = 4000;

%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end

%% EB_040714_C cell 1

fin

%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% L2/3 PY cell. Data are faily clean.  I think this is from area PM. Adding
% it to the list for population analysis.
%




%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'EB_040714_C';      % The mouse's name
params.cellNum = 1;    % The neuron number that day
params.DCsteps = '2014_04_25_0001';    % DC steps for Rin and cell identification
params.photo = 'cell_1_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.validCh = 'HS2_';    % 'HS2_' or 'HS1_'
params.files = {'2014_04_25_0004';...
                '2014_04_25_0006';...
                '2014_04_25_0007';...
                '2014_04_25_0011';...
                '2014_04_25_0012';...
                '2014_04_25_0013'};      % File names of the raw data. <Nx1> cell array
params.skipSweeps = {}; % In case I need to ignore certain sweeps
params.vHold = -85 .* ones(size(params.files));     % The holding potential for vClamp experiments. One value for each expt. <Nx1>
params.stimLoc = [0 0;...
                -235 5;...
                -235 5';...
                -211 181;...
                -211 181;...
                -211 181];    % The (x,y) coordinates of the objective at the locations stimulated with the LED. <Nx2> matrix
params.legTxt = {'soma';...
                'pos 1 20 Hz';...
                'pos 1 40 Hz';...
                'pos 2 5 Hz';...
                'pos 2 40 Hz';...
                'pos 2 20 Hz'};     % Text that will appear in figures to annotate each data file. <Nx1> cell array
params.tags = {};
params.subtractexp = false;
params.filter = 4000;


%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end
%% EB_040714_C cell 2

fin

%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% A small FS cell from L2/3 in PM. Good agreement between the stimulation
% locations. I wish there were more data. Added to pop analysis text file.




%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'EB_040714_C';      % The mouse's name
params.cellNum = 2;    % The neuron number that day
params.DCsteps = '2014_04_25_0015';    % DC steps for Rin and cell identification
params.photo = 'cell_2_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.validCh = 'HS2_';    % 'HS2_' or 'HS1_'
params.files = {'2014_04_25_0019';...
                '2014_04_25_0025';...
                '2014_04_25_0026'};      % File names of the raw data. <Nx1> cell array
params.skipSweeps = {[] [] [18:25]}; % In case I need to ignore certain sweeps
params.vHold = [-85 -85 -85];     % The holding potential for vClamp experiments. One value for each expt. <Nx1>
params.stimLoc = [0 0; -258 10; -258 10];    % The (x,y) coordinates of the objective at the locations stimulated with the LED. <Nx2> matrix
params.legTxt = {'Soma' 'pos 1 40 Hz', 'pos1 20 Hz'};     % Text that will appear in figures to annotate each data file. <Nx1> cell array
params.tags = {};
params.subtractexp = false;
params.filter = 4000;


%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end
%% EB_040714_C cell 3

fin

%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Acess resistance might be a problem. Other than that, pretty standard
% data. Two different locations. L2/3 PY cell in. I think that this cell
% was located in PM based off the nick in the axon bundle below L6 and the
% location of the harp string scars. 
%



%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'EB_040714_C';      % The mouse's name
params.cellNum = 3;    % The neuron number that day
params.DCsteps = '2014_04_25_0028';    % DC steps for Rin and cell identification
params.photo = 'cell_3_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.validCh = 'HS2_';    % 'HS2_' or 'HS1_'
params.files = {'2014_04_25_0030';...
                '2014_04_25_0031';...
                '2014_04_25_0033';...
                '2014_04_25_0034';...
                '2014_04_25_0035';...
                '2014_04_25_0038';...
                '2014_04_25_0039';...
                '2014_04_25_0040';...
                '2014_04_25_0041'};      % File names of the raw data. <Nx1> cell array
params.skipSweeps = {}; % In case I need to ignore certain sweeps
params.vHold = -85.*ones(size(params.files));     % The holding potential for vClamp experiments. One value for each expt. <Nx1>
params.stimLoc = [0 0;...
                0 0;...
                -261 260;...
                -261 260';...
                -261 260';...
                -366 394;...
                -366 394;...
                -366 394;...
                -366 394'];    % The (x,y) coordinates of the objective at the locations stimulated with the LED. <Nx2> matrix
params.legTxt = {'soma 20';...
                'soma 40';...
                'pos1 20';...
                'pos1 40';...
                'pos1 5';...
                'pos2 20';...
                'pos2 40';...
                'pos2 5';...
                'pos2 10'};     % Text that will appear in figures to annotate each data file. <Nx1> cell array
params.tags = {};
params.subtractexp = false;
params.filter = 4000;


%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end
%% EB_040714_D Cell 2

fin

%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Neat data. This neuron seems like an FS cell on the basis of AP width,
% and there's very little spike rate adaptation. During intense bouts of
% APs, the baseline Vm increases. There may be some sort of fast AHP, and
% there's definitely a whopping Ih current activated near -120 mV. 
%
% For all TFs tested (only at the soma unfortunately) this neuron exhibited
% PPfacilitation. Facilitation was strongest for high TFs. 



%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'EB_040714_D';      % The mouse's name
params.cellNum = 2;    % The neuron number that day
params.DCsteps = '2014_04_24_0008';    % DC steps for Rin and cell identification
params.photo = 'cell_2_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.validCh = 'HS2_';    % 'HS2_' or 'HS1_'
params.files = {'2014_04_24_0011';...
                '2014_04_24_0012';...
                '2014_04_24_0013';...
                '2014_04_24_0014'};      % File names of the raw data. <Nx1> cell array
params.skipSweeps = {}; % In case I need to ignore certain sweeps
params.vHold = -85.*ones(size(params.files));     % The holding potential for vClamp experiments. One value for each expt. <Nx1>
params.stimLoc = [0 0;...
                  0 0;...
                  0 0;...
                  0 0];    % The (x,y) coordinates of the objective at the locations stimulated with the LED. <Nx2> matrix
params.legTxt = {'soma 20 hz', 'soma 40 hz', 'soma 5 hz', 'soma 10 hz', };     % Text that will appear in figures to annotate each data file. <Nx1> cell array
params.tags = {'FS cell, facilitation, TF STP'};
params.subtractexp = false;
params.filter = 4000;


%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end



%% EB_040714_D Cell 3

fin

%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PY cell from PM



%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'EB_040714_D';      % The mouse's name
params.cellNum = 3;    % The neuron number that day
params.DCsteps = '2014_04_24_0016';    % DC steps for Rin and cell identification
params.photo = 'cell_3_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.validCh = 'HS2_';    % 'HS2_' or 'HS1_'
params.files = {'2014_04_24_0018';...
                '2014_04_24_0019';...
                '2014_04_24_0024';...
                '2014_04_24_0025';...
                '2014_04_24_0026';...
                '2014_04_24_0027';...
                '2014_04_24_0028';...
                '2014_04_24_0029';...
                '2014_04_24_0030';...
                '2014_04_24_0031';...
                '2014_04_24_0032';...
                '2014_04_24_0033';...
                '2014_04_24_0035'};      % File names of the raw data. <Nx1> cell array
params.skipSweeps = {}; % In case I need to ignore certain sweeps
params.vHold = -85.*ones(size(params.files));     % The holding potential for vClamp experiments. One value for each expt. <Nx1>
params.stimLoc = [0 0;...
                0 0;...
                -241 15;...
                -241 15;...
                -241 15;...
                -334 26;...
                -334 26;...
                -334 26;...
                -401 33;...
                -401 33;...
                -401 33;...
                -505 40;...
                0 0];    % The (x,y) coordinates of the objective at the locations stimulated with the LED. <Nx2> matrix
params.legTxt = {'Soma 20 Hz';...
                'Soma 40 Hz';...
                'Pos1 5 Hz';...
                'Pos1 20 Hz';...
                'Pos1 40 Hz';...
                'Pos2 40 Hz';...
                'Pos2 20 Hz';...
                'Pos2 5 Hz';...
                'Pos3 5 Hz';...
                'Pos3 20 Hz';...
                'Pos3 40 Hz';...
                'Pos4 20 Hz';...
                'Soma 20 Hz'};     % Text that will appear in figures to annotate each data file. <Nx1> cell array
params.tags = {};
params.subtractexp = false;
params.filter = 4000;


%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end


%% EB_041414_A cell 1


fin

%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PY cell likely from PM. Adding to population analysis.


%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'EB_041414_A';      % The mouse's name
params.cellNum = 1;    % The neuron number that day
params.DCsteps = '2014_05_02_0001';    % DC steps for Rin and cell identification
params.photo = 'cell_1_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.validCh = 'HS2_';    % 'HS2_' or 'HS1_'
params.files = {'2014_05_02_0003';...
                '2014_05_02_0005';...
                '2014_05_02_0007';...
                '2014_05_02_0008';...
                '2014_05_02_0009';...
                '2014_05_02_0010';...
                '2014_05_02_0011';...
                '2014_05_02_0012';...
                '2014_05_02_0013'};      % File names of the raw data. <Nx1> cell array
params.skipSweeps = {}; % In case I need to ignore certain sweeps
params.vHold = -85 .* ones(size(params.files));     % The holding potential for vClamp experiments. One value for each expt. <Nx1>
params.stimLoc = [0 0;...
                  -206 182;...
                  -206 182;...
                  -206 182;...
                  -275 290;...
                  -275 290;...
                  -275 290;...
                  -275 290;...
                  -358 372];    % The (x,y) coordinates of the objective at the locations stimulated with the LED. <Nx2> matrix
params.legTxt = {};     % Text that will appear in figures to annotate each data file. <Nx1> cell array
params.tags = {'axon stim', 'synaptic blockers'};
params.subtractexp = false;
params.filter = 4000;


%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end


%% EB_041414_A Cell 3

fin

%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Online notes say this was a small cell body. Definitely inhibitory
% interneuron, but I'm not sure what type. Maybe an FS cell (little spike
% frequency accomodation), a small Ih current. I'm putting it in the FS
% catagory... Access resistance is poor.


%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'EB_041414_A';      % The mouse's name
params.cellNum = 3;    % The neuron number that day
params.DCsteps = '2014_05_02_0022';    % DC steps for Rin and cell identification
params.photo = 'cell_3_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.validCh = 'HS2_';    % 'HS2_' or 'HS1_'
params.files = {'2014_05_02_0025';...
                '2014_05_02_0029';...
                '2014_05_02_0030';...
                '2014_05_02_0031'};      % File names of the raw data. <Nx1> cell array
params.skipSweeps = {}; % In case I need to ignore certain sweeps
params.vHold = -85.*ones(size(params.files));     % The holding potential for vClamp experiments. One value for each expt. <Nx1>
params.stimLoc = [0 0;...
                  -133 104;...
                  -133 104;...
                  -133 104];    % The (x,y) coordinates of the objective at the locations stimulated with the LED. <Nx2> matrix
params.legTxt = {};     % Text that will appear in figures to annotate each data file. <Nx1> cell array
params.tags = {'Axon stim', 'synaptic blockers', 'retro infect'};
params.subtractexp = false;
params.filter = 4000;


%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end


%% EB_041414_B Cell 1

fin

%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% An interneuron from a L2/3 PM. No Ih, and looks like a delayed onset.
% Adding to the population analysis for PM.


%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'EB_041414_B';      % The mouse's name
params.cellNum = 1;    % The neuron number that day
params.DCsteps = '2014_05_01_0001';    % DC steps for Rin and cell identification
params.photo = 'cell_1_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.validCh = 'HS2_';    % 'HS2_' or 'HS1_'
params.files = {'2014_05_01_0008';...
                '2014_05_01_0011';...
                '2014_05_01_0012';...
                '2014_05_01_0013';...
                '2014_05_01_0018'};      % File names of the raw data. <Nx1> cell array
params.skipSweeps = {}; % In case I need to ignore certain sweeps
params.vHold = -85.*ones(size(params.files));     % The holding potential for vClamp experiments. One value for each expt. <Nx1>
params.stimLoc = [0 0;...
                  -315 156;...
                  -315 156;...                  
                  -315 156;...
                  -341 227];    % The (x,y) coordinates of the objective at the locations stimulated with the LED. <Nx2> matrix
params.legTxt = {};     % Text that will appear in figures to annotate each data file. <Nx1> cell array
params.tags = {};
params.subtractexp = false;
params.filter = 4000;


%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end


%% EB_041414_B Cell 3

fin

%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% L2/3 PY cell. Adding to pop analysis for PM.
%

%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'EB_041414_B';      % The mouse's name
params.cellNum = 3;    % The neuron number that day
params.DCsteps = '2014_05_01_0031';    % DC steps for Rin and cell identification
params.photo = 'cell_3_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.validCh = 'HS2_';    % 'HS2_' or 'HS1_'
params.files = {'2014_05_01_0033';...
                '2014_05_01_0035';...
                '2014_05_01_0036';...
                '2014_05_01_0038';...
                '2014_05_01_0041';...
                '2014_05_01_0042';...
                '2014_05_01_0043'};      % File names of the raw data. <Nx1> cell array
params.skipSweeps = {}; % In case I need to ignore certain sweeps
params.vHold = -85.*ones(size(params.files));     % The holding potential for vClamp experiments. One value for each expt. <Nx1>
params.stimLoc = [0 0;...
                  -198 187;...
                  -198 187;...
                  -198 187;...
                  -307 187;...
                  -307 187;...
                  -307 187];    % The (x,y) coordinates of the objective at the locations stimulated with the LED. <Nx2> matrix
params.legTxt = {};     % Text that will appear in figures to annotate each data file. <Nx1> cell array
params.tags = {'Axon stim', 'retro infect'};
params.subtractexp = false;
params.filter = 4000;


%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end


%% EB_041414_C Cell 1

fin

%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Probably an FS cell. Tough to locate. Lower L2/3 upper layer 4. Also,
% could be on the border between AL and LM. Adding to the population
% analysis for AL, but will make a note that it's hard to tell.

%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'EB_041414_C';      % The mouse's name
params.cellNum = 1;    % The neuron number that day
params.DCsteps = '2014_04_28_0001';    % DC steps for Rin and cell identification
params.photo = 'cell_1_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.validCh = 'HS2_';    % 'HS2_' or 'HS1_'
params.files = {'2014_04_28_0006';...
                '2014_04_28_0013';...
                '2014_04_28_0014';...
                '2014_04_28_0015';...
                '2014_04_28_0018';...
                '2014_04_28_0019';...
                '2014_04_28_0020'};      % File names of the raw data. <Nx1> cell array
params.skipSweeps = {}; % In case I need to ignore certain sweeps
params.vHold = -85.*ones(size(params.files));     % The holding potential for vClamp experiments. One value for each expt. <Nx1>
params.stimLoc = [0 0;...
                  -301 197;...
                  -301 197;...
                  -301 197;...
                  -188 154;...
                  -188 154;...
                  -188 154];    % The (x,y) coordinates of the objective at the locations stimulated with the LED. <Nx2> matrix
params.legTxt = {};     % Text that will appear in figures to annotate each data file. <Nx1> cell array
params.tags = {};
params.subtractexp = false;
params.filter = 4000;


%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end


%% EB_041414_C Cell 2

fin

%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PY cell in L2/3 of AL. Some issues with series resistance. Adding to
% population analysis.



%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'EB_041414_C';      % The mouse's name
params.cellNum = 2;    % The neuron number that day
params.DCsteps = '2014_04_28_0023';    % DC steps for Rin and cell identification
params.photo = 'cell_2_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.validCh = 'HS2_';    % 'HS2_' or 'HS1_'
params.files = {'2014_04_28_0026';...
                '2014_04_28_0027';...
                '2014_04_28_0030';...
                '2014_04_28_0031';...
                '2014_04_28_0032';...
                '2014_04_28_0035';...
                '2014_04_28_0036';...
                '2014_04_28_0037'};      % File names of the raw data. <Nx1> cell array
params.skipSweeps = {}; % In case I need to ignore certain sweeps
params.vHold = -85.*ones(size(params.files));     % The holding potential for vClamp experiments. One value for each expt. <Nx1>
params.stimLoc = [0 0;...
                  0 0;...
                  -120 251;...
                  -120 251;...
                  -120 251;...
                  -204 381;...
                  -204 381;...
                  -204 381];    % The (x,y) coordinates of the objective at the locations stimulated with the LED. <Nx2> matrix
params.legTxt = {};     % Text that will appear in figures to annotate each data file. <Nx1> cell array
params.tags = {'Axon stim'};
params.subtractexp = false;
params.filter = 4000;


%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end


%%  EB_041414_D Cell 1

fin

%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Pretty clean looking data from a L2/3 PY cell in PM i think. Series
% resistance is consistent through out.




%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'EB_041414_D';      % The mouse's name
params.cellNum = 1;    % The neuron number that day
params.DCsteps = '2014_04_30_0001';    % DC steps for Rin and cell identification
params.photo = 'cell_1_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.validCh = 'HS2_';    % 'HS2_' or 'HS1_'
params.files = {'2014_04_30_0004';...
                '2014_04_30_0006';...
                '2014_04_30_0007';...
                '2014_04_30_0008';...
                '2014_04_30_0009';...
                '2014_04_30_0010';...
                '2014_04_30_0011'};      % File names of the raw data. <Nx1> cell array
params.skipSweeps = {}; % In case I need to ignore certain sweeps
params.vHold = -85.*ones(size(params.files));     % The holding potential for vClamp experiments. One value for each expt. <Nx1>
params.stimLoc = [0 0;...
                  -302 79;...
                  -302 79;...
                  -302 79;...
                  -370 288;...
                  -370 288;...
                  -370 288];    % The (x,y) coordinates of the objective at the locations stimulated with the LED. <Nx2> matrix
params.legTxt = {};     % Text that will appear in figures to annotate each data file. <Nx1> cell array
params.tags = {};
params.subtractexp = false;
params.filter = 4000;


%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end

%% EB_041414_D cell 1 (light diffusion control)

fin

%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'EB_041414_D';      % The mouse's name
params.cellNum = 1;    % The neuron number that day
params.DCsteps = '2014_04_30_0001';    % DC steps for Rin and cell identification
params.photo = 'cell_1_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.validCh = 'HS2_';    % 'HS2_' or 'HS1_'
params.files = {'2014_04_30_0012';...
                '2014_04_30_0013';...
                '2014_04_30_0014';...
                '2014_04_30_0015'};      % File names of the raw data. <Nx1> cell array
params.skipSweeps = {}; % In case I need to ignore certain sweeps
params.vHold = -85.*ones(size(params.files));     % The holding potential for vClamp experiments. One value for each expt. <Nx1>
params.stimLoc = [-605 136;...
                  -550 186
                  -481 247
                  -419 314];    % The (x,y) coordinates of the objective at the locations stimulated with the LED. <Nx2> matrix
params.legTxt = {};     % Text that will appear in figures to annotate each data file. <Nx1> cell array
params.tags = {'Light spread', 'control'};
params.subtractexp = false;
params.filter = 4000;


%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

%% EB_041414_D cell 2

fin

%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Ok data from a L2/3 PY cell in AL. Ra is pretty variable (goes up then
% down). Adding to pop analysis for AL.


%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'EB_041414_D';      % The mouse's name
params.cellNum = 2;    % The neuron number that day
params.DCsteps = '2014_04_30_0021';    % DC steps for Rin and cell identification
params.photo = 'cell_2_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.validCh = 'HS2_';    % 'HS2_' or 'HS1_'
params.files = {'2014_04_30_0023';...
                '2014_04_30_0025';...
                '2014_04_30_0026';...
                '2014_04_30_0027';...
                '2014_04_30_0028';...
                '2014_04_30_0029';...
                '2014_04_30_0030'};      % File names of the raw data. <Nx1> cell array
params.skipSweeps = {}; % In case I need to ignore certain sweeps
params.vHold = -85.*ones(size(params.files));     % The holding potential for vClamp experiments. One value for each expt. <Nx1>
params.stimLoc = [0 0;...
                  -111 258;...
                  -111 258;...
                  -111 258;...
                  -220 287;...
                  -220 287;...
                  -220 287];    % The (x,y) coordinates of the objective at the locations stimulated with the LED. <Nx2> matrix
params.legTxt = {};     % Text that will appear in figures to annotate each data file. <Nx1> cell array
params.tags = {};
params.subtractexp = false;
params.filter = 4000;


%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end


%% EB_041414_E Cell 1

fin

%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% No picture saved for this cell, so I'm not adding it to the population
% analysis. Not much data anyways...


%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'EB_041414_E';      % The mouse's name
params.cellNum = 1;    % The neuron number that day
params.DCsteps = '2014_04_29_0001';    % DC steps for Rin and cell identification
params.photo = '';      % To assess where the light stimulus was, and the HOA that contains each cell
params.validCh = 'HS2_';    % 'HS2_' or 'HS1_'
params.files = {'2014_04_29_0008';...
                '2014_04_29_0011'};      % File names of the raw data. <Nx1> cell array
params.skipSweeps = {[], [8]}; % In case I need to ignore certain sweeps
params.vHold = -85.*ones(size(params.files));     % The holding potential for vClamp experiments. One value for each expt. <Nx1>
params.stimLoc = [0 0;...
                  -204 170];    % The (x,y) coordinates of the objective at the locations stimulated with the LED. <Nx2> matrix
params.legTxt = {};     % Text that will appear in figures to annotate each data file. <Nx1> cell array
params.tags = {};
params.subtractexp = false;
params.filter = 4000;


%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

%% EB_041414_E Cell 2

fin

%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Stable Ra, only one distal stimulus location. Deep L2/3 PY cell. Not 100%
% sure where this cell is from. It could be from PM or from an area
% directly anterior.

%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'EB_041414_E';      % The mouse's name
params.cellNum = 2;    % The neuron number that day
params.DCsteps = '2014_04_29_0013';    % DC steps for Rin and cell identification
params.photo = 'cell_2_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.validCh = 'HS2_';    % 'HS2_' or 'HS1_'
params.files = {'2014_04_29_0015';...
                '2014_04_29_0017';...
                '2014_04_29_0018';...
                '2014_04_29_0019'};      % File names of the raw data. <Nx1> cell array
params.skipSweeps = {}; % In case I need to ignore certain sweeps
params.vHold = -85.*ones(size(params.files));     % The holding potential for vClamp experiments. One value for each expt. <Nx1>
params.stimLoc = [0 0;...
                  -243 205;...
                  -243 205;...
                  -243 205];    % The (x,y) coordinates of the objective at the locations stimulated with the LED. <Nx2> matrix
params.legTxt = {};     % Text that will appear in figures to annotate each data file. <Nx1> cell array
params.tags = {};
params.subtractexp = false;
params.filter = 4000;


%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end


%% EB_041414_E Cell 3

fin

%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Very good and very stable Ra. FS cell from L2/3 in AL. Adding to pop
% analysis

%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'EB_041414_E';      % The mouse's name
params.cellNum = 3;    % The neuron number that day
params.DCsteps = '2014_04_29_0027';    % DC steps for Rin and cell identification
params.photo = 'cell_3_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.validCh = 'HS2_';    % 'HS2_' or 'HS1_'
params.files = {'2014_04_29_0031';...
                '2014_04_29_0034';...
                '2014_04_29_0035';...
                '2014_04_29_0036';...
                '2014_04_29_0037'};      % File names of the raw data. <Nx1> cell array
params.skipSweeps = {}; % In case I need to ignore certain sweeps
params.vHold = -85.*ones(size(params.files));     % The holding potential for vClamp experiments. One value for each expt. <Nx1>
params.stimLoc = [0 0;...
                  -233 185;...
                  -233 185;...
                  -233 185;...
                  -356 241];    % The (x,y) coordinates of the objective at the locations stimulated with the LED. <Nx2> matrix
params.legTxt = {};     % Text that will appear in figures to annotate each data file. <Nx1> cell array
params.tags = {'axon stim', 'spatial slight spread', 'retro infect'};
params.subtractexp = false;
params.filter = 4000;


%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end


%% EB_042114_A Cell 1

fin

%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'EB_042114_A';      % The mouse's name
params.cellNum = 1;    % The neuron number that day
params.DCsteps = '2014_05_09_0001';    % DC steps for Rin and cell identification
params.photo = 'cell_1_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.validCh = 'HS2_';    % 'HS2_' or 'HS1_'
params.files = {'2014_05_09_0003';...
                '2014_05_09_0005';...
                '2014_05_09_0006';...
                '2014_05_09_0007';...
                '2014_05_09_0008'};      % File names of the raw data. <Nx1> cell array
params.skipSweeps = {}; % In case I need to ignore certain sweeps
params.vHold = -85.*ones(numel(params.files));     % The holding potential for vClamp experiments. One value for each expt. <Nx1>
params.stimLoc = [0 0;...
                  20 230;...
                  20 230;...
                  20 230;...
                  50 384];    % The (x,y) coordinates of the objective at the locations stimulated with the LED. <Nx2> matrix
params.legTxt = {};     % Text that will appear in figures to annotate each data file. <Nx1> cell array
params.tags = {};
params.subtractexp = false;
params.filter = 4000;


%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end


%% EB_042114_A Cell 2

fin

%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'EB_042114_A';      % The mouse's name
params.cellNum = 2;    % The neuron number that day
params.DCsteps = '2014_05_09_0013';    % DC steps for Rin and cell identification
params.photo = 'cell_2_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.validCh = 'HS2_';    % 'HS2_' or 'HS1_'
params.files = {'2014_05_09_0012';...
                '2014_05_09_0015';...
                '2014_05_09_0016';...
                '2014_05_09_0017';...
                '2014_05_09_0018';...
                '2014_05_09_0019';...
                '2014_05_09_0020';...
                '2014_05_09_0021';...
                '2014_05_09_0022'};      % File names of the raw data. <Nx1> cell array
params.skipSweeps = {}; % In case I need to ignore certain sweeps
params.vHold = -85 .* ones(numel(params.files));     % The holding potential for vClamp experiments. One value for each expt. <Nx1>
params.stimLoc = [0 0;...
                  -83 228;...
                  -83 228;...
                  -83 228;...
                  -125 363;...
                  -125 363;...
                  -125 363;...
                  -152 480;...
                  -152 480];    % The (x,y) coordinates of the objective at the locations stimulated with the LED. <Nx2> matrix
params.legTxt = {};     % Text that will appear in figures to annotate each data file. <Nx1> cell array
params.tags = {};
params.subtractexp = false;
params.filter = 4000;


%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end


%% EB_042114_C Cell 1

fin

%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Interneuron in L2/3 of AL. The histolgy didn't work (PBS instead of PFA?)
% but there were good online notes, and some histology photos to suggest
% that this neuron is in AL. Adding to population analysis.

%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'EB_042114_C';      % The mouse's name
params.cellNum = 1;    % The neuron number that day
params.DCsteps = '2014_05_07_0001';    % DC steps for Rin and cell identification
params.photo = 'cell_1_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.validCh = 'HS2_';    % 'HS2_' or 'HS1_'
params.files = {'2014_05_07_0007';...
                '2014_05_07_0013';...
                '2014_05_07_0014';...
                '2014_05_07_0015';...
                '2014_05_07_0016';...
                '2014_05_07_0017';...
                '2014_05_07_0018';...
                '2014_05_07_0019';...
                '2014_05_07_0020';...
                '2014_05_07_0021';...
                '2014_05_07_0022'};      % File names of the raw data. <Nx1> cell array
params.skipSweeps = {}; % In case I need to ignore certain sweeps
params.vHold = -85.*ones(size(params.files));     % The holding potential for vClamp experiments. One value for each expt. <Nx1>
params.stimLoc = [0 0;...
                  131 245;...
                  131 245;...
                  131 245;...
                  131 245;...
                  162 322;...
                  162 322;...
                  162 322;...
                  162 322;...
                  198 420;...
                  198 420];    % The (x,y) coordinates of the objective at the locations stimulated with the LED. <Nx2> matrix
params.legTxt = {};     % Text that will appear in figures to annotate each data file. <Nx1> cell array
params.tags = {};
params.subtractexp = false;
params.filter = 4000;


%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end


%% EB_042114_C Cell 2

fin

%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% These seem like variable data. There's one file where the neuron switches
% b/w facilitation and depression (_0032). Otherwise, good Ra. L2/3 cell in
% AL/LM. Adding to pop anlysis with notes on data quality and ambiguous
% HVA.

%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'EB_042114_C';      % The mouse's name
params.cellNum = 2;    % The neuron number that day
params.DCsteps = '2014_05_07_0025';    % DC steps for Rin and cell identification
params.photo = 'cell_2_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.validCh = 'HS2_';    % 'HS2_' or 'HS1_'
params.files = {'2014_05_07_0027';...
                '2014_05_07_0029';...
                '2014_05_07_0030';...
                '2014_05_07_0031';...
                '2014_05_07_0032';...
                '2014_05_07_0033';...
                '2014_05_07_0034';...
                '2014_05_07_0035';...
                '2014_05_07_0036';...
                '2014_05_07_0037';...
                '2014_05_07_0038'};      % File names of the raw data. <Nx1> cell array
params.skipSweeps = {}; % In case I need to ignore certain sweeps
params.vHold = -85.*ones(size(params.files));     % The holding potential for vClamp experiments. One value for each expt. <Nx1>
params.stimLoc = [0 0;...
                  -63 257;...
                  -63 257;...
                  -63 257;...
                  -222 355;...
                  -222 355;...
                  -222 355;...
                  -63 257;...
                  -63 257;...
                  0 0;...
                  0 0];    % The (x,y) coordinates of the objective at the locations stimulated with the LED. <Nx2> matrix
params.legTxt = {};     % Text that will appear in figures to annotate each data file. <Nx1> cell array
params.tags = {};
params.subtractexp = false;
params.filter = 4000;


%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end


%% EB_150326_A cell 2

fin

%
% Notes
%
%%%%%%%%%%%%%%%%%
% 
%

%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'EB_150326_A';      % The mouse's name
params.cellNum = 2;    % The neuron number that day
params.photo = 'EB_150326_A_cell2';      % To assess where the light stimulus was, and the HOA that contains each cell
params.files = {'2015_04_17_', [4]};  % <file name prefix, suffix>
params.groups = {'recovery', [4]};
params.excludeHS1 = {};
params.excludeHS2 = {}; 


% stuff for E/I and AMPA/NMDA ratios
% key for isolatedCurrents = {<current><group><Vhold><Erev>}
% Erev is to calculate driving force for conversion from pA to pS
params.isolatedCurrents = {'excit', 'control', -85, 17};
params.tags = {};
params.filter = 2e3;

HS1loc = [];
HS2loc = [];
Pialoc1 = [];
Pialoc2 = [];
params.celldepth = [norm(HS1loc-Pialoc1), norm([HS2loc-Pialoc2])];


%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_avgOuterleave};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end

%% EB_150427_C cell 1

fin

%
% NOTES (cell type, location, experiment type, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% GIN cell recorded in response to ChR2 activation at different LED
% positions and pulse train frequencies. 


%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%
params.mouse = 'EB_150427_C';      % The mouse's name
params.cellNum = 1;    % The neuron number that day
params.DCsteps = '2015_05_14_0000';    % DC steps for Rin and cell identification
params.photo = 'EB_150427_C_cell1_tdTomato';      % To assess where the light stimulus was, and the HOA that contains each cell
params.validCh = 'HS2_';    % 'HS2_' or 'HS1_'
params.files = {'2015_05_14_0001';...
                '2015_05_14_0002';...
                '2015_05_14_0003';...
                '2015_05_14_0004';...
                '2015_05_14_0005';...
                '2015_05_14_0006';...
                '2015_05_14_0007'};      % File names of the raw data. <Nx1> cell array
params.skipSweeps = {}; % In case I need to ignore certain sweeps
params.vHold = -87 .* ones(numel(params.files));     % The holding potential for vClamp experiments. One value for each expt. <Nx1>
stimSite_raw = [56 101;...
                56 101;...
                -168 -130;...
                -22 -13;...
                -22 -13;...
                -168 -138;...
                -168 -138];
params.stimLoc = bsxfun(@minus, stimSite_raw, [-22 -13]);
params.legTxt = {'loc2 20hz o 10V';...
                 'loc2 40hz o 10V';...
                 'loc3 40hz o 10V';...
                 'loc1 40hz o 10V';...
                 'loc1 40hz o 5V';...
                 'loc3 40hz o 5V';...
                 'loc3 40hz c 10V'};     % Text that will appear in figures to annotate each data file. <Nx1> cell array
params.tags = {};
params.subtractexp = false;
params.filter = 4000;


%
% ANALYZE OR ADD TO PARAMSDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GL_SUPPRESS_ANALYSIS', 'var') || ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if exist('GL_ADD_TO_MDB', 'var') && GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end




%% BOOKMARK FOR EB MICE
% Added to population analyses up to here




