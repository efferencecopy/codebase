%% DEFINE THE TERMS OF THE ANALYSIS

global GL_ADD_TO_MDB GL_SUPPRESS_ANALYSIS

GL_ADD_TO_MDB = false;
GL_SUPPRESS_ANALYSIS = false;

%% TEMPLATE VERSION (mouse name and cell num)

fin

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
params.validCh = 'HS2_';    % 'HS2_' or 'HS1_'
params.files = {};      % File names of the raw data. <Nx1> cell array
params.skipSweeps = {}; % In case I need to ignore certain sweeps
params.vHold = nan;     % The holding potential for vClamp experiments. One value for each expt. <Nx1>
params.stimLoc = [];    % The (x,y) coordinates of the objective at the locations stimulated with the LED. <Nx2> matrix
params.legTxt = {};     % Text that will appear in figures to annotate each data file. <Nx1> cell array
params.tags = {};
params.subtractexp = GL_ADD_TO_MDB;
params.filter = 4000;

% 
% %
% % ANALYZE OR ADD TO MDB
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if ~GL_SUPPRESS_ANALYSIS
%     params.fxns = {@anlyMod_pulseTrains_stimLoc};
%     params = invitroAnalysisOverview(params);
% end
% 
% if GL_ADD_TO_MDB
%     addPopAnlyParamsToMDB(params);
% end


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
params.subtractexp = GL_ADD_TO_MDB;
params.filter = 4000;

%
% ANALYZE OR ADD TO MDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if GL_ADD_TO_MDB
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
params.subtractexp = GL_ADD_TO_MDB;
params.filter = 4000;

%
% ANALYZE OR ADD TO MDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if GL_ADD_TO_MDB
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
params.subtractexp = GL_ADD_TO_MDB;
params.filter = 4000;

%
% ANALYZE OR ADD TO MDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if GL_ADD_TO_MDB
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
params.subtractexp = GL_ADD_TO_MDB;
params.filter = 4000;

%
% ANALYZE OR ADD TO MDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if GL_ADD_TO_MDB
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
params.subtractexp = GL_ADD_TO_MDB;
params.filter = 4000;

%
% ANALYZE OR ADD TO MDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if GL_ADD_TO_MDB
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
params.subtractexp = GL_ADD_TO_MDB;
params.filter = 4000;

%
% ANALYZE OR ADD TO MDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if GL_ADD_TO_MDB
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
params.subtractexp = GL_ADD_TO_MDB;
params.filter = 4000;

%
% ANALYZE OR ADD TO MDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if GL_ADD_TO_MDB
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
params.subtractexp = GL_ADD_TO_MDB;
params.filter = 4000;


%
% ANALYZE OR ADD TO MDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if GL_ADD_TO_MDB
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
params.subtractexp = GL_ADD_TO_MDB;
params.filter = 4000;


%
% ANALYZE OR ADD TO MDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if GL_ADD_TO_MDB
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
params.subtractexp = GL_ADD_TO_MDB;
params.filter = 4000;


%
% ANALYZE OR ADD TO MDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~GL_SUPPRESS_ANALYSIS
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
params.subtractexp = GL_ADD_TO_MDB;
params.filter = 4000;


%
% ANALYZE OR ADD TO MDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~GL_SUPPRESS_ANALYSIS
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
% increases... Facilitation at soma, depression otherwise, accept at high
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
params.subtractexp = GL_ADD_TO_MDB;
params.filter = 4000;


%
% ANALYZE OR ADD TO MDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if GL_ADD_TO_MDB
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
params.subtractexp = GL_ADD_TO_MDB;
params.filter = 4000;


%
% ANALYZE OR ADD TO MDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end


%% %% CH_050614_B cell 1

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
params.subtractexp = GL_ADD_TO_MDB;
params.filter = 4000;


%
% ANALYZE OR ADD TO MDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if GL_ADD_TO_MDB
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
params.subtractexp = GL_ADD_TO_MDB;
params.filter = [1000];


%
% ANALYZE OR ADD TO MDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if GL_ADD_TO_MDB
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
params.subtractexp = GL_ADD_TO_MDB;
params.filter = 500;


%
% ANALYZE OR ADD TO MDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end


%% BOOKMARK FOR CH MICE
% data added to population analyses up to here


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
params.subtractexp = GL_ADD_TO_MDB;
params.filter = 4000;


%
% ANALYZE OR ADD TO MDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~GL_SUPPRESS_ANALYSIS
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
params.subtractexp = GL_ADD_TO_MDB;
params.filter = 4000;

%
% ANALYZE OR ADD TO MDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~GL_SUPPRESS_ANALYSIS
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
params.subtractexp = GL_ADD_TO_MDB;
params.filter = 4000;

%
% ANALYZE OR ADD TO MDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if GL_ADD_TO_MDB
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
params.subtractexp = GL_ADD_TO_MDB;
params.filter = 4000;

%
% ANALYZE OR ADD TO MDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if GL_ADD_TO_MDB
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
params.subtractexp = GL_ADD_TO_MDB;
params.filter = 4000;

%
% ANALYZE OR ADD TO MDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if GL_ADD_TO_MDB
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
params.subtractexp = GL_ADD_TO_MDB;
params.filter = 4000;

%
% ANALYZE OR ADD TO MDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if GL_ADD_TO_MDB
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
params.subtractexp = GL_ADD_TO_MDB;
params.filter = 4000;

%
% ANALYZE OR ADD TO MDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if GL_ADD_TO_MDB
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
params.subtractexp = GL_ADD_TO_MDB;
params.filter = 4000;

%
% ANALYZE OR ADD TO MDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if GL_ADD_TO_MDB
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
params.subtractexp = GL_ADD_TO_MDB;
params.filter = 4000;

%
% ANALYZE OR ADD TO MDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if GL_ADD_TO_MDB
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
params.subtractexp = GL_ADD_TO_MDB;
params.filter = 4000;

%
% ANALYZE OR ADD TO MDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if GL_ADD_TO_MDB
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
params.subtractexp = GL_ADD_TO_MDB;
params.filter = 4000;

%
% ANALYZE OR ADD TO MDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if GL_ADD_TO_MDB
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
params.subtractexp = GL_ADD_TO_MDB;
params.filter = 4000;

%
% ANALYZE OR ADD TO MDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if GL_ADD_TO_MDB
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
params.subtractexp = GL_ADD_TO_MDB;
params.filter = 4000;

%
% ANALYZE OR ADD TO MDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if GL_ADD_TO_MDB
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
params.subtractexp = GL_ADD_TO_MDB;
params.filter = 4000;

%
% ANALYZE OR ADD TO MDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if GL_ADD_TO_MDB
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
params.subtractexp = GL_ADD_TO_MDB;
params.filter = 4000;

%
% ANALYZE OR ADD TO MDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if GL_ADD_TO_MDB
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
params.subtractexp = GL_ADD_TO_MDB;
params.filter = 4000;

%
% ANALYZE OR ADD TO MDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if GL_ADD_TO_MDB
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
params.subtractexp = GL_ADD_TO_MDB;
params.filter = 4000;


%
% ANALYZE OR ADD TO MDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if GL_ADD_TO_MDB
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
params.subtractexp = GL_ADD_TO_MDB;
params.filter = 4000;


%
% ANALYZE OR ADD TO MDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if GL_ADD_TO_MDB
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
params.subtractexp = GL_ADD_TO_MDB;
params.filter = 4000;


%
% ANALYZE OR ADD TO MDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if GL_ADD_TO_MDB
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
params.subtractexp = GL_ADD_TO_MDB;
params.filter = 4000;


%
% ANALYZE OR ADD TO MDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if GL_ADD_TO_MDB
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
params.subtractexp = GL_ADD_TO_MDB;
params.filter = 4000;


%
% ANALYZE OR ADD TO MDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if GL_ADD_TO_MDB
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
params.subtractexp = GL_ADD_TO_MDB;
params.filter = 4000;


%
% ANALYZE OR ADD TO MDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if GL_ADD_TO_MDB
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
params.subtractexp = GL_ADD_TO_MDB;
params.filter = 4000;


%
% ANALYZE OR ADD TO MDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if GL_ADD_TO_MDB
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
params.subtractexp = GL_ADD_TO_MDB;
params.filter = 4000;


%
% ANALYZE OR ADD TO MDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if GL_ADD_TO_MDB
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
params.subtractexp = GL_ADD_TO_MDB;
params.filter = 4000;


%
% ANALYZE OR ADD TO MDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if GL_ADD_TO_MDB
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
params.subtractexp = GL_ADD_TO_MDB;
params.filter = 4000;


%
% ANALYZE OR ADD TO MDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if GL_ADD_TO_MDB
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
params.subtractexp = GL_ADD_TO_MDB;
params.filter = 4000;


%
% ANALYZE OR ADD TO MDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if GL_ADD_TO_MDB
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
params.subtractexp = GL_ADD_TO_MDB;
params.filter = 4000;


%
% ANALYZE OR ADD TO MDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if GL_ADD_TO_MDB
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
params.subtractexp = GL_ADD_TO_MDB;
params.filter = 4000;


%
% ANALYZE OR ADD TO MDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~GL_SUPPRESS_ANALYSIS
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
params.subtractexp = GL_ADD_TO_MDB;
params.filter = 4000;


%
% ANALYZE OR ADD TO MDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if GL_ADD_TO_MDB
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
params.subtractexp = GL_ADD_TO_MDB;
params.filter = 4000;


%
% ANALYZE OR ADD TO MDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~GL_SUPPRESS_ANALYSIS
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
params.subtractexp = GL_ADD_TO_MDB;
params.filter = 4000;


%
% ANALYZE OR ADD TO MDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if GL_ADD_TO_MDB
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
params.subtractexp = GL_ADD_TO_MDB;
params.filter = 4000;


%
% ANALYZE OR ADD TO MDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if GL_ADD_TO_MDB
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
params.subtractexp = GL_ADD_TO_MDB;
params.filter = 4000;


%
% ANALYZE OR ADD TO MDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if GL_ADD_TO_MDB
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
params.subtractexp = GL_ADD_TO_MDB;
params.filter = 4000;


%
% ANALYZE OR ADD TO MDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if GL_ADD_TO_MDB
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
params.subtractexp = GL_ADD_TO_MDB;
params.filter = 4000;


%
% ANALYZE OR ADD TO MDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if GL_ADD_TO_MDB
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
params.subtractexp = GL_ADD_TO_MDB;
params.filter = 4000;


%
% ANALYZE OR ADD TO MDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~GL_SUPPRESS_ANALYSIS
    params.fxns = {@anlyMod_pulseTrains_stimLoc};
    params = invitroAnalysisOverview(params);
end

if GL_ADD_TO_MDB
    addPopAnlyParamsToMDB(params);
end

%% BOOKMARK FOR EB MICE
% Added to population analyses up to here




