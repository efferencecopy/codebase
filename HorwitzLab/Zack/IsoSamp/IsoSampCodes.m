function [C, trialParamsHeader, exptParamsHeader, additionalRasterFields, badLogic, otherInstructions] = IsoSampCodes()

% This is a collection of constant definitions for use with analysis
% programs.  These codes are specific to the WhiteNoise paradigm

% The offset that is added to the value ecodes that follow the
% identifier ecodes (7000-9000 range).  Note, this value must be
% identical to the definition of VALODFFSET in rig1.h on the REX machine.
C.VALOFFSET = 4000;

% Ecodes that identify the header parameters
% (Parameters that do not change on a trial to trial basis
% and so are dropped once, at the begining of the run)
C.FIXXCD =		   8000;
C.FIXYCD =		   8001;
C.FPSIZECD =	   8002;
C.FPRCD =		   8003;
C.FPGCD =		   8004;
C.FPBCD =		   8005;
C.EYEWINXCD =	   8006;
C.EYEWINYCD =	   8007;
C.RFXCD = 		   8008;
C.RFYCD =		   8009;
C.MONKEYCD =       8010;
C.NUMSTIMSCD =     8011;
C.BKGNDRGBCD =     8012; % double, 1
C.GAMMATABLECD =   8013; % double, 1
C.MONSPDCD =       8014; % double, 1
C.FRAMERATECD =	   8015; % double, 0
C.FUNDAMENTALSCD = 8016; % double, 1
C.PIXPERDEGCD =	   8017; % double, 0
C.MODELPARAMSCD =  8018; % double, 1
C.PLOTLIMITSCD =   8019; % double, 1
C.HDRCOMPLETECD =  8998;

% Ecodes that identify the trial parameters
% (Parameters that change on a trial to trial basis
% and so are dropped at the end of each trial)
C.GABFLASHTIMECD = 7000;
C.GABPHICD =       7001;
C.GABTHETACD =     7002;
C.GABSFCD =        7003;
C.GABGAMMACD =     7004;
C.GABTFCD =        7005;
C.GABNSIGMASCD =   7006;
C.GABSIGMACD =     7007;
C.GABLCD =         7008;
C.GABMCD =         7009;
C.GABSCD =         7010;

% Ecodes that are used to mark the times of various trial events
C.FPONCD =      1004;
C.FPOFFCD =     1006;
C.ABORTCD =     1013;
C.REWCD =       1015;
C.EOTCD =       1016;
C.CORRECTCD =   1018;
C.STIMONCD =    1030;
C.STIMOFFCD =	1031;
C.FPACQCD =     1032;

%******************%
%   for nex2stro   %
%******************%

% includes the trial column identifiers, the ecode id's, the data type, and
% bookending info.
trialParamsHeader = {
    'fpon_t',       'time',   0, C.FPONCD
    'fpacq_t',      'time',   0, C.FPACQCD
    'fpoff_t',      'time',   0, C.FPOFFCD
    'stimon_t',     'time',   0, C.STIMONCD
    'stimoff_t',    'time',   0, C.STIMOFFCD
    'rew_t',        'time',   0, C.REWCD
	'flash_time',   'int',    0, C.GABFLASHTIMECD
    'phi',          'float',  0, C.GABPHICD
    'theta',        'double', 0, C.GABTHETACD
    'sf',           'double', 0, C.GABSFCD
    'gamma',        'float',  0, C.GABGAMMACD
    'tf',           'float',  0, C.GABTFCD
    'sigmas_n',     'float',  0, C.GABNSIGMASCD
    'sigma',        'float',  0, C.GABSIGMACD
    'stim_l',       'double', 0, C.GABLCD
    'stim_m',       'double', 0, C.GABMCD
    'stim_s',       'double', 0, C.GABSCD
}';

%for the stro.sum.exptParams array (includes things that get dropped only
%once):
exptParamsHeader = {
    'fp_x',            'int', 0, C.FIXXCD
    'fp_y',            'int', 0, C.FIXYCD
    'fp_size',         'int', 0, C.FPSIZECD
    'fp_r',            'int', 0, C.FPRCD
    'fp_g',            'int', 0, C.FPGCD
    'fp_b',            'int', 0, C.FPBCD
    'eyewin_x',        'int', 0, C.EYEWINXCD
    'eyewin_y',        'int', 0, C.EYEWINYCD
    'rf_x',            'int', 0, C.RFXCD
    'rf_y',            'int', 0, C.RFYCD
    'bkgndrgb',     'double', 1, C.BKGNDRGBCD
    'gamma_table',  'double', 1, C.GAMMATABLECD
    'mon_spd',      'double', 1, C.MONSPDCD
    'framerate',    'double', 0, C.FRAMERATECD
    'fundamentals', 'double', 1, C.FUNDAMENTALSCD
    'pixperdeg',    'double', 0, C.PIXPERDEGCD
    'modelparams',  'double', 1, C.MODELPARAMSCD
    'plotlimits',   'double', 1, C.PLOTLIMITSCD
    'monkey_id',       'int', 0, C.MONKEYCD
    'nstims_block',    'int', 0, C.NUMSTIMSCD
}';

additionalRasterFields = {};

%define how good and bad trials are to be distinguished
badLogic = 'find(trialCodes(:,2) == C.ABORTCD, 1)';

%define other routines as an array of function handles. the first row of
%the array should be a name of the handle, and the second row should be the
%actual function handle. The functions should be included as non-nested
%subfunctions in this file. DTStimCodes.m serves as a template
otherInstructions = {}; %none for now

end
