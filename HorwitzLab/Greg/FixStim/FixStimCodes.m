function [C, trialParamsHeader, exptParamsHeader, additionalRasterFields, badLogic, otherInstructions] = FixStimCodes()

% This is a collection of constant definitions for use with analysis
% programs.  These codes are specific to the WhiteNoise paradigm

% The offset that is added to the value ecodes that follow the
% identifier ecodes (7000-9000 range).  Note, this value must be
% identical to the definition of VALODFFSET in rig1.h on the REX machine.
C.VALOFFSET = 4000;

% Ecodes that identify the header parameters
% (Parameters that do not change on a trial to trial basis
% and so are dropped once, at the begining of the run)
C.FIXXCD            = 8000;
C.FIXYCD            = 8001;
C.FPSIZECD          = 8002;
C.FPRCD             = 8003;
C.FPGCD             = 8004;
C.FPBCD             = 8005;
C.EYEWINXCD         = 8006;
C.EYEWINYCD         = 8007;
C.TARGSIZECD        = 8008;
C.TARGRCD           = 8009;
C.TARGGCD           = 8010;
C.TARGBCD           = 8011;
C.TARGWINXCD        = 8012;
C.TARGWINYCD        = 8013;
C.RFXCD             = 8014;
C.RFYCD             = 8015;
C.NTRLSPERCONDCD	= 8016;
C.BKGNDRGBCD        = 8017;
C.GAMMATABLECD      = 8018;
C.FUNDAMENTALSCD    = 8019;
C.MONSPDCD          = 8020;

% Ecodes that identify the trial parameters
% (Parameters that change on a trial to trial basis
% and so are dropped at the end of each trial)
C.TARGSHOWNCD =     7000;
C.STIMTYPECD =      7001;
C.ELECSTIMFREQCD =	7002;
C.LASERPOWERCD =    7003;
C.CONTRASTCD =      7004;
C.STIMHEIGHT =      7005;
C.STIMWIDTH =       7006;

% Ecodes that are used to mark the times of various trial events
C.FPONCD =      1004;
C.FPOFFCD =     1006;
C.TARGONCD =	1008;
C.TARGOFFCD =	1009;
C.ABORTCD =     1013;
C.REWCD =       1015;
C.EOTCD =       1016;
C.STIMONCD =    1030;
C.STIMOFFCD =	1031;
C.FPACQCD =     1032;

%******************%
%   for nex2stro   %
%******************%

% includes the trial column identifiers, the ecode id's, the data type, and
% bookending info.
trialParamsHeader = {
    'fpon_t',       'time', 0, C.FPONCD
    'fpacq_t',      'time', 0, C.FPACQCD
    'fpoff_t',      'time', 0, C.FPOFFCD
    'stimon_t',     'time', 0, C.STIMONCD
    'stimoff_t',    'time', 0, C.STIMOFFCD
    'targon_t',     'time', 0, C.TARGONCD
    'targoff_t',    'time',	0, C.TARGOFFCD
    'targ_shown',    'int',	0, C.TARGSHOWNCD
    'stim_type',     'int',	0, C.STIMTYPECD
    'elec_stimfreq', 'int',	0, C.ELECSTIMFREQCD
    'laser_power',   'int', 0, C.LASERPOWERCD
    'contrast',      'float', 0, C.CONTRASTCD
    }';

%for the stro.sum.expParams array (includes things that get dropped only
%once):
exptParamsHeader = {
    'fp_x',        'int', 0, C.FIXXCD
    'fp_y',        'int', 0, C.FIXYCD
    'fp_size',     'int', 0, C.FPSIZECD
    'fp_r',        'int', 0, C.FPRCD
    'fp_g',        'int', 0, C.FPGCD
    'fp_b',        'int', 0, C.FPBCD
    'eyewin_x',    'int', 0, C.EYEWINXCD
    'eyewin_y',    'int', 0, C.EYEWINYCD
    'targ_size',   'int', 0, C.TARGSIZECD
    'targ_r',      'int', 0, C.TARGRCD
    'targ_g',      'int', 0, C.TARGGCD
    'targ_b',      'int', 0, C.TARGBCD
    'targwin_x',   'int', 0, C.TARGWINXCD
    'targwin_y',   'int', 0, C.TARGWINYCD
    'rf_x',        'int', 0, C.RFXCD
    'rf_y',        'int', 0, C.RFYCD
    'fundamentals','double',1, C.FUNDAMENTALSCD
    'trlspercond', 'int', 0, C.NTRLSPERCONDCD
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
