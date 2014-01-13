function [C, trialParamsHeader, exptParamsHeader, additionalRasterFields, badLogic, otherInstructions] = GaborEdgeCodes()

% This is a collection of constant definitions for use with analysis
% programs.  These codes are specific to the WhiteNoise paradigm

% The offset that is added to the value ecodes that follow the
% identifier ecodes (7000-9000 range).  Note, this value must be
% identical to the definition of VALODFFSET in rig1.h on the REX machine.
C.VALOFFSET = 4000;

% Ecodes that identify the header parameters
% (Parameters that do not change on a trial to trial basis
% and so are dropped once, at the begining of the run)
C.FIXXCD = 8000;
C.FIXYCD = 8001;
C.FPSIZECD = 8002;
C.FPRCD = 8003;
C.FPGCD = 8004;
C.FPBCD = 8005;
C.EYEWINXCD = 8006;
C.EYEWINYCD = 8007;
C.RFXCD = 8008;
C.RFYCD = 8009;
C.BKGNDRGBCD = 8010; %double,1
C.GAMMATABLECD = 8011; %double,1
C.MONSPDCD = 8012; %double,1
C.FRAMERATECD = 8013; %double, 1
C.FUNDAMENTALSCD = 8014; %double, 1
C.TRLSPERCONDCD = 8015;
C.GAUSSLIMCD = 8016;
C.NFRAMESPERDISPLAYCD = 8017;
C.GABORRAWRGBCD = 8018; %double,1
C.PIXPERDEGCD = 8019; % double, 1
C.EPCOMPENSATION = 8020;
C.NEPSAMPLESTOAVERAGE = 8021;

% Ecodes that identify the trial parameters
% (Parameters that change on a trial to trial basis
% and so are dropped at the end of each trial)
C.GABORTHETACD = 7000;
C.GABORLAMBDACD = 7001;
C.GABORPHICD = 7002;
C.GABORSIGMACD = 7003;
C.GABORGAMMACD = 7004;
C.GABORXOFFSETCD = 7005;
C.GABORYOFFSETCD = 7006;
C.GABORCONTRASTCD = 7007;
C.GABORCOLORTYPECD = 7008;
C.GABORRCD = 7009;
C.GABORGCD = 7010;
C.GABORBCD = 7011;
C.GABORLCD = 7012;
C.GABORMCD = 7013;
C.GABORSCD = 7014;
C.EDGETHETACD = 7015;
C.EDGEDISPLACEMENTCD = 7016;
C.EDGECONTRASTCD = 7017;
C.EDGECOLORTYPECD = 7018;
C.EDGERCD = 7019;
C.EDGEGCD = 7020;
C.EDGEBCD = 7021;
C.EDGELCD = 7022;
C.EDGEMCD = 7023;
C.EDGESCD = 7024;

% Ecodes that are used to mark the times of various trial events
C.FPONCD = 1004;
C.FPACQCD = 1032;
C.STIMONCD = 1030;
C.STIMOFFCD = 1031;
C.FPOFFCD = 1006;
C.REWCD = 1015;
C.EOTCD = 1016;
C.ABORTCD = 1013;

%******************%
%   for nex2stro   %
%******************%

% includes the trial column identifiers, the ecode id's, the data type, and
% bookending info.
trialParamsHeader = {'fp_on', 'fp_acq', 'stim_on', 'stim_off', 'fp_off', 'g_theta', 'g_lambda', 'g_phi', 'g_sigma', 'g_gamma', 'g_xoff', 'g_yoff', 'g_cont', 'g_colortype', 'g_r',   'g_g',   'g_b',   'g_l',   'g_m',   'g_s',   'e_theta', 'e_disp', 'e_cont',  'e_colortype', 'e_r',   'e_g',   'e_b',    'e_l',  'e_m',   'e_s';...
                     'time',  'time',   'time',    'time',     'time',   'float',   'float',    'float', 'float',   'float',   'float',  'float',  'float',  'int',         'double','double','double','double','double','double','float',   'float',  'float',   'int',         'double','double','double','double','double','double';...
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;  ...
                       C.FPONCD, C.FPACQCD, C.STIMONCD, C.STIMOFFCD, C.FPOFFCD, C.GABORTHETACD, C.GABORLAMBDACD, C.GABORPHICD, C.GABORSIGMACD, C.GABORGAMMACD, C.GABORXOFFSETCD, C.GABORYOFFSETCD, C.GABORCONTRASTCD, C.GABORCOLORTYPECD, C.GABORRCD, C.GABORGCD, C.GABORBCD, C.GABORLCD, C.GABORMCD, C.GABORSCD, C.EDGETHETACD, C.EDGEDISPLACEMENTCD, C.EDGECONTRASTCD, C.EDGECOLORTYPECD, C.EDGERCD, C.EDGEGCD, C.EDGEBCD, C.EDGELCD, C.EDGEMCD, C.EDGESCD};

%for the stro.sum.expParams array (includes things that get dropped only
%once):
exptParamsHeader = {'fp_x', 'fp_y', 'fp_size', 'fp_r', 'fp_g', 'fp_b', 'eyewin_x', 'eyewin_y', 'rf_x', 'rf_y', 'bkgndrgb', 'gamma_table', 'mon_spd', 'framerate', 'fundamentals', 'trialspercond','gausslim','nframes','gaborrawrgb', 'pixperdeg', 'epcomp', 'nepsamp';...
                    'int',  'int',  'int',     'int',  'int',  'int',  'int',      'int',      'int',  'int',  'double',   'double',      'double',  'double',    'double',       'int',          'int',     'int',    'double',      'double',    'int',    'int';...
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0;...
                    C.FIXXCD, C.FIXYCD, C.FPSIZECD, C.FPRCD, C.FPGCD, C.FPBCD, C.EYEWINXCD, C.EYEWINYCD, C.RFXCD, C.RFYCD, C.BKGNDRGBCD, C.GAMMATABLECD, C.MONSPDCD, C.FRAMERATECD, C.FUNDAMENTALSCD, C.TRLSPERCONDCD, C.GAUSSLIMCD, C.NFRAMESPERDISPLAYCD, C.GABORRAWRGBCD, C.PIXPERDEGCD, C.EPCOMPENSATION, C.NEPSAMPLESTOAVERAGE};

additionalRasterFields = {};


%define how good and bad trials are to be distinguished
badLogic = 'find(trialCodes(:,2) == C.ABORTCD, 1)';

%define other routines as an array of function handles. the first row of
%the array should be a name of the handle, and the second row should be the
%actual function handle. The functions should be included as non-nested
%subfunctions in this file. DTStimCodes.m serves as a template
otherInstructions = {}; %none for now

end