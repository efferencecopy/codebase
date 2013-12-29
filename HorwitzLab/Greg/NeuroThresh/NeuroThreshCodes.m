function [C, trialParamsHeader, exptParamsHeader, additionalRasterFields, badLogic, otherInstructions] = NeuroThreshCodes()

% This is a collection of constant definitions for use with analysis
% programs.  These codes are specific to the Grating paradigm

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
C.NFRAMESPLATCD = 8010;
C.NFRAMESRAMPCD = 8011;
C.NSIGMASCD = 8012;
C.BKGNDRGBCD = 8013; % double,1
C.GAMMATABLECD = 8014; %double,1
C.MONSPDCD = 8015; %double,1
C.FRAMERATECD = 8016; %double, 0
C.FUNDAMENTALSCD = 8017; %double, 1
C.PIXPERDEGCD = 8018; %double, 0
C.INITLMS1CD = 8019; %float, 1
C.INITLMS2CD = 8020; %float, 1
C.INITLMS3CD = 8021; %float, 1
C.INITSTEPSIZECD = 8022; %float, 0
C.STEPSIZESCALECD = 8023; %float, 0
C.NREVERSALSCD = 8024;
C.NNEWDIRSCD = 8025;
C.LINPREDTOLCD = 8026;% float, 0
C.THRESHCD = 8027;% float, 0
C.PREFLMSCD = 8028;% float, 1
C.LATENCYCD = 8029;% float, 0
C.OOGSCALECD = 8030;% float, 0
C.MONOPOLARCD = 8031;
C.NPGRIDCD = 8032;
C.WHICHKERNELCD = 8033;
C.RANDSTIMSELCD = 8034;
C.LHSTIMSELCD = 8035;
C.ADAPTSTIMSELCD = 8036;

% Ecodes that identify the trial parameters
% (Parameters that change on a trial to trial basis
% and so are dropped at the end of each trial)
C.GABORTHETACD = 7000;   % float
C.GABORSFCD = 7001;      % float
C.GABORPHICD = 7002;     % float
C.GABORSIGMACD = 7003;
C.GABORGAMMACD = 7004;
C.GABORLCD = 7005;
C.GABORMCD = 7006;
C.GABORSCD = 7007;
C.COLORIDXCD = 7008;
C.REVERSALCD = 7009;
C.STEPSIZECD = 7010;
C.OVERTHRESHCD = 7011;
C.LEVELIDXCD = 7012;

% Ecodes that are used to mark the times of various trial events
C.FPONCD = 1004;
C.FPACQCD = 1032;
C.STIMONCD = 1030;
C.STIMOFFCD = 1031;
C.FPOFFCD = 1006;
C.CORRECTCD = 1018;
C.REWCD = 1015;
C.EOTCD = 1016;
C.ABORTCD = 1013;

%******************%
%   for nex2stro   %
%******************%

% includes the trial column identifiers, the ecode id's, the data type, and
% bookending info.
trialParamsHeader = {'fp_on', 'fp_acq', 'stim_on', 'stim_off', 'eot', 'lcont', 'mcont', 'scont', 'coloridx', 'reversal', 'stepsize', 'overthresh', 'levelidx', 'sf',    'orient';...
                      'time',  'time',   'time',    'time',     'time','float', 'float', 'float', 'int',      'int',      'float',    'int',        'int',      'float', 'float';...
                      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;...
                      C.FPONCD, C.FPACQCD, C.STIMONCD, C.STIMOFFCD, C.EOTCD, C.GABORLCD, C.GABORMCD, C.GABORSCD, C.COLORIDXCD, C.REVERSALCD, C.STEPSIZECD, C.OVERTHRESHCD, C.LEVELIDXCD, C.GABORSFCD, C.GABORTHETACD};
 
%for the stro.sum.expParams array (includes things that get dropped only
%once):
exptParamsHeader = {'fp_x', 'fp_y', 'fp_size', 'fp_r', 'fp_g', 'fp_b', 'eyewin_x', 'eyewin_y', 'rf_x', 'rf_y', 'nframesplateau', 'nframesramp', 'nreversals', 'stepsize', 'scale', 'threshold', 'latency', 'nnewdirsperlev','framerate', 'gamma_table', 'mon_spd', 'fundamentals', 'initLMS1', 'initLMS2', 'initLMS3', 'bkgndrgb', 'preflms', 'linpredtol', 'oogscale', 'monopolar';...
                    'int',  'int',  'int',     'int',  'int',  'int',  'int',      'int',      'int',  'int',  'int',            'int',         'int',        'float',    'float', 'float',     'float',   'int'            'double',    'double',      'double',  'double',       'float',    'float',    'float',    'double',   'float',   'float',        'float',    'int'; ...
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0; ...
                    C.FIXXCD, C.FIXYCD, C.FPSIZECD, C.FPRCD, C.FPGCD, C.FPBCD, C.EYEWINXCD, C.EYEWINYCD, C.RFXCD, C.RFYCD, C.NFRAMESPLATCD, C.NFRAMESRAMPCD, C.NREVERSALSCD, C.INITSTEPSIZECD, C.STEPSIZESCALECD, C.THRESHCD, C.LATENCYCD, C.NNEWDIRSCD, C.FRAMERATECD, C.GAMMATABLECD, C.MONSPDCD, C.FUNDAMENTALSCD, C.INITLMS1CD, C.INITLMS2CD, C.INITLMS3CD, C.BKGNDRGBCD, C.PREFLMSCD, C.LINPREDTOLCD, C.OOGSCALECD, C.MONOPOLARCD};

additionalRasterFields = {};
 
%define how good and bad trials are to be distinguished
badLogic = 'find(trialCodes(:,2) == C.ABORTCD, 1)';

%define other routines as an array of function handles. the first row of
%the array should be a name of the handle, and the second row should be the
%actual function handle. The functions should be included as non-nested
%subfunctions in this file. DTStimCodes.m serves as a template
otherInstructions = {}; %none for now

end