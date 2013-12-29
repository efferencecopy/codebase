function [C, trialParamsHeader, exptParamsHeader, additionalRasterFields, badLogic, otherInstructions] = GratingCodes()

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
C.BKGNDRGBCD = 8010; % double,1
C.GAMMATABLECD = 8011; %double,1
C.MONSPDCD = 8012; %double,1
C.FRAMERATECD = 8013; %double, 1
C.FUNDAMENTALSCD = 8014; %double, 1
C.INITPROTOCOLCD = 8015;
C.NTRLSPERCONDCD = 8016;
C.INITDIAMCD = 8017;
C.INITSFCD = 8018;
C.INITTFCD = 8019;
C.INITORIENTCD = 8020;
C.INITPHASECD = 8021;
C.INITLCONTCD = 8022;
C.INITMCONTCD = 8023;
C.INITSCONTCD = 8024;
C.INITNCYCLESCD = 8025;
C.NSIGMASCD = 8026;
C.NCYCLESPLATCD = 8027;
C.NCYCLESRAMPCD = 8028;
C.MAXCONTSCALEFACTCD = 8029;
C.MINCONTSCALEFACTCD = 8030;
C.NGABORREPSCD = 8031;
C.NCONTRASTSCD = 8032;

% Ecodes that identify the trial parameters
% (Parameters that change on a trial to trial basis
% and so are dropped at the end of each trial)
C.DIAMCD = 7000;
C.SFCD = 7001;
C.TFCD = 7002;
C.LCONTCD = 7003;
C.MCONTCD = 7004;
C.SCONTCD = 7005;
C.ORIENTCD = 7006;
C.PHASECD = 7007;
C.NCYCLESCD = 7008;
C.NFRAMESCD = 7009;
C.PROTOCOLCD = 7010;
C.STIMTYPECD = 7011;

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
trialParamsHeader = {'fp_on', 'fp_acq', 'stim_on', 'stim_off', 'fp_off', 'diam', 'sf', 'tf', 'lcont', 'mcont', 'scont', 'orient', 'phase', 'ncycles', 'nframes', 'protocol', 'stimtype';...
                     'time', 'time', 'time', 'time', 'time', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'int',  'int',      'int';...
                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; ...
                     C.FPONCD, C.FPACQCD, C.STIMONCD, C.STIMOFFCD, C.FPOFFCD, C.DIAMCD, C.SFCD, C.TFCD, C.LCONTCD, C.MCONTCD, C.SCONTCD, C.ORIENTCD, C.PHASECD, C.NCYCLESCD, C.NFRAMESCD, C.PROTOCOLCD, C.STIMTYPECD};


%for the stro.sum.expParams array (includes things that get dropped only
%once):
exptParamsHeader = {'fp_x', 'fp_y', 'fp_size', 'fp_r', 'fp_g', 'fp_b', 'eyewin_x', 'eyewin_y', 'rf_x', 'rf_y', 'protocol', 'ntrlspercond', 'framerate', 'gamma_table', 'mon_spd', 'fundamentals', 'nsigmas', 'ncyclesplat', 'ncyclesramp', 'maxcontrastfactor', 'mincontrastfactor', 'ngaborreps', 'ngaborcontrasts'; ...
                    'int', 'int', 'int', 'int', 'int', 'int', 'int', 'int', 'int', 'int', 'int', 'int', 'double', 'double', 'double', 'double', 'float', 'float', 'float', 'float', 'float', 'int', 'int'; ...
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0; ...
                    C.FIXXCD, C.FIXYCD, C.FPSIZECD, C.FPRCD, C.FPGCD, C.FPBCD, C.EYEWINXCD, C.EYEWINYCD, C.RFXCD, C.RFYCD, C.INITPROTOCOLCD, C.NTRLSPERCONDCD, C.FRAMERATECD, C.GAMMATABLECD, C.MONSPDCD, C.FUNDAMENTALSCD, C.NSIGMASCD, C.NCYCLESPLATCD, C.NCYCLESRAMPCD, C.MAXCONTSCALEFACTCD, C.MINCONTSCALEFACTCD, C.NGABORREPSCD, C.NCONTRASTSCD};

additionalRasterFields = {};


%define how good and bad trials are to be distinguished
badLogic = 'find(trialCodes(:,2) == C.ABORTCD, 1)';


%define other routines as an array of function handles. the first row of
%the array should be a name of the handle, and the second row should be the
%actual function handle. The functions should be included as non-nested
%subfunctions in this file. DTStimCodes.m serves as a template
otherInstructions = {}; %none for now

end