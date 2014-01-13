function [C, trialParamsHeader, exptParamsHeader, additionalRasterFields, badLogic, otherInstructions] = NM2SCodes()

% This is a collection of constant definitions for use with analysis
% programs.  These codes are specific to the "nonmatch to sample" (NM2S) paradigm

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
C.TARGSIZECD = 8006;
C.TARGRCD = 8007;
C.TARGGCD = 8008;
C.TARGBCD = 8009;
C.TARGDISTPROPCD = 8010;
C.FIXWINSIZECD = 8011;
C.TARGWINSIZECD = 8012;
C.RFXCD = 8013;
C.RFYCD = 8014;
C.BKGNDRGBCD = 8015; %double,1
C.GAMMATABLECD = 8016; %double,1
C.MONSPDCD = 8017; %double,1
C.FRAMERATECD = 8018; %double, 1
C.FUNDAMENTALSCD = 8019; %double, 1
C.TRLSPERCONDCD = 8020;
C.NFRAMESSAMPLECD = 8021;
C.NFRAMESDISCCD = 8022;
C.EPCOMPCD = 8023;
C.NEPSAMPAVGCD = 8024;
C.THETACD = 8025;
C.LAMBDACD = 8026;
C.PHICD = 8027;
C.XOFFSETCD = 8028;
C.YOFFSETCD = 8029;

% Ecodes that identify the trial parameters
% (Parameters that change on a trial to trial basis
% and so are dropped at the end of each trial)
C.SAMPLELCD = 7000;
C.SAMPLEMCD = 7001;
C.SAMPLESCD = 7002;
C.T1LCD = 7003;
C.T1MCD = 7004;
C.T1SCD = 7005;
C.T2LCD = 7006;
C.T2MCD = 7007;
C.T2SCD = 7008;
C.GRAD1LCD = 7009;
C.GRAD1MCD = 7010;
C.GRAD1SCD = 7011;
C.GRAD2LCD = 7012;
C.GRAD2MCD = 7013;
C.GRAD2SCD = 7014;
C.GRADPOSCD = 7015;
C.NONMATCHSIDECD = 7016;
C.CORRECTCD = 7017;
C.MSTIMCD = 7018;
C.STIMCONFIGCD = 7019;

% Ecodes that are used to mark the times of various trial events
C.FPONCD = 1004;
C.FPACQCD = 1032;
C.SAMPLEONCD = 1010;
C.SAMPLEOFFCD = 1011;
C.STIMONCD = 1030;
C.STIMOFFCD = 1031;
C.MSTIMONCD = 1033;
C.FPOFFCD = 1006;
C.TARGONCD = 1008;
C.SACMADCD = 1026;
C.ERRCD = 1014;
C.REWCD = 1015;
C.EOTCD = 1016;
C.ABORTCD = 1013;

%******************%
%   for nex2stro   %
%******************%

% includes the trial column identifiers, the ecode id's, the data type, and
% bookending info.
trialParamsHeader = {'fp_on', 'fp_acq', 'sample_on','sample_off', 'stim_on', 'stim_off', 'mstim_on', 'fp_off', 'saccade', 'reward', 'T1_Lcc', 'T1_Mcc', 'T1_Scc', 'T2_Lcc', 'T2_Mcc', 'T2_Scc', 'Grad1_Lcc', 'Grad1_Mcc', 'Grad1_Scc', 'Grad2_Lcc', 'Grad2_Mcc', 'Grad2_Scc', 'mstim', 'gradient_pos', 'nonmatchside', 'correct', 'stimconfig';...
                     'time',  'time',   'time',     'time',       'time',    'time',     'time',     'time',   'time',    'time',   'float',  'float',  'float',  'float',  'float',  'float',  'float',     'float',     'float',     'float',     'float',     'float',     'int'    'float',        'int',          'int',     'int';...
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;  ...
                       C.FPONCD, C.FPACQCD, C.SAMPLEONCD, C.SAMPLEOFFCD, C.STIMONCD, C.STIMOFFCD, C.MSTIMONCD, C.FPOFFCD, C.SACMADCD, C.REWCD, C.T1LCD, C.T1MCD, C.T1SCD, C.T2LCD, C.T2MCD, C.T2SCD, C.GRAD1LCD, C.GRAD1MCD, C.GRAD1SCD, C.GRAD2LCD, C.GRAD2MCD, C.GRAD2SCD, C.MSTIMCD, C.GRADPOSCD, C.NONMATCHSIDECD, C.CORRECTCD, C.STIMCONFIGCD};

%for the stro.sum.expParams array (includes things that get dropped only
%once):
exptParamsHeader = {'fp_x', 'fp_y', 'fp_size', 'fp_r', 'fp_g', 'fp_b', 'targsize', 'targ_r', 'targ_g', 'targ_b', 'targdistprop', 'fixwinsize', 'targwinsize', 'rf_x', 'rf_y', 'bkgndrgb', 'gamma_table', 'mon_spd', 'framerate', 'fundamentals', 'trialspercond', 'nframes_sample','nframes_discriminanda', 'ep_compensation', 'nepsampavg', 'theta', 'lambda', 'phi', 'xoffset', 'yoffset';...
                    'int',  'int',  'int',     'int',  'int',  'int',  'int',      'int',    'int',    'int',    'float',        'int',        'int',         'int',  'int',  'double',   'double',      'double',  'double',    'double',       'int',           'int',           'int',                   'int',             'int',        'float', 'float',  'float','float',  'float';...
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;...
                    C.FIXXCD, C.FIXYCD, C.FPSIZECD, C.FPRCD, C.FPGCD, C.FPBCD, C.TARGSIZECD, C.TARGRCD, C.TARGGCD, C.TARGBCD, C.TARGDISTPROPCD, C.FIXWINSIZECD, C.TARGWINSIZECD, C.RFXCD, C.RFYCD, C.BKGNDRGBCD, C.GAMMATABLECD, C.MONSPDCD, C.FRAMERATECD, C.FUNDAMENTALSCD, C.TRLSPERCONDCD, C.NFRAMESSAMPLECD, C.NFRAMESDISCCD, C.EPCOMPCD, C.NEPSAMPAVGCD, C.THETACD, C.LAMBDACD, C.PHICD, C.XOFFSETCD, C.YOFFSETCD};

additionalRasterFields = {};


%define how good and bad trials are to be distinguished
badLogic = 'find(trialCodes(:,2) == C.ABORTCD, 1)';

%define other routines as an array of function handles. the first row of
%the array should be a name of the handle, and the second row should be the
%actual function handle. The functions should be included as non-nested
%subfunctions in this file. DTStimCodes.m serves as a template
otherInstructions = {}; %none for now

end