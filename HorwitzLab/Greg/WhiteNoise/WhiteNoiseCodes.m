function [C, trialParamsHeader, exptParamsHeader, additionalRasterFields, badLogic, otherInstructions] = WhiteNoiseCodes()

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
C.NPIXPERSTIXCD = 8010;
C.NSTIXPERSIDECD = 8011;
C.GAMMATABLECD = 8012; %double,1
C.MONSPDCD = 8013; %double,1
C.FRAMERATECD = 8014; %double, 1
C.GAUSSLOCUTCD = 8015;
C.GAUSSHICUTCD = 8016;
C.FUNDAMENTALSCD = 8017; %double, 1

% Ecodes that identify the trial parameters
% (Parameters that change on a trial to trial basis
% and so are dropped at the end of each trial)
C.SEEDCD = 7000;
C.MU1CD = 7001;
C.MU2CD = 7002;
C.MU3CD = 7003;
C.SIGMA1CD = 7004;
C.SIGMA2CD = 7005;
C.SIGMA3CD = 7006;
C.NFRAMESCD = 7007;
C.BKGNDRCD = 7008;
C.BKGNDGCD = 7009;
C.BKGNDBCD = 7010;
C.SYNTHIMAGECD = 7011;
C.NOISETYPECD = 7012;

% Ecodes that are used to mark the times of various trial events
C.FPONCD = 1004;
C.STIMONCD = 1030;
C.STIMOFFCD = 1031;
C.ALLOFFCD = 1012;
C.ERRCD = 1014;
C.REWCD = 1015;
C.EOTCD = 1016;
C.ABORTCD = 1013; %CAH

%******************%
%   for nex2stro   %
%******************%

% includes the trial column identifiers, the ecode id's, the data type, and
% bookending info.
trialParamsHeader = {'fp_on', 'stim_on', 'stim_off', 'all_off', 'seed', 'mu1', 'mu2', 'mu3', 'sigma1', 'sigma2', 'sigma3', 'num_frames', 'bkgnd_r', 'bkgnd_g', 'bkgnd_b', 'noise_type';...
                     'time', 'time', 'time', 'time', 'long', 'int', 'int', 'int', 'int', 'int', 'int', 'int', 'int', 'int', 'int', 'int';...
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; ...
                       C.FPONCD, C.STIMONCD, C.STIMOFFCD, C.ALLOFFCD, C.SEEDCD, C.MU1CD, C.MU2CD, C.MU3CD, C.SIGMA1CD, C.SIGMA2CD, C.SIGMA3CD, C.NFRAMESCD, C.BKGNDRCD, C.BKGNDGCD, C.BKGNDBCD, C.NOISETYPECD};


%for the stro.sum.expParams array (includes things that get dropped only
%once):
exptParamsHeader = {'fp_x', 'fp_y', 'fp_size', 'fp_r', 'fp_g', 'fp_b', 'eyewin_x', 'eyewin_y', 'rf_x', 'rf_y', 'npixperstix', 'nstixperside', 'gamma_table', 'mon_spd', 'fundamentals', 'framerate', 'gauss_locut', 'gauss_hicut'; ...
                    'int', 'int', 'int', 'int', 'int', 'int', 'int', 'int', 'int', 'int', 'int', 'int', 'double', 'double', 'double', 'double', 'int', 'int'; ...
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0; ...
                    C.FIXXCD, C.FIXYCD, C.FPSIZECD, C.FPRCD, C.FPGCD, C.FPBCD, C.EYEWINXCD, C.EYEWINYCD, C.RFXCD, C.RFYCD, C.NPIXPERSTIXCD, C.NSTIXPERSIDECD, C.GAMMATABLECD, C.MONSPDCD, C.FUNDAMENTALSCD, C.FRAMERATECD, C.GAUSSLOCUTCD, C.GAUSSHICUTCD};

additionalRasterFields = {'synth_image' ; ...
                          'long' ; ...
                          1 ; ...
                          C.SYNTHIMAGECD};


%define how good and bad trials are to be distinguished
badLogic = 'find(trialCodes(:,2) == C.ABORTCD, 1)';

%define other routines as an array of function handles. the first row of
%the array should be a name of the handle, and the second row should be the
%actual function handle. The functions should be included as non-nested
%subfunctions in this file. DTStimCodes.m serves as a template
otherInstructions = {}; %none for now

end




