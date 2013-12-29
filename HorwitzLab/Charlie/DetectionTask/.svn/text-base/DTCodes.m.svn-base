function [C, trialParamsHeader, exptParamsHeader, additionalRasterFields, badLogic, otherInstructions] = DTCodes()
% All the relavent codes and header info for the detection task returned to
% the global structure 'c'
%
% CAH 04/08
% CAH 07/08 did away with most of the training stuff and updated for use
%           with adaptive proceedures
% CAH 09/08 updated to work with floats natively. Added codes for frame rates
% CAH 12/08 changed some of the cast types for better integration with
%           online gabor fitting
% CAH 03/09 added codes for QUEST

C.VALOFFSET = 4000;

%***************************%
%  EXPERIMENT CODES         %
%***************************%
C.FPPOSXCD =        8000;
C.FPPOSYCD =        8001;
C.FPWINSIZECD =     8002;
C.FPSIZECD =        8003;
C.FLASHLENCD =      8004;
C.FLASHSIZECD =     8005;
C.FRAMEOFFSETCD =   8006;
C.TARGSIZECD =      8007;
C.GAMMATABLECD =    8008;
C.MONSPECTCD =      8009;
C.MMTXCD =          8010;
C.BKGNDRCD =        8011;       %double (normalized intensity)
C.BKGNDGCD =        8012;       %double (normalized intensity)
C.BKGNDBCD =        8013;       %double (normalized intensity)
C.RFPOSXCD =        8014;       %float
C.RFPOSYCD =        8015;       %float
C.PIXPERDEGCD =     8016;       %DOUBLE
C.EXPTMETHCD =      8017;       %int
C.FRAMERATECD =     8018;       %DOUBLE
C.COLORDIRSCD =     8019;       %float
C.LOWSCALARSCD =    8020;       %floats (the number that each LMS get multiplied by)
C.NCONTRASTSCD =    8021;       %int
C.QUESTSIGMACD =    8022;       %float
C.QUESTBETACD =     8023;       %float
C.QUESTRANGECD =    8024;       %doubles (1x9)

%***************************%
%       TRIAL CODES         %
%***************************%
C.CORRECTCD =       7000; %int
C.FLASHRCD =        7001; %long
C.FLASHGCD =        7002; %long
C.FLASHBCD =        7003; %long
C.FLASHPOSXCD =     7004; %float
C.FLASHPOSYCD =     7005; %float
C.GABORSIGMACD =    7006; %float
C.GABORLAMBDACD =   7007; %float
C.GABORSPEEDCD =    7008; %float
C.GABORTHETACD =    7009; %float
C.GABORBGAMMACD =   7010; %float
C.NUMFRAMESCD =     7011; %long
C.ACTRGUNCD =       7012; %long
C.ACTGGUNCD =       7013; %long
C.ACTBGUNCD =       7014; %long
C.COLORDIRCD =      7015; %long
C.CNTRSTLEVCD =     7016; %long
C.QUESTTRLTHRESHCD =7017; %float
C.TRIALTYPECD =     7018; %int   0=> detection, 1=>fixation
C.FRAMESPRESENTCD = 7019; %int

%***************************%
%    REX TIMING CODES       %
%***************************%
C.TRLSTARTCD =      1000;
C.PAUSECD =         1003;
C.FPONCD =          1004;
C.TARGONCD =        1008;
C.ABORTCD =         1013;
C.ERRCD =           1014;
C.REWCD =           1015;
C.EOTCD =           1016;
C.FRAMEONCD =       1050;
C.STIMONCD =        1051;
C.DECISIONCD =      1052;
C.STIMONREPORTCD =  1053;
C.STIMOFFCD =       1054;

%***************************%
%    STUFF FOR NEX2STRO     %
%***************************%

% this stuff will define the indicies and unpacking of stro.trial
trialParamsHeader = {
    'trial_type'        'int'      0    C.TRIALTYPECD
    'fp_on'             'time'     0    C.FPONCD
    'frames_present'    'int'      0    C.FRAMESPRESENTCD
    'frame_on'          'time'     0    C.FRAMEONCD
    'flash_on'          'time'     0    C.STIMONCD
    'rep_flash_on'      'time'     0    C.STIMONREPORTCD
    'flash_off'         'time'     0    C.STIMOFFCD
    'targ_on'           'time'     0    C.TARGONCD
    'choice_time'       'time'     0    C.DECISIONCD
    'rew_time'          'time'     0    C.REWCD
    'correct'           'int'      0    C.CORRECTCD
    'flash_x'           'float'    0    C.FLASHPOSXCD
    'flash_y'           'float'    0    C.FLASHPOSYCD
    'flash_R'           'long'     0    C.FLASHRCD
    'flash_G'           'long'     0    C.FLASHGCD
    'flash_B'           'long'     0    C.FLASHBCD
    'act_flash_R'       'long'     0    C.ACTRGUNCD
    'act_flash_G'       'long'     0    C.ACTGGUNCD
    'act_flash_B'       'long'     0    C.ACTBGUNCD
    'gabor_sigma'       'float'    0    C.GABORSIGMACD
    'gabor_lambda'      'float'    0    C.GABORLAMBDACD
    'gabor_theta'       'float'    0    C.GABORTHETACD
    'gabor_gamma'       'float'    0    C.GABORBGAMMACD
    'gabor_speed'       'float'    0    C.GABORSPEEDCD
    'numframes'         'long'     0    C.NUMFRAMESCD
    'color_dir'         'long'     0    C.COLORDIRCD
    'cntrst_lev'        'long'     0    C.CNTRSTLEVCD
    'quest_thresh'      'float'    0    C.QUESTTRLTHRESHCD
    }';

%determines the indicies and unpacking of stro.sum.exptParams
exptParamsHeader = {
    'fp_posX'         'int'       0    C.FPPOSXCD
    'fp_posY'         'int'       0    C.FPPOSYCD
    'fp_winsize'      'int'       0    C.FPWINSIZECD
    'fp_size'         'int'       0    C.FPSIZECD
    'flash_length'    'int'       0    C.FLASHLENCD
    'flash_size'      'int'       0    C.FLASHSIZECD
    'frame_offset'    'int'       0    C.FRAMEOFFSETCD
    'targ_size'       'int'       0    C.TARGSIZECD
    'RF_colors'       'float'     1    C.COLORDIRSCD
    'low_scalars'     'float'     1    C.LOWSCALARSCD
    'nContrasts'      'int'       0    C.NCONTRASTSCD
    'gamma_table'     'double'    1    C.GAMMATABLECD
    'mon_spect'       'double'    1    C.MONSPECTCD
    'm_mtx'           'double'    1    C.MMTXCD
    'bkgnd_r'         'double'    0    C.BKGNDRCD
    'bkgnd_g'         'double'    0    C.BKGNDGCD
    'bkgnd_b'         'double'    0    C.BKGNDBCD
    'rf_x'            'float'     0    C.RFPOSXCD
    'rf_y'            'float'     0    C.RFPOSYCD
    'pixperdeg'       'double'    0    C.PIXPERDEGCD
    'expt_meth'       'int'       0    C.EXPTMETHCD
    'frame_rate'      'double'    0    C.FRAMERATECD
    'quest_sigma'     'float'     0    C.QUESTSIGMACD
    'quest_beta'      'float'     0    C.QUESTBETACD
    'quest_ranges'    'double'    1    C.QUESTRANGECD
    }';

%determines the additional fields to stro.ras (other than spikes and analog
%signals acquired online
additionalRasterFields = {};  %none for this paradigm yet

%deterimes how to identify bad trials
badLogic = 'find(trialCodes(:,2) == C.ABORTCD, 1)';

%define other routines as an array of function handles. the first row of
%the array should be a name of the handle, and the second row should be the
%actual function handle. The functions should be included as non-nested
%subfunctions in this file. DTStimCodes.m serves as a template
otherInstructions = {}; %none for now

end