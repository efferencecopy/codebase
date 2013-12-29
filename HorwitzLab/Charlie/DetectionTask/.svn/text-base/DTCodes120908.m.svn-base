function [C, trialParamsHeader, exptParamsHeader, additionalRasterFields, badLogic] = DTCodes();


% All the relavent codes and header info for the detection task returned to
% the global structure 'c'
%
% CAH 04/08
% CAH 07/08 did away with most of the training stuff and updated for use
%           with adaptive proceedures
% CAH 09/08 updated to work with floats natively. Added codes for frame
%           rates

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
C.RFPOSXCD =        8014;       %int
C.RFPOSYCD =        8015;       %int
C.PIXPERDEGCD =     8016;       %DOUBLE
C.EXPTMETHCD =      8017;       %int
C.FRAMERATECD =     8018;       %DOUBLE
C.COLORDIRSCD =     8019;       %float

%***************************%
%       TRIAL CODES         %
%***************************%
C.CORRECTCD =       7000;  %int
C.FLASHRCD =        7001; %long
C.FLASHGCD =        7002; %long
C.FLASHBCD =        7003; %long
C.FLASHPOSXCD =     7004; %double
C.FLASHPOSYCD =     7005; %double
C.GABORSIGMACD =    7006; %int
C.GABORLAMBDACD =   7007; %int
C.GABORSPEEDCD =    7008; %double
C.GABORTHETACD =    7009; %double
C.GABORBGAMMACD =   7010; %double
C.NUMFRAMESCD =     7011; %long
C.ACTRGUNCD =       7012; %long
C.ACTGGUNCD =       7013;%long
C.ACTBGUNCD =       7014; %long
C.COLORDIRCD =      7015; %long
C.CNTRSTLEVCD =     7016; %long


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



%***************************%
%    STUFF FOR NEX2STRO     %
%***************************%

% this stuff will define the indicies and unpacking of stro.trial
trialParamsHeader = { 'fp_on', 'flash_on', 'rep_flash_on', 'targ_on', 'choice_time', 'correct', 'flash_x', 'flash_y', 'flash_R', 'flash_G', 'flash_B', 'act_flash_R', 'act_flash_G', 'act_flash_B', 'gabor_sigma', 'gabor_lambda', 'gabor_theta', 'gabor_gamma', 'gabor_speed', 'numframes', 'color_dir', 'cntrst_lev';...
                      'time', 'time', 'time', 'time', 'time', 'int', 'double', 'double', 'long', 'long', 'long', 'long', 'long', 'long', 'int', 'double', 'double', 'double', 'double', 'long', 'long', 'long';...
                      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
                      C.FPONCD, C.STIMONCD, C.STIMONREPORTCD, C.TARGONCD, C.DECISIONCD, C.CORRECTCD, C.FLASHPOSXCD, C.FLASHPOSYCD, C.FLASHRCD, C.FLASHGCD, C.FLASHBCD, C.ACTRGUNCD, C.ACTGGUNCD, C.ACTBGUNCD, C.GABORSIGMACD, C.GABORLAMBDACD, C.GABORTHETACD, C.GABORBGAMMACD, C.GABORSPEEDCD, C.NUMFRAMESCD, C.COLORDIRCD, C.CNTRSTLEVCD};

%determines the indicies and unpacking of stro.sum.exptParams
exptParamsHeader = {'fp_posX', 'fp_posY', 'fp_winsize', 'fp_size', 'flash_length', 'flash_size', 'frame_offset', 'targ_size', 'RF_colors', 'gamma_table', 'mon_spect', 'm_mtx', 'bkgnd_r', 'bkgnd_g', 'bkgnd_b', 'rf_x', 'rf_y', 'pixperdeg', 'expt_meth', 'frame_rate'; ...
                    'int', 'int', 'int', 'int', 'int', 'int', 'int', 'int', 'float', 'double', 'double', 'double', 'double', 'double', 'double', 'int', 'int', 'double', 'int', 'double';...
                    0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0;...
                    C.FPPOSXCD, C.FPPOSYCD, C.FPWINSIZECD, C.FPSIZECD, C.FLASHLENCD, C.FLASHSIZECD, C.FRAMEOFFSETCD, C.TARGSIZECD, C.COLORDIRSCD, C.GAMMATABLECD, C.MONSPECTCD, C.MMTXCD, C.BKGNDRCD, C.BKGNDGCD, C.BKGNDBCD, C.RFPOSXCD, C.RFPOSYCD, C.PIXPERDEGCD, C.EXPTMETHCD, C.FRAMERATECD};
                

%determines the additional fields to stro.ras (other than spikes and analog
%signals acquired online
additionalRasterFields = {};  %none for this paradigm yet


%deterimes how to identify bad trials
badLogic = 'find(trialCodes(:,2) == C.ABORTCD, 1)';


