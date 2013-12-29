function [C, trialParamsHeader, exptParamsHeader, additionalRasterFields, badLogic, otherInstructions] = DTCodesTraining()


% All the relavent codes and header info for the detection task returned to
% the global structure 'c'
%
% CAH 

C.VALOFFSET = 4000;


%***************************%
%  EXPERIMENT CODES         %
%***************************%
C.FPPOSXCD = 8000;
C.FPPOSYCD = 8001;
C.FPWINSIZECD = 8002;
C.FPSIZECD = 8003;
C.FLASHLENCD = 8004;
C.FLASHSIZECD = 8005;
C.FRAMEOFFSETCD = 8006;
C.FRAMEMODECD = 8007;
C.TARGDISTCD = 8008;
C.TARGSIZE = 8009;
C.RFCOLORLCD = 8010;
C.RFCOLORMCD = 8011;
C.RFCOLORSCD = 8012;
C.GAMMATABLECD = 8013;
C.MONSPECTCD = 8014;
C.MMTXCD = 8015;
C.BKGNDRCD = 8016;       %double (normalized intensity)
C.BKGNDGCD = 8017;       %double (normalized intensity)
C.BKGNDBCD = 8018;       %double (normalized intensity)
C.NCONTRASTSCD = 8019; %int
C.RFPOSXCD = 8020;   %int
C.RFPOSYCD = 8021;   %int

%***************************%
%       TRIAL CODES         %
%***************************%
C.CORRECTCD = 7000;  %int
C.FLASHCONTRASTCD = 7001;
C.FLASHRCD = 7002; %DOUBLE
C.FLASHGCD = 7003; %DOUBLE
C.FLASHBCD = 7004; %DOUBLE
C.FLASHPOSXCD = 7005; %int
C.FLASHPOSYCD = 7006; %int
C.GABORSIGMACD = 7007; %int
C.GABORLAMBDACD = 7008; %int
C.GABORSPEEDCD = 7009; %double
C.GABORTHETACD = 7010; %double
C.GABORBGAMMACD = 7011; %double
C.NUMFRAMESCD = 7012; %long
C.ACTRGUNCD = 7013; %DOUBLE
C.ACTGGUNCD = 7014;%DOUBLE
C.ACTBGUNCD =7015; %DOUBLE


%***************************%
%    REX TIMING CODES       %
%***************************%
C.TRLSTARTCD = 1000;
C.PAUSECD = 1003;
C.FPONCD = 1004;
C.TARGONCD = 1008;
C.ABORTCD = 1013;
C.ERRCD = 1014;
C.REWCD = 1015;
C.EOTCD = 1016;
C.FRAMEONCD = 1050;
C.STIMONCD = 1051;
C.DECISIONCD = 1052;
C.STIMONREPORTCD = 1053;



%***************************%
%    STUFF FOR NEX2STRO     %
%***************************%

% this stuff will define the indicies and unpacking of stro.trial
trialParamsHeader = { 'fp_on', 'flash_on', 'flash_on_rep', 'targ_on', 'choice_time', 'correct', 'flash_x', 'flash_y', 'flash_R', 'flash_G', 'flash_B', 'act_flash_R', 'act_flash_G', 'act_flash_B', 'gabor_sigma', 'gabor_lambda', 'gabor_theta', 'gabor_gamma', 'gabor_speed', 'numframes';...
                      'time', 'time', 'time', 'time', 'time', 'int', 'int', 'int', 'double', 'double', 'double', 'double', 'double', 'double', 'int', 'int', 'double', 'double', 'double', 'long';...
                      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
                      C.FPONCD, C.STIMONCD, C.STIMONREPORTCD, C.TARGONCD, C.DECISIONCD, C.CORRECTCD, C.FLASHPOSXCD, C.FLASHPOSYCD, C.FLASHRCD, C.FLASHGCD, C.FLASHBCD, C.ACTRGUNCD, C.ACTGGUNCD, C.ACTBGUNCD, C.GABORSIGMACD, C.GABORLAMBDACD, C.GABORTHETACD, C.GABORBGAMMACD, C.GABORSPEEDCD, C.NUMFRAMESCD};

%determines the indicies and unpacking of stro.sum.exptParams
exptParamsHeader = {'fp_posX', 'fp_posY', 'fp_winsize', 'fp_size', 'flash_length', 'flash_size', 'frame_offset', 'frame_mode', 'targ_dist', 'targ_size', 'RF_color_L', 'RF_color_M', 'RF_color_S', 'gamma_table', 'mon_spect', 'm_mtx', 'bkgnd_r', 'bkgnd_g', 'bkgnd_b', 'nContrasts', 'rf_x', 'rf_y'; ...
                    'int', 'int', 'int', 'int', 'int', 'int', 'int', 'int', 'double', 'int', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'int', 'int', 'int';...
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0;...
                    C.FPPOSXCD, C.FPPOSYCD, C.FPWINSIZECD, C.FPSIZECD, C.FLASHLENCD, C.FLASHSIZECD, C.FRAMEOFFSETCD, C.FRAMEMODECD, C.TARGDISTCD, C.TARGSIZE, C.RFCOLORLCD, C.RFCOLORMCD, C.RFCOLORSCD, C.GAMMATABLECD, C.MONSPECTCD, C.MMTXCD, C.BKGNDRCD, C.BKGNDGCD,C.BKGNDBCD, C.NCONTRASTSCD, C.RFPOSXCD, C.RFPOSYCD};
               

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
