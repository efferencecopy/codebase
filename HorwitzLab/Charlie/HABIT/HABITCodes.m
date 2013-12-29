function [C, trialParamsHeader, exptParamsHeader, additionalRasterFields, badLogic, otherInstructions] = HABITCodes()

% codes necessary to unpack HABIT data files
%
% CAH 12/09  Adding Codes


C.VALOFFSET = 4000;

%***************************%
%  EXPERIMENT CODES         %
%***************************%
C.GAMMATABLE =      8001; %DOUBLE 1
C.MONSPECT =        8002; %DOUBLE 1
C.MMTX =            8003; %DOUBLE 1
C.BKGNDRGB =        8004; %DOUBLE (gun intensity b/w 0&1) 1
C.PIXPERDEG =       8005; %DOUBLE 0
C.FRAMERATE =       8006; %DOUBLE 0
C.FPPOSX =          8007; %INT 0
C.FPPOSY =          8008; %INT 0
C.FPWINSIZE =       8009; %INT 0
C.FPSIZE =          8010; %INT 0
C.ALPHAGUESS =      8011; %FLOATS 1
C.GABORCOLOR =      8012; %FLOATS 1
C.GABORSIGMA =      8013; %FLOATS 0
C.GABORTHETA =      8014; %FLOATS 0
C.GABORSIZE =       8015; %FLOATS 0
C.GABORGAMMA =      8016; %FLOATS 0
C.GABORSPEED =      8017; %FLOATS 0
C.GABORLENGTH =     8018; %INT 0
C.RFPOSX =          8019; %FLOATS 0
C.RFPOSY =          8020; %FLOATS 0
C.HABITORIENT =     8021; %FLOATS 0
C.HABITSF =         8022; %FLOATS 0
C.HABITSPEED =      8023; %FLOATS 0
C.HABITLMS =        8024; %FLOATS 1
C.HABITLENGTH =     8025; %INT 0
C.TARGSIZE =        8026; %INT 0
C.T1DIST =          8027; %FLOATS 0
C.T2DIST =          8028; %FLOATS 0
C.TARGWINSIZE =     8029; %INT 0
C.QUESTSIGMA =      8030; %FLOATS 0
C.QUESTBETA =       8031; %FLOATS 0
C.QUESTDOMAIN =     8032; %DOUBLE 1

%***************************%
%       TRIAL ECODES        %
%***************************%
C.CORRECTCD =       7000; %INT
C.GABORRGUNCD =     7001; %LONG
C.GABORGGUNCD =     7002; %LONG
C.GABORBGUNCD =     7003; %LONG
C.ACTRGUNCD =       7004; %LONG
C.ACTGGUNCD =       7005; %LONG
C.ACTBGUNCD =       7006; %LONG
C.GABORPOSXCD =     7007; %FLOAT
C.GABORPOSYCD =     7008; %FLOAT
C.GABORLAMBDACD =   7009; %FLOAT
C.COLORDIRCD =      7010; %LONG
C.GABORFRAMESCD =   7011; %LONG
C.HABITFRAMESCD =   7012; %LONG


%***************************%
%   TRIAL TIMING CODES      %
%***************************%
C.TRIALSTARTCD =    1001;
C.FPONCD =          1002;
C.HABITONCD =       1003;
C.REPHABITONCD =    1004;
C.HABITLASTFRAMECD =1005;
C.GABORONCD =       1006;
C.REPGABORONCD =    1007;
C.GABORLASTFRAMECD =1008;
C.TARGONCD =        1009;
C.CHOICETIMECD =    1010;
C.EOTCD =           1011;
C.REWARDCD =        1012;
C.ABORTCD =         1013;
C.PAUSECD =         1014;
C.EOEXPTCD =        1015;


%***************************%
%    STUFF FOR NEX2STRO     %
%***************************%

% this stuff will define the indicies and unpacking of stro.trial
trialParamsHeader = {'trial_start', 'fp_on', 'habit_on', 'rep_habit_on', 'habit_lastframe', 'gabor_on', 'rep_gabor_on', 'gabor_lastframe', 'targ_on', 'choice_time', 'reward_time', 'trial_correct', 'R_gun', 'G_gun', 'B_gun', 'Act_R_gun', 'Act_G_gun', 'Act_B_gun', 'gabor_pos_x', 'gabor_pos_y', 'gabor_lambda', 'colordir', 'gabor_frames', 'habit_frames';...
                     'time',  'time',  'time',  'time',  'time',  'time',  'time',  'time',  'time',  'time',  'time', 'int', 'long', 'long', 'long', 'long', 'long', 'long', 'float', 'float', 'float', 'long', 'long', 'long';...
                     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0; ...
                     C.TRIALSTARTCD, C.FPONCD, C.HABITONCD, C.REPHABITONCD, C.HABITLASTFRAMECD, C.GABORONCD, C.REPGABORONCD, C.GABORLASTFRAMECD, C.TARGONCD, C.CHOICETIMECD, C.REWARDCD, C.CORRECTCD, C.GABORRGUNCD, C.GABORGGUNCD, C.GABORBGUNCD, C.ACTRGUNCD, C.ACTGGUNCD, C.ACTBGUNCD, C.GABORPOSXCD, C.GABORPOSYCD, C.GABORLAMBDACD, C.COLORDIRCD, C.GABORFRAMESCD, C.HABITFRAMESCD};

%determines the indicies and unpacking of stro.sum.exptParams
exptParamsHeader = {'gammaTable', 'monSpect', 'Mmtx', 'bkgndrgb', 'pixperdeg', 'frameRate', 'fp_posX', 'fp_posY', 'fp_winsize', 'fp_size', 'alphaGuess', 'gaborColor', 'gaborSigma', 'gaborTheta', 'gaborSize', 'gaborGamma', 'gaborSpeed', 'gaborLength','RFPosX', 'RFPosY', 'habitOrient', 'habitSF', 'habitSpeed', 'habitLMS', 'habitLength', 'targSize', 'T1dist', 'T2dist','targWinsize', 'questSigma', 'questBeta', 'questDomain';...
                    'double','double','double','double','double','double','int','int','int','int','float','float','float','float','float','float','float','int','float','float','float','float','float','float','int','int','float','float','int','float','float', 'double';...
                    1,1,1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1;...
                    C.GAMMATABLE, C.MONSPECT, C.MMTX, C.BKGNDRGB, C.PIXPERDEG, C.FRAMERATE, C.FPPOSX, C.FPPOSY, C.FPWINSIZE, C.FPSIZE, C.ALPHAGUESS, C.GABORCOLOR, C.GABORSIGMA, C.GABORTHETA, C.GABORSIZE, C.GABORGAMMA, C.GABORSPEED, C.GABORLENGTH, C.RFPOSX, C.RFPOSY, C.HABITORIENT, C.HABITSF, C.HABITSPEED, C.HABITLMS, C.HABITLENGTH, C.TARGSIZE, C.T1DIST, C.T2DIST, C.TARGWINSIZE, C.QUESTSIGMA, C.QUESTBETA, C.QUESTDOMAIN};
                

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
