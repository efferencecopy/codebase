function [C, trialParamsHeader, exptParamsHeader, additionalRasterFields, badLogic, otherInstructions] = DTStimCodes()


% All the relavent codes and header info for the optical stimulation
% version of the detection task returned to the global structure 'c'
%
% CAH 05/11

C.VALOFFSET = 4000;


%***************************%
%  EXPERIMENT CODES         %
%***************************%
C.FPPOSXCD =				8000; %int
C.FPPOSYCD =				8001; %int
C.FPWINSIZECD =				8002; %int
C.FPSIZECD =				8003; %int
C.RFPOSXCD =				8004; %float
C.RFPOSYCD =				8005; %float
C.GABORLENGTHCD =			8006; %int
C.GABORGAMMACD =			8007; %float
C.GABORSIGMACD =			8008; %float
C.GABORTHETACD =			8009; %float
C.GABORSPEEDCD =			8010; %float
C.GABORSIZECD =				8011; %int
C.STIMCOLORCD =				8012; %float
C.LOWSCALARCD =				8013; %float
C.NCONTRASTSCD =			8014; %int
C.TARGSIZECD =				8015; %int
C.T1DISTCD =				8016; %float
C.T2DISTCD =				8017; %float
C.GAMMATABLECD =			8018; %double
C.MONSPECTCD =				8019; %double
C.MMTXCD =					8020; %double
C.BKGNDRCD =				8021; %double (normalized intensity)
C.BKGNDGCD =				8022; %double (normalized intensity)
C.BKGNDBCD =				8023; %double (normalized intensity)
C.PIXPERDEGCD =				8024; %double
C.FRAMERATECD =				8025; %double
C.BEHAVIORTYPECD =			8026; %int
C.TARGRGBCD =               8027; %long (1x3);
C.PARADIGMIDCD =			8999;


%***************************%
%       TRIAL CODES         %
%***************************%
C.CORRECTCD =       		7000; %int
C.GABORRGUNCD = 			7001; %long
C.GABORGGUNCD = 			7002; %long
C.GABORBGUNCD = 			7003; %long
C.ACTRGUNCD =           	7004; %long
C.ACTGGUNCD =               7005; %long
C.ACTBGUNCD =           	7006; %long
C.GABORPOSXCD = 			7007; %float
C.GABORPOSYCD =        		7008; %float
C.GABORLAMBDACD =           7009; %float
C.NUMFRAMESCD = 			7010; %long
C.COLORDIRCD =          	7011; %long
C.CNTRSTLEVCD = 			7012; %long
C.LASERWIDTHCD =            7013; %int
C.LASERINTPULSEINTCD =  	7014; %int
C.PREGABORDELAYCD = 		7015; %int
C.PREGABORLASERCD = 		7016; %int
C.POSTGABORLASERCD =    	7017; %int
C.LASERCATCHTRIALCD =       7018; %int  (i.e., fixation trials)
C.LASERTRIALCD =            7019; %int  (laser present as part of behavioral trial)



%***************************%
%    REX TIMING CODES       %
%***************************%
C.TRLSTARTCD = 				1000;
C.PAUSECD = 				1001;
C.FPONCD = 					1002;
C.GABORONCD = 				1003;
C.REPORTEDGABORONCD = 		1004;
C.GABOROFFCD = 				1005;
C.LASERONCD = 				1006;
C.LASEROFFCD = 				1007;
C.TARGONCD = 				1008;
C.REPORTEDTARGONCD =        1009;
C.TARGACQUIREDCD = 			1010;
C.ABORTCD = 				1011;
C.REWCD = 					1012;
C.EOTCD = 					1013;


%***************************%
%    STUFF FOR NEX2STRO     %
%***************************%

% this stuff will define the indicies druing construction of stro.trial.
% I'm constructing this as an <N x 4> cell array for clarity, but I need to
% transpose the array at the end to make the dimensionality agree with
% nex2stro
trialParamsHeader = {'trialStart',       'time',  0, C.TRLSTARTCD
                     'fpOn',             'time',  0, C.FPONCD
                     'gaborOn',          'time',  0, C.GABORONCD
                     'repGaborOn',       'time',  0, C.REPORTEDGABORONCD
                     'gaborOff',         'time',  0, C.GABOROFFCD
                     'targOn',           'time',  0, C.TARGONCD
                     'repTargOn',        'time',  0, C.REPORTEDTARGONCD
                     'targAcq',          'time',  0, C.TARGACQUIREDCD
                     'rewOn',            'time',  0, C.REWCD
                     'correct',          'int',   0, C.CORRECTCD
                     'Rgun',             'long',  0, C.GABORRGUNCD
                     'Ggun',             'long',  0, C.GABORGGUNCD
                     'Bgun',             'long',  0, C.GABORBGUNCD
                     'actRgun',          'long',  0, C.ACTRGUNCD
                     'actGgun',          'long',  0, C.ACTGGUNCD
                     'actBgun',          'long',  0, C.ACTBGUNCD
                     'gaborPosX',        'float', 0, C.GABORPOSXCD
                     'gaborPosY',        'float', 0, C.GABORPOSYCD
                     'gaborLambda',      'float', 0, C.GABORLAMBDACD
                     'numFrames',        'long',  0, C.NUMFRAMESCD
                     'colordir',         'long',  0, C.COLORDIRCD
                     'cntrstLev',        'long',  0, C.CNTRSTLEVCD
                     'laserWidth',       'int',   0, C.LASERWIDTHCD
                     'laserIPI',         'int',   0, C.LASERINTPULSEINTCD
                     'preGaborDelay',    'int',   0, C.PREGABORDELAYCD
                     'preGaborLaser',    'int',   0, C.PREGABORLASERCD
                     'postGaborLaser',   'int',   0, C.POSTGABORLASERCD
                     'laserCatchTrial',  'int',   0, C.LASERCATCHTRIALCD
                     'laserTrial',       'int',   0, C.LASERTRIALCD
                     }';
                 
%determines the indicies and unpacking of stro.sum.exptParams
exptParamsHeader = {'fixPosX',      'int',      0, C.FPPOSXCD
                    'fixPosy',      'int',      0, C.FPPOSYCD
                    'fixWinSize',   'int',      0, C.FPWINSIZECD
                    'fixSize',      'int',      0, C.FPSIZECD
                    'rfPosX',       'float',    0, C.RFPOSXCD
                    'rfPosY',       'float',    0, C.RFPOSYCD
                    'gaborLength',  'int',      0, C.GABORLENGTHCD
                    'gaborGamma',   'float',    0, C.GABORGAMMACD
                    'gaborSigma',   'float',    0, C.GABORSIGMACD
                    'gaborTheta',   'float',    0, C.GABORTHETACD
                    'gaborSpeed',   'float',    0, C.GABORSPEEDCD
                    'gaborSize',    'int',      0, C.GABORSIZECD
                    'gaborColors',  'float',    1, C.STIMCOLORCD
                    'lowScalars',   'float',    1, C.LOWSCALARCD
                    'nContrasts',   'int',      0, C.NCONTRASTSCD
                    'targSize',     'int',      0, C.TARGSIZECD
                    'targRGB',      'long',     1, C.TARGRGBCD
                    'T1dist',       'float',    0, C.T1DISTCD
                    'T2dist',       'float',    0, C.T2DISTCD
                    'gammaTable',   'double',   1, C.GAMMATABLECD
                    'monSpect',     'double',   1, C.MONSPECTCD
                    'mMtx',         'double',   1, C.MMTXCD
                    'bkgndr',       'double',   0, C.BKGNDRCD
                    'bkgndg',       'double',   0, C.BKGNDGCD
                    'bkgndb',       'double',   0, C.BKGNDBCD
                    'pixperdeg',    'double',   0, C.PIXPERDEGCD
                    'frameRate',    'double',   0, C.FRAMERATECD
                    'behaviorType', 'int',      0, C.BEHAVIORTYPECD
                    }';

%determines the additional fields to stro.ras (other than spikes and analog
%signals acquired online
additionalRasterFields = {};  %none for this paradigm yet


%deterimes how to identify bad trials
badLogic = 'find(trialCodes(:,2) == C.ABORTCD, 1)';

%define other routines as an array of function handles. the first row of
%the array should be a name of the handle, and the second row should be the
%actual function handle. The functions should be included as non-nested
%subfunctions in this file.
otherInstructions = {'laserOnTime', 'laserOffTime';...
                     @getLaserOnset, @getLaserOffset};

end

%gets the onset time of each laser pulse
function laserOnTimes = getLaserOnset(trialCodes, C)
    list = trialCodes(:,2) == C.LASERONCD;
    if any(list)
        laserOnTimes = trialCodes(list,1);
    else
        laserOnTimes = [];
    end
end


%gets the offset time of each laser pulse (if present)
function laserOffTimes = getLaserOffset(trialCodes, C)
    list = trialCodes(:,2) == C.LASEROFFCD;
    if any(list)
        laserOffTimes = trialCodes(list,1);
    else
        laserOffTimes = [];
    end
end






