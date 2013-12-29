function [C, trialParamsHeader, exptParamsHeader, additionalRasterFields, badLogic, otherInstructions] = DTScotCodes()
C.VALOFFSET = 4000;

% header codes
C.FPPOSXCD        = 8000; % ints
C.FPPOSYCD        = 8001;
C.FPWINSIZECD     = 8002;
C.FPSIZECD        = 8003;
C.STIMTIMEONCD    = 8004;
C.STIMSIZECD      = 8005;
C.TARGSIZECD      = 8006;
C.GAMMATABLECD    = 8007; % double, bookended
C.MONSPECTCD      = 8008; % double, bookended
C.BKGNDRCD        = 8009; % double
C.BKGNDGCD        = 8010; % double
C.BKGNDBCD        = 8011; % double
C.STIMPOSXCD      = 8012; % float
C.STIMPOSYCD      = 8013; % float
C.PIXPERDEGCD     = 8014; % double
C.REFRESHRATECD   = 8015; % double
C.STIMCOLORDIRSCD = 8016; % float, bookended
C.THRESHESTCD     = 8017; % float, bookended
C.QUESTSIGMACD    = 8018; % float
C.QUESTBETACD     = 8019; % float
C.QUESTRANGECD    = 8020; % double, bookended

% trial codes
C.CORRECTCD        = 7000; % int
C.STIMRGUNCD       = 7001; % float
C.STIMGGUNCD       = 7002; % float
C.STIMBGUNCD       = 7003; % float
C.STIMPOSACTUALXCD = 7004; % float
C.STIMPOSACTUALYCD = 7005; % float
C.NUMFRAMESCD      = 7006; % int
C.COLORDIRCD       = 7007; % int
C.STIMLOCCD        = 7008; % int
C.QUESTMODECD      = 7009; % float

% timing codes
C.TRLSTARTCD = 1000;
C.PAUSECD    = 1003;
C.FPONCD     = 1004;
C.TARGONCD   = 1008;
C.ABORTCD    = 1013;
C.ERRCD      = 1014;
C.REWCD      = 1015;
C.EOT        = 1016;
C.STIMONCD   = 1030;
C.STIMOFFCD  = 1031;
C.DECISIONCD = 1050;

% this stuff will define the indicies and unpacking of stro.trial
trialParamsHeader = {
    'fp_on'             'time'     0    C.FPONCD
    'stim_on'           'time'     0    C.STIMONCD
    'stim_off'          'time'     0    C.STIMOFFCD
    'targ_on'           'time'     0    C.TARGONCD
    'choice_time'       'time'     0    C.DECISIONCD
    'rew_time'          'time'     0    C.REWCD
    'correct'           'int'      0    C.CORRECTCD
    'stim_x'            'float'    0    C.STIMPOSACTUALXCD
    'stim_y'            'float'    0    C.STIMPOSACTUALYCD
    'stim_r'            'float'    0    C.STIMRGUNCD
    'stim_g'            'float'    0    C.STIMGGUNCD
    'stim_b'            'float'    0    C.STIMBGUNCD
    'numframes'         'int'      0    C.NUMFRAMESCD
    'color_dir'         'int'      0    C.COLORDIRCD
    'stim_loc'          'int'      0    C.STIMLOCCD
    'quest_mode'        'float'    0    C.QUESTMODECD
    }';

% determines the indicies and unpacking of stro.sum.exptParams
exptParamsHeader = {
    'fp_x'            'int'       0    C.FPPOSXCD
    'fp_y'            'int'       0    C.FPPOSYCD
    'fp_winsize'      'int'       0    C.FPWINSIZECD
    'fp_size'         'int'       0    C.FPSIZECD
    'stim_timeon'     'int'       0    C.STIMTIMEONCD
    'stim_size'       'int'       0    C.STIMSIZECD
    'targ_size'       'int'       0    C.TARGSIZECD
    'stim_colors'     'float'     1    C.STIMCOLORDIRSCD
    'gamma_table'     'double'    1    C.GAMMATABLECD
    'mon_spd'         'double'    1    C.MONSPECTCD
    'bkgnd_r'         'double'    0    C.BKGNDRCD
    'bkgnd_g'         'double'    0    C.BKGNDGCD
    'bkgnd_b'         'double'    0    C.BKGNDBCD
    'stim_x'          'float'     0    C.STIMPOSXCD
    'stim_y'          'float'     0    C.STIMPOSYCD
    'pixperdeg'       'double'    0    C.PIXPERDEGCD
    'refreshrate'     'double'    0    C.REFRESHRATECD
    'quest_sigma'     'float'     0    C.QUESTSIGMACD
    'quest_beta'      'float'     0    C.QUESTBETACD
    'quest_ranges'    'double'    1    C.QUESTRANGECD
    }';

% determines the additional fields to stro.ras (other than spikes and analog
% signals acquired online
additionalRasterFields = {};  % none for this paradigm yet

% deterimes how to identify bad trials
badLogic = 'find(trialCodes(:,2) == C.ABORTCD, 1)';

% define other routines as an array of function handles. the first row of
% the array should be a name of the handle, and the second row should be the
% actual function handle. The functions should be included as non-nested
% subfunctions in this file. DTStimCodes.m serves as a template
otherInstructions = {}; % none for now

end