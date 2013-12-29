%%%%%%%%%%%%%%%%%
% HeaderCodes.m %
%%%%%%%%%%%%%%%%%
% This is a collection of constant definitions for use with analysis
% programs.

% The offset that is added to the value ecodes that follow the
% identifier ecodes (7000-9000 range).  Note, this value must be
% identical to the definition of VALODFFSET in rig1.h on the REX machine.
VALOFFSET = 4000;

% Ecodes that identify the header parameters
% (Parameters that do not change on a trial to trial basis
% and so are dropped once, at the begining of the run)
FIXXCD = 8000;
FIXYCD = 8001;
EYEWINXCD = 8002;
EYEWINYCD = 8003;
BKGNDRGBCD = 8004;
DIMRGBCD = 8005;
FPSIZECD = 8006;
ABORTENBLCD = 8007;
MONKSTARTCD = 8008;

% Ecodes that identify the trial parameters
% (Parameters that change on a trial to trial basis
% and so are dropped at the end of each trial)
CORRECTCD = 7000;

% Ecodes that are used to mark the times of various trial events
PAUSECD = 1003;
FPONCD = 1004;
FPDIMCD = 1005;
ERRCD = 1014;
REWCD = 1015;
EOTCD = 1016;
TOUCHCD = 1024;