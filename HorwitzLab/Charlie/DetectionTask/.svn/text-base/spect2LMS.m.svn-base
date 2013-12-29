function triplets = spect2LMS

%unpacks a bunch of spectral data from previous measurments using the PR650
%and then outputs the cone contrast triplets for each of the spectra. I'm
%using a function b/c it allows for subfunctions...
%
% CAH 10/09

%find the appropriate path and cal structure
presDir = pwd;
if isunix
    cd('/Users/charliehass/LabStuff/Data/DTcalibrationData')
    calData = load('/Users/charliehass/Desktop/MatLab Code/Slave/Monitor Calibration/Monitor data/Dell 4/Dell4BitsCal.mat');
elseif ispc
    cd('C:\Documents and Settings\Administrator\My Documents\Rex Code\DTcalibrationData')
    calData = load('C:\Documents and Settings\Administrator\Desktop\MatlabCode\Slave\Monitor Calibration\Monitor data\Dell 4\Dell4BitsCal.mat');
else
    error('path is unspecified')
end


%load in the fundamentals:
t = load ('T_cones_smj.mat');
fundamentals = t.T_cones_smj';
calData = calData.cals{end};
monspd = SplineSpd(calData.S_device, calData.P_device, t.S_cones_smj);

% % data from the center of the screen
% disp('  CENTER OF THE SCREEN:')
% bkgndspd = load('bkgndCenter.mat');
% lspd = load('lisoCenter.mat');
% mspd = load('misoCenter.mat');
% sspd = load('sisoCenter.mat');

% % data from the Lower Left quadrant
% disp('    LOWER LEFT QUAD')
% bkgndspd = load('bkgndLowLeft.mat');
% lspd = load('lisoLowLeft.mat');
% mspd = load('misoLowLeft.mat');
% sspd = load('sisoLowLeft.mat');

% data from the Upper Right quadrant
disp('    UPPER RIGHT QUAD')
bkgndspd = load('bkgndUpRight.mat');
lspd = load('lisoUpRight.mat');
mspd = load('misoUpRight.mat');
sspd = load('sisoUpRight.mat');

%print out the triplet values to the workspace:
triplets = [];
spd2lms(mean(bkgndspd.spd,1), mean(lspd.spd,1), 'L_iso');
spd2lms(mean(bkgndspd.spd,1), mean(mspd.spd,1), 'M_iso');
spd2lms(mean(bkgndspd.spd,1), mean(sspd.spd,1), 'S_iso');

%be nice and cd back to the original dir
cd(presDir);


    function spd2lms(bkgndspd, testspd, description);
        bkgndspd = SplineSpd([380, 4, 101], bkgndspd(:), t.S_cones_smj);
        testspd = SplineSpd([380, 4, 101], testspd(:), t.S_cones_smj);
        bkgndlms = bkgndspd' * fundamentals;
        testlms = testspd' * fundamentals;
        LMS = (testlms - bkgndlms) ./ bkgndlms;
        triplets = setfield(triplets, description, LMS);
        fprintf('%s => [%.3f, %.3f, %.3f] CC\n', char(description), LMS)
    end
end