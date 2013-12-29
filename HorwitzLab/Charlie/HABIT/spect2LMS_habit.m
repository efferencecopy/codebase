function triplets = spect2LMS_habit

%unpacks a bunch of spectral data from previous measurments using the PR650
%and then outputs the cone contrast triplets for each of the spectra. I'm
%using a function b/c it allows for subfunctions...
%
% CAH 10/09

%find the appropriate path and cal structure
presDir = pwd;
cd('N:\NexFiles\Charlie\Paradigm Test Data\Habit')

%load in the fundamentals:
t = load ('T_cones_smj.mat');
fundamentals = t.T_cones_smj';

% data from the gabor
disp('    GABOR')
bk_gabor = load('bkgnd_gabor.mat');
lspd = load('l_gabor.mat');
lspd = lspd.l_gabor;
mspd = load('m_gabor.mat');
mspd = mspd.m_gabor;
sspd = load('s_gabor.mat');
sspd = sspd.s_gabor;

%print out the triplet values to the workspace:
triplets = [];
spd2lms(bk_gabor.bkgnd_gabor, lspd, 'L_iso');
spd2lms(bk_gabor.bkgnd_gabor, mspd, 'M_iso');
spd2lms(bk_gabor.bkgnd_gabor, sspd, 'S_iso');


% data from the habituator
disp('    HABITUATOR')
bk_habit = load('bkgnd_habit.mat');
lspd = load('l_habit.mat');
lspd = lspd.l_habit;
mspd = load('m_habit.mat');
mspd = mspd.m_habit;
sspd = load('s_habit.mat');
sspd = sspd.s_habit;

%print out the triplet values to the workspace:
triplets = [];
spd2lms(bk_habit.bkgnd_habit, lspd, 'L_iso');
spd2lms(bk_habit.bkgnd_habit, mspd, 'M_iso');
spd2lms(bk_habit.bkgnd_habit, sspd, 'S_iso');


%differences b/w gabor and habit bkgnd?
disp('    HABIT VS GABOR BKGND')
spd2lms(bk_habit.bkgnd_habit, bk_gabor.bkgnd_gabor, 'bkgnd');

%repeated measurements of lIso gabor
disp('    REPEATED L ISO GABORS')
l_reps = load('l_gabor_repeats.mat');
l_reps = l_reps.l_gabor_repeats;
spd2lms(bk_gabor.bkgnd_gabor, l_reps, 'l_reps')

%repeated measurements of lIso gabor
disp('    REPEATED S ISO GABORS')
s_reps = load('s_gabor_repeats.mat');
s_reps = s_reps.s_gabor_repeats;
spd2lms(bk_gabor.bkgnd_gabor, s_reps, 's_reps')

%be nice and cd back to the original dir
cd(presDir);


    function spd2lms(bkgndspd, testspd, description);
        bkgndspd = SplineSpd([380, 4, 101], bkgndspd(:), t.S_cones_smj);
        for a = 1:size(testspd,2)
            tmp_spd = SplineSpd([380, 4, 101], testspd(:,a), t.S_cones_smj);
            bkgndlms = bkgndspd' * fundamentals;
            testlms = tmp_spd' * fundamentals;
            LMS(a,:) = (testlms - bkgndlms) ./ bkgndlms;
        end
        if size(LMS,1) == 1
            triplets = setfield(triplets, description, LMS);
            fprintf('%s => [%.3f, %.3f, %.3f] CC\n', char(description), LMS)
        elseif size(LMS,1) > 1
            disp(LMS)
        end
    end
end