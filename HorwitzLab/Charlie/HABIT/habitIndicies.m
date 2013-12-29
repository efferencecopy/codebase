% indicies for habituation paradigm

actGaborRInd = strmatch('Act_R_gun', [stro.sum.trialFields(1,:)]);
actGaborGInd = strmatch('Act_G_gun', [stro.sum.trialFields(1,:)]);
actGaborBInd = strmatch('Act_B_gun', [stro.sum.trialFields(1,:)]);
colorDirInd = strmatch('colordir', [stro.sum.trialFields(1,:)]);
correctInd = strmatch('trial_correct', [stro.sum.trialFields(1,:)]);
gaborFramesInd = strmatch('gabor_frames', [stro.sum.trialFields(1,:)]);
gaborLambdaInd = strmatch('gabor_lambda', [stro.sum.trialFields(1,:)]);
gaborRInd = strmatch('R_gun', [stro.sum.trialFields(1,:)]);
gaborGInd = strmatch('G_gun', [stro.sum.trialFields(1,:)]);
gaborBInd = strmatch('B_gun', [stro.sum.trialFields(1,:)]);
gaborOnInd = strmatch('gabor_on', [stro.sum.trialFields(1,:)]);
gaborXInd =  strmatch('gabor_pos_x', [stro.sum.trialFields(1,:)]);
gaborYInd =  strmatch('gabor_pos_y', [stro.sum.trialFields(1,:)]);
habitOnInd = strmatch('habit_on', [stro.sum.trialFields(1,:)]);
habitFramesInd = strmatch('habit_frames', [stro.sum.trialFields(1,:)]);
repHabitOnInd = strmatch('rep_habit_on', [stro.sum.trialFields(1,:)]);
repHabitOffInd = strmatch('habit_lastframe', [stro.sum.trialFields(1,:)]);
repGaborOnInd = strmatch('rep_gabor_on', [stro.sum.trialFields(1,:)]);
repGaborOffInd = strmatch('gabor_lastframe', [stro.sum.trialFields(1,:)]);

