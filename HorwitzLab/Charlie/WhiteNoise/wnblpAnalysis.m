%   A BUNCH OF SCRIPTS TO LOOK AT WHITE NOISE DATA FROM THE ARRAY

%% DEFINE AND OPEN A DATAFILE
fin
global wn
fileName = 'A051111001';
wn = wnobj(fileName);

%% RAW LFP AND BLP SYNCHED TO REWARD ONSET

clear
global wn

%define useful trial events
nTrials = size(wn.trial,1);
oneVec = ones(nTrials,1);
rewOn = mat2cell(wn.trial
stimOn = mat2cell(wn.trial(:, wn.idx.stim_on), oneVec);
% nFrames = mat2cell(wn.trial(:, wn.idx.num_frames), oneVec);
% stimOff = cellfun(@(x,y,z) (x/y)+
% stimOff = mat2cell(wn.trial(:, wn.idx.stim_off, oneVec);


nLFPs = sum(wn.idx.lfp);
for a = 1:nLFPs;
    
end

%% BLP FOR DIFFERENT COLOR DIRECTIONS