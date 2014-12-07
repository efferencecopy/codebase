function makeStimulusFile(params)

% params should have
%
% params.type     =>  'train', 'pulse'
% params.si       =>  the sample INTERVAL (needs to be an iteger)
% params.tStart   =>  the time of the first pulse
% params.pAmp     =>  A vector of amplitudes for the pulse height
% params.pFreq    =>  A vector of frequencies for the pulse train
% params.nReps    =>  Number of repeates each stimulus should be presented



header{1,1} = 'ATF';
header{1,2} = 1;


fileID = fopen('celldata.atf','w');
formatSpec = '%s \t %d \t %2.1f \t %s\n';
[nrows,ncols] = size(C);
for row = 1:nrows
    fprintf(fileID,formatSpec,C{row,:});
end
fclose(fileID);