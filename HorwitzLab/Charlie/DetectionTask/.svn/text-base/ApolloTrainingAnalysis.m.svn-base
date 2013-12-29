%% Analysis script to look at Apollo's behavioral data

fin % clear out the junk from the matlab command window

%specify a file name. The .nex file must be located in a place matlab can
%find it... in this case, in the NexFiles directory
filename = 'A121211001.nex';
DT = dtobj(filename); %open the file

%% (1)  Print out the psychometric functions

%run a small analysis routine. It should print out psychometric
%functions... "cellData" will likely be empty if no neruons are recorded.
[monkData, cellData, exptParams] = DTunpack(DT,1);

%psychometric thresholds are stored in monkData.alpha. Each row is a
%different color, each column is a different spatial frequency:
disp('Psychometric Thresholds: ')
monkData.alpha

disp('Color Directions: ')
exptParams.standColors %printing to command window

disp('Spatial Frequencies: ')
exptParams.sfs' %printing to command window


%% (2) Analysis to look for bias in choices 

clc

% define some important binary vectors
l_correct = DT.trial(:, DT.idx.correct);
l_inRF = DT.trial(:, DT.idx.flashX) == DT.sum.exptParams.rf_x;
l_T1choice = (l_inRF & l_correct) | (~l_inRF & ~l_correct);

% define some other things that will be helpful
sfs = 1./(DT.trial(:, DT.idx.gaborLambda) ./ DT.sum.exptParams.pixperdeg);
signedCntrstLev = DT.trial(:, DT.idx.cntrstLev);
signedCntrstLev(~l_inRF) = signedCntrstLev(~l_inRF) * -1; %flip the sign of the contrast for the T2 location
signedCntrstLev(signedCntrstLev == -1) = 1; % +/- 1 both correspond to zero contrast

% access the data and loop through each color/spatial freq/contrast
% condition. Store the p(T1 choice) data in a cell array with dimensions
% (Num Colors) x (Num Spatial Freqs)
nColors = size(exptParams.standColors,1);
nSfs = numel(exptParams.sfs);
pT1 = {};
signedContrasts = {};
colors = [];
sptFreqs = [];
trialCounter = 0;

for clr = 1:nColors
    l_clr = DT.trial(:, DT.idx.colorDir) == clr;
    
    for sf = 1:nSfs
        l_sf = softEq(exptParams.sfs(sf), sfs, 3, 'rows');
        sptFreqs(1,sf) = unique(sfs(l_sf));
        signedContrasts{clr, sf} = [fliplr(-exptParams.norms{clr}(2:end)), 0, exptParams.norms{clr}(2:end)];
        uniqueSignedContrasts = unique(signedCntrstLev);
        
        pT1{clr, sf} = nan(1, numel(uniqueSignedContrasts));
        for cnt = 1:numel(uniqueSignedContrasts);
            l_cnt = signedCntrstLev == uniqueSignedContrasts(cnt);
            
            if uniqueSignedContrasts(cnt) == 1
                tList = l_cnt;
            else
                tList = l_clr & l_sf & l_cnt;
                if sum(tList) > 0;
                    norms = sqrt(sum(DT.LMS(tList,:).^2, 2));
                    units = DT.LMS(tList,:) ./ repmat(norms, 1,3);
                    colors(clr,:) = unique(units, 'rows');
                end
            end
            
            if uniqueSignedContrasts(cnt) ~= 1
                trialCounter = trialCounter + sum(tList); %making sure each trial only gets counted once...
            end
            
            if sum(tList)>0
                pT1{clr, sf}(cnt) = sum(l_T1choice&tList) ./ sum(tList);
            end
        end
    end
end
trialCounter = trialCounter + (sum(signedCntrstLev == 1)); %add the zero contrast trials only once...
fprintf('total trials: %d, trials analyzed: %d\n', size(DT.trial,1),trialCounter) %making sure that each trial gets included once and only once.



%plot the results
for clr = 1:nColors
    for sf = 1:nSfs
        figure, hold on,
        plot(signedContrasts{clr, sf}, pT1{clr, sf}, 'o-')
        plot(0, 0.5, 'r+', 'markersize', 14)
        xlabel('Signed Contrast')
        ylabel('p(T1 choice)')
        title(sprintf('color: [%.3f, %.3f, %.3f], Spt Freq: %.3f',...
            colors(clr,1), colors(clr,2), colors(clr,3), sptFreqs(sf)));
    end
end


