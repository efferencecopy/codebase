%% rudamentary DTstim analysis

%open a file
% SEDNA
% {'S082811002','S082811003','S082811005','S082811006','S082811007','S082811008'}  %DTstim with laser
% {'S082811004','S082811009'} %DTstim w/o laser
% {'S082611005','S082611006'} %DTstim w/laser x3 colors

% KALI 11/27
% {'K112711004', 'K112711006', 'K112711013'}; %pref orientation with stim
% {'K112711005','K112711011','K112711012'}; %orthogonal orientation with stim
% {'K112711007', 'K112711008'}; % different axis of action
% {'K112711009'}; %delay fields with same stim protocol as DT
% {'K112711014'}; %delay fields 100ms stim protocol
% {'K112711016'}; %delay fields 200ms stim protocol

% KALI 11/29
%{'K112911003', 'K112911006', 'K112911009', 'K112911011', 'K112911012', 'K112911014'} % Null orientation
%{'K112911004', 'K112911007', 'K112911008', 'K112911010', 'K112911015', 'K112911016'} % Pref orientation


fin
filename = {'K112711004', 'K112711006', 'K112711013','K112711005','K112711011','K112711012'};
global DS;


if numel(filename)>1
    for a = 1:numel(filename)
        tmp{a} = dsobj(filename{a});
    end
    DS = strocat(tmp);
else
    DS = dsobj(filename{1});
end

%% Eye-pos and Rasters (Guided Saccades and Detection)

% clear out any potential junk in the workspace:
clear, clc, close all
global DS; %#ok<*REDEF>

% CONSTANTS %
PLOT = 'psth'; %'psth', 'raster'
binSize = 0.001; % in sec;
SYNC = 'gaborOn';   % other options: 'laserOn', 'laserOff', 'gaborOn', 'gaborOff', 'targon', 'targacq'
PRETIME = 0.200;    % in sec
POSTTIME = .800;   % in sec

%select a spike channel (if there's more than one)
spkChIdx = strncmpi('sig*', DS.sum.rasterCells,3);
if sum(spkChIdx)>1
    str = DS.sum.rasterCells(spkChIdx);
    idx = listdlg('PromptString', 'Select a spike channel', 'SelectionMode', 'single', 'ListString', str);
    spkCh = strcmpi(str{idx}, DS.sum.rasterCells);
else
    spkCh = find(DS.idx.spikes,1, 'first');
end

l_inRF = (DS.trial(:, DS.idx.gaborPosX) == DS.sum.exptParams.rfPosX) & (DS.trial(:, DS.idx.gaborPosY) == DS.sum.exptParams.rfPosY);
l_laserTrial = DS.trial(:, DS.idx.laserTrial) == 1;
counts.val = {};
counts.title = {};
figure
maxY = 0; %for psth
for a = 1:4;
    
    %pick the trials to analyze
    if a == 1;
        tList = l_inRF & l_laserTrial;
        titleString = sprintf('In RF + Laser');
    elseif a == 2;
        tList = ~l_inRF & l_laserTrial;
        titleString = sprintf('Out RF + Laser');
    elseif a == 3;
        tList = l_inRF & ~l_laserTrial;
        titleString = sprintf('In RF, No Laser');
    elseif a == 4;
        tList = ~l_inRF & ~l_laserTrial;
        titleString = sprintf('Out RF, No Laser');
    end
    
    %bail if there weren't any of a particular trial type
    if ~any(tList)
        continue
    end
    
    %find the synch time.
    switch lower(SYNC)
        case 'laseron';
            if sum(l_laserTrial&tList) == 0;
                continue
            end
            syncTimes = cellfun(@(x) x(1), DS.other(tList,1), 'uniformoutput', 0);
        case 'laseroff';
            if sum(l_laserTrial&tList) == 0;
                continue
            end
            syncTimes = cellfun(@(x) x(end), DS.other(tList,2), 'uniformoutput', 0);
        case 'gaboron'
            syncTimes = mat2cell(DS.trial(tList, DS.idx.repGaborOn), ones(sum(tList),1));
        case 'gaboroff'
            syncTimes = mat2cell(DS.trial(tList, DS.idx.gaborOff), ones(sum(tList),1));
        case 'targon'
            syncTimes = mat2cell(DS.trial(tList, DS.idx.targOn), ones(sum(tList),1));
        case 'targacq'
            syncTimes = mat2cell(DS.trial(tList, DS.idx.targAcq), ones(sum(tList),1));
    end
    
    %find the spike times
    spikeTimes = cellfun(@(x,y) (x-y), DS.ras(tList, spkCh), syncTimes, 'uniformoutput', 0);
    
    if sum(tList & l_laserTrial)>0;
        %make a time series representing the laser signal
        laserOnset = cellfun(@(x,y) x(1)-y, DS.other(tList,1), syncTimes, 'uniformoutput', 0);
        laserOffset = cellfun(@(x,y) x(end)-y, DS.other(tList,2), syncTimes, 'uniformoutput', 0);
    end
    
    %make a time series representing the targets
    tmp = DS.trial(tList, DS.idx.targOn)-cell2mat(syncTimes);
    targOn = mat2cell(tmp, ones(length(tmp),1));
    tmp = DS.trial(tList, DS.idx.targAcq)-cell2mat(syncTimes);
    targAcq = mat2cell(tmp, ones(length(tmp),1));
    
    %find the fp on time
    tmp = DS.trial(tList, DS.idx.fpOn) - cell2mat(syncTimes);
    fpon = mat2cell(tmp, ones(length(tmp),1));
    
    %Plot when the gabor comes on and off
    tmp = DS.trial(tList, DS.idx.repGaborOn)-cell2mat(syncTimes);
    gabOn = mat2cell(tmp, ones(length(tmp),1));
    tmp = DS.trial(tList, DS.idx.gaborOff)-cell2mat(syncTimes);
    gabOff = mat2cell(tmp, ones(length(tmp),1));
    
    %determine the number of spikes during the gabor interval:
    counts.val{a} = cellfun(@(x,y,z) sum((x>=y)&(x<z)), spikeTimes,gabOn,gabOff);
    counts.title{a} = titleString;
    
    %do the plotting
    switch lower(PLOT)
        case 'raster'
            subplot(2,2,a)
            hold on,
            counter = mat2cell([0:sum(tList)-1]', ones(sum(tList), 1), 1);
            cellfun(@(x, y) plot([x, x]', [zeros(1,length(x))+y; [ones(1, length(x)).*0.8 + y]], 'k'), spikeTimes, counter);
            cellfun(@(x,y) plot(x,y,'g*'), targOn, counter) % go sig (targs on)
            cellfun(@(x,y) plot(x,y,'c*'), targAcq, counter) % targs acq
            cellfun(@(x,y) plot(x,y,'k*'), fpon, counter) % fp onset
            cellfun(@(x,y) plot(x,y,'mx'), gabOn, counter) %gabor onset (approx)
            cellfun(@(x,y) plot(x,y,'mx'), gabOff, counter) %gabor offset (approx)
            if sum(tList & l_laserTrial)>0;
                cellfun(@(x,y) plot(x,y,'b*'), laserOnset, counter) %laser on
                cellfun(@(x,y) plot(x,y,'b*'), laserOffset, counter) %laser on
            end
            xlim([-PRETIME, POSTTIME])
            ylim([0, counter{end}+1])
            xlabel('Time (sec)')
            ylabel('Trial Number')
            title(titleString)
            hold off
            
        case 'psth'
            edges = repmat({-PRETIME:binSize:POSTTIME}, size(spikeTimes));
            psth = cellfun(@histc, spikeTimes, edges, 'uniformoutput', 0);
            psth = cellfun(@(x) x(:)', psth, 'uniformoutput', 0);
            psth = vertcat(psth{:}); %convert to mtx.
            psth = mean(psth,1)./binSize; %in sp/sec
            subplot(2,2,a)
            hold on,
            
            bar(edges{1}, psth, 'type', 'histc')
            xlim([-PRETIME, POSTTIME])
            xlabel('Time (sec)')
            ylabel('sp/sec')
            title(titleString)
            hold off
            maxY = max([maxY, get(gca, 'ylim')]);
            set(gca, 'linewidth', 2)
    end
end

switch lower(PLOT)
    case 'psth'
        for a = 1:4
            subplot(2,2,a)
            set(gca,'ylim', [0, maxY])
        end
end


%plot the average spike count during the gabor for the different laser
%conditions
avg = cellfun(@mean, counts.val);
sem = cellfun(@(x) std(x)./sqrt(length(x)), counts.val);
figure, hold on,
bar(1:4, avg, 'facecolor', 'w')
errorbar(1:4, avg, sem, 'k.', 'markersize', 0.1)
for a = 1:4
    text(a-.25, min(avg)*.6, counts.title{a})
end
set(gcf, 'position', [-1179         135         941         488]);
set(gca, 'xticklabel', [])
ylabel('spike count')


%% DOUBLE SIDED PSYCHOMETRIC FUNCTIONS (Detection)
% LINE OF NO RETURN!!!
% clear out any potential junk in the workspace:
clear, clc
global DS;

% determine which colors were used, and turn them into unit vectors
colors = reshape(DS.sum.exptParams.gaborColors, 3,3)';
blanks =sum(abs(colors),2) == 0;
colors(blanks,:) = [];
norms = sqrt(sum(colors.^2, 2));
colors = bsxfun(@rdivide, colors, norms);
nColors = size(colors,1);
norms = sqrt(sum(DS.LMS.^2,2));
unitVecs = bsxfun(@rdivide, DS.LMS, norms);

%determine the spatial frequencies.
trlSfs = 1./(DS.trial(:, DS.idx.gaborLambda) ./ DS.sum.exptParams.pixperdeg);  %in cpd
sptFreqs = unique(trlSfs);  %in cpd
nSfs = numel(sptFreqs);

% determine which trials were in the RF and also which were T1 choices
l_inRF = (sign(DS.trial(:, DS.idx.gaborPosX)) == sign(DS.sum.exptParams.rfPosX)) & (sign(DS.trial(:, DS.idx.gaborPosY)) == sign(DS.sum.exptParams.rfPosY));
corrects = DS.trial(:, DS.idx.correct);
l_T1choices = (l_inRF & corrects) | (~l_inRF & ~corrects);

%set up a vector of signed contrasts
signedCntrstLev = DS.trial(:, DS.idx.cntrstLev);
signedCntrstLev(~l_inRF) = signedCntrstLev(~l_inRF) * -1; %flip the sign of the contrast for the T2 location
signedCntrstLev(signedCntrstLev == -1) = 1; % +/- 1 both correspond to zero contrast
uniqueSignedCntrstLev = unique(signedCntrstLev);

%allocate space in the data matricies
dat.RTrepmat = repmat({nan(1,numel(uniqueSignedCntrstLev))}, [nColors, nSfs,2,2]);
dat.RTsem = repmat({nan(1,numel(uniqueSignedCntrstLev))}, [nColors, nSfs,2,2]);
dat.pT1choice = repmat({nan(1,numel(uniqueSignedCntrstLev))}, [nColors, nSfs,2]);
dat.pT1sbp = repmat({nan(1,numel(uniqueSignedCntrstLev))}, [nColors, nSfs,2]);
dat.nTrials = repmat({nan(1,numel(uniqueSignedCntrstLev))}, [nColors, nSfs,2]);
dat.norms = repmat({nan(1,numel(uniqueSignedCntrstLev))}, [nColors, nSfs,2]);
trialCounter = 0;
for clr = 1:nColors
    l_color = DS.trial(:, DS.idx.colordir) == clr;
    
    for sf = 1:nSfs
        l_sf = trlSfs == sptFreqs(sf);
        dat.sptFreq(1,sf) = sptFreqs(sf); %explictly make this a row vector.
        
        for cntrst = 1:numel(uniqueSignedCntrstLev);
            cntrstType = uniqueSignedCntrstLev(cntrst);
            l_contrast = signedCntrstLev == cntrstType;
            
            
            for lsr = 1:2;
                l_laser = DS.trial(:, DS.idx.laserTrial) == (lsr-1);
                
                %identify thet trials that share the common
                %characteristics. Treat the zero contrast condition
                %specially
                if cntrstType == 1;
                    tList = l_contrast & l_laser;
                else
                    tList = l_color & l_sf & l_contrast & l_laser;
                    %define the color direction here. Doing so earlier
                    %might requre using softEq (which I hate) but
                    %means that the code defines the color more often
                    %than it needs to.
                    if sum(tList) > 0
                        tmp = unique(unitVecs(tList,:), 'rows');
                        if size(tmp,1)>1
                            warning('more than one color per condition??')
                            disp(tmp)
                        end
                        dat.color(clr,:) = mean(tmp,1);
                    end
                end
                
                nTrials = sum(tList);
                if nTrials>0
                    nT1choices = sum(l_T1choices(tList));
                    pT1choice = nT1choices./nTrials;
                else
                    nT1choices = nan;
                    pT1choice = nan;
                end
                dat.pT1choice{clr,sf,lsr}(cntrst) = pT1choice;
                dat.pT1sbp{clr,sf,lsr}(cntrst) = std(binornd(nTrials, pT1choice, 1, 5000)./nTrials);
                dat.nTrials{clr,sf,lsr}(cntrst) = nTrials;
                tmp = unique(norms(l_color&l_contrast)).*sign(cntrstType);
                if numel(tmp)>1
                    warning('more than one norm for this condition!!')
                    disp(tmp)
                end
                dat.norms{clr,sf,lsr}(cntrst) = mean(tmp);
                
                % Calculate RTs for correct and incorrect trials separately
                for c = 1:2
                    RT_tList = tList & (corrects == (c-1));
                    RTs = DS.trial(RT_tList, DS.idx.targAcq) - DS.trial(RT_tList, DS.idx.targOn);
                    dat.RT{clr,sf,lsr,c}(cntrst)  = mean(RTs);
                    dat.RTsem{clr,sf,lsr,c}(cntrst) = std(RTs)./sqrt(numel(RTs));
                end
                
                %update the trial counter, but only for non-zero cntrst
                if cntrstType ~= 1;
                    trialCounter = trialCounter + sum(tList);
                end
                
            end %for laser
        end %for contrasts
    end %for sfs
end %for color

%compare the number of trials analyzed to the total number of trials.
trialCounter = trialCounter + sum(signedCntrstLev == 1);
fprintf('Total trials: %d, trials analyzed = %d\n', size(DS.trial,1), trialCounter);

%plot the 2sided psychometric functions
plotColors = {'k', 'r', 'b', 'g', 'm'};
for sf = 1:nSfs;
    for clr = 1:nColors;
        figure
        hold on,
        lsrStyles = {'-', '--'};
        for lsr = 1:2
            e = errorbar(dat.norms{clr, sf, lsr}, dat.pT1choice{clr, sf, lsr}, dat.pT1sbp{clr, sf, lsr}, [lsrStyles{lsr}, 'o'], 'color', plotColors{clr});
            set(e, 'linewidth', 2)
        end
        plot(0, 0.5, 'r+', 'markersize', 14)
        xx = horzcat(dat.norms{clr, sf, :});
        xlim([min(xx)*1.05, max(xx)*1.05])
        ylim([-.02, 1.02])
        ylabel('P(T1 Choice)');
        xlabel('Cone Contrast');
        legend('No Laser', 'Laser', 'location', 'southeast')
        legend boxoff
        title(sprintf('color: [%.3f, %.3f, %.3f], Spt Freq: %.3f',...
            dat.color(clr,1), dat.color(clr,2), dat.color(clr,3), dat.sptFreq(sf)));
        hold off
    end
end

%plot the 2sided chronometric functions
plotColors = {'k', 'b', 'g', 'm'};
for sf = 1:nSfs;
    for clr = 1:nColors;
        figure
        hold on,
        lsrStyles = {'-', ':'};
        for crt = 1:2
            for lsr = 1:2
                %p = errorbar(dat.norms{clr, sf, lsr}, dat.RT{clr, sf, lsr, crt}, dat.RTsem{clr, sf, lsr, crt}, [lsrStyles{lsr}, 'o'], 'color', plotColors{crt});
                p = plot(dat.norms{clr, sf, lsr}, dat.RT{clr, sf, lsr, crt}, [lsrStyles{lsr}, 'o'], 'color', plotColors{crt});
                set(p, 'linewidth', 2)
            end
        end %for correct
        xx = horzcat(dat.norms{clr, sf, :});
        xlim([min(xx)*1.05, max(xx)*1.05])
        ylabel('Reaction Times (sec)');
        xlabel('Cone Contrast');
        legend('Incorrect No Laser', 'Incorrect Laser', 'Correct No Laser', 'Correct Laser', 'location', 'southeast')
        legend boxoff
        title(sprintf('color: [%.3f, %.3f, %.3f], Spt Freq: %.3f',...
            dat.color(clr,1), dat.color(clr,2), dat.color(clr,3), dat.sptFreq(sf)));
        hold off
    end %for clr
end %for sf

%plot the difference in RTs b/w stim and no stim
for sf = 1:nSfs;
    for clr = 1:nColors;
        figure
        hold on,
        RTdiff = dat.RT{clr, sf, 1} - dat.RT{clr,sf,2};
        plot(dat.norms{clr,sf,1}, RTdiff*1000, '-b.')
        xx = horzcat(dat.norms{clr, sf, :});
        xlim([min(xx)*1.05, max(xx)*1.05])
        ylabel('Reaction Times (stim - no stim');
        xlabel('Cone Contrast');
        legend boxoff
        title(sprintf('color: [%.3f, %.3f, %.3f], Spt Freq: %.3f',...
            dat.color(clr,1), dat.color(clr,2), dat.color(clr,3), dat.sptFreq(sf)));
        hold off
    end
end
%% GUIDED SACCADE ANALYSIS (DELAY FIELDS)

% clear out any potential junk in the workspace:
clear, clc
global DS;


l_inRF = (DS.trial(:, DS.idx.gaborPosX) == DS.sum.exptParams.rfPosX) & (DS.trial(:, DS.idx.gaborPosY) == DS.sum.exptParams.rfPosY);
l_laserTrial = DS.trial(:, DS.idx.laserTrial) == 1;
for a = 1:4;
    
    %pick the trials to analyze
    if a == 1;
        tList = l_inRF & l_laserTrial;
        titleString = sprintf('In RF + Laser');
    elseif a == 2;
        tList = ~l_inRF & l_laserTrial;
        titleString = sprintf('Out RF + Laser');
    elseif a == 3;
        tList = l_inRF & ~l_laserTrial;
        titleString = sprintf('In RF, No Laser');
    elseif a == 4;
        tList = ~l_inRF & ~l_laserTrial;
        titleString = sprintf('Out RF, No Laser');
    end
    
    %bail if there weren't any of a particular trial type
    if ~any(tList)
        continue
    end
    
    
    % find the sacccade RTs
    RT.dat{a} = DS.trial(tList,DS.idx.targAcq) - DS.trial(tList,DS.idx.targOn);
    RT.title{a} = titleString;
end


% plot histos of RTs
figure
maxRT = max(cellfun(@max, RT.dat));
edges = linspace(0,maxRT,40);
for a = 1:4
    subplot(2,2,a), hold on,
    n=histc(RT.dat{a}, edges);
    bar(edges,n, 'type', 'histc')
    ylims = get(gca, 'ylim');
    plot(mean(RT.dat{a}), ylims(2)*0.98, 'mv', 'markerfacecolor', 'm')
    title(RT.title{a})
    xlim([0, maxRT])
    xlabel('time (sec)')
    ylabel('count')
end


% plot the means and SDs for inRF vs. laser trials
figure, hold on,
errorbar([1,2], [mean(RT.dat{strcmpi('In RF + Laser', RT.title)}), mean(RT.dat{strcmpi('Out RF + Laser', RT.title)})],...
                [stdErr(RT.dat{strcmpi('In RF + Laser', RT.title)}), stdErr(RT.dat{strcmpi('Out RF + Laser', RT.title)})], '-ko')
errorbar([1,2], [mean(RT.dat{strcmpi('In RF, No Laser', RT.title)}), mean(RT.dat{strcmpi('Out RF, No Laser', RT.title)})],...
                [stdErr(RT.dat{strcmpi('In RF, No Laser', RT.title)}), stdErr(RT.dat{strcmpi('Out RF, No Laser', RT.title)})],'-ro')
set(gca, 'xtick', [1,2], 'xticklabel', {'In RF', 'Out RF'})
ylabel('time (sec)')
legend('Laser', 'No Laser')

%run an anova to see if the values are significant
tmpDat = [];
tmpPos = [];
tmpLaser = [];
for a = 1:length(RT.dat)
    tmpDat = [tmpDat; RT.dat{a}(:)];
    
    if ~isempty(regexpi(RT.title{a}, 'No Laser'))
        tmpLaser = [tmpLaser; repmat('No Laser', numel(RT.dat{a}), 1)];
    else
        tmpLaser = [tmpLaser; repmat('__Laser_', numel(RT.dat{a}), 1)];
    end
    
    if ~isempty(regexpi(RT.title{a}, 'In RF'))
        tmpPos = [tmpPos; repmat('_In RF', numel(RT.dat{a}), 1)];
    else
        tmpPos = [tmpPos; repmat('Out RF', numel(RT.dat{a}), 1)];
    end
end

[P,T,STATS,] = anovan(tmpDat, {tmpPos, tmpLaser}, 'model', 'full');


%% SINGLE SIDED PSYCHOMETRIC FUNCTIONS 

% clear out any potential junk in the workspace:
clear, clc, close all
global DS;

norms = sqrt(sum(DS.LMS.^2,2));
units = bsxfun(@rdivide, DS.LMS, norms);
nColors = numel(unique(DS.trial(:, DS.idx.colordir)));
nSfs = numel(unique(DS.trial(:, DS.idx.gaborLambda)));
lambdas = unique(DS.trial(:, DS.idx.gaborLambda));
sptFreqs = 1./(lambdas ./ DS.sum.exptParams.pixperdeg);
nContrasts = numel(unique(DS.trial(:, DS.idx.cntrstLev)));


dat = [];
for clr = 1:nColors
    l_color = DS.trial(:, DS.idx.colordir) == clr;
    
    for sf = 1:nSfs
        l_sf = DS.trial(:, DS.idx.gaborLambda) == lambdas(sf);
        
        for lsr = 1:2
            l_laser = DS.trial(:, DS.idx.laserTrial) == (lsr-1);
            
            for cnt = 1:nContrasts
                l_contrast = DS.trial(:, DS.idx.cntrstLev) == cnt;
                
                if cnt == 1
                    tList = l_contrast & l_laser;
                else
                    tList = l_contrast & l_color & l_sf & l_laser;
                end
                
                if sum(tList)<= 0
                    continue
                end
                
                nTrials = sum(tList);
                nCorrect = sum(DS.trial(tList, DS.idx.correct));
                dat.pCorrect{sf, clr, lsr}(cnt) = nCorrect./nTrials;
                dat.norms{sf, clr, lsr}(cnt) = unique(norms(tList));
                dat.nTrials{sf, clr, lsr}(cnt) = nTrials;
            end
        end
    end
end


for clr = 1:nColors
    for sf = 1:nSfs
        figure
        hold on,
        semilogx(dat.norms{clr,sf,1}, dat.pCorrect{clr, sf, 1}, '-ko');
        if size(dat.pCorrect,3)>1 %i.e. there are laser trials present
            semilogx(dat.norms{clr, sf, 2}, dat.pCorrect{clr, sf, 2}, ':ko')
            legend('No Laser', 'Laser')
        end
    end
end






