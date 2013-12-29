function DTonline

%******** GLOBALS AND CONSTANTS ***********%

udpCom.sock = [];
udpCom.port = 6665;
udpCom.mac = '192.168.1.122';
udpCom.rexip = '192.168.1.120';
udpCom.plexip = '192.168.1.121';
udpCom.rexbuff = 8000;
maxNumNeurons = 3; %hard coding the number of neurons to analyze

[C, trialParamsHeader, exptParamsHeader] = DTCodes();
%*******************************************%

%Initialize some variables
dControls('init');


%SETTING UP THE EXTERNAL EQUIPMENT
[udpCom.sock, success] = pnetStart(udpCom.port); %open a connectionless udp socket
plxServer = PL_InitClient(0);  %open a connection with the plexon MAP box
dControls('newTrial'); %setup the trial parameters on the plexon server

%MAIN LOOP
allDone = 0;
while ~allDone
    allDone = checkESCkey(allDone);
    allDone = dealWithUDPMessages(udpCom.sock, allDone);
    
    %for debugging.
    if allDone ==2;
        keyboard
        allDone = 0;
    end
    
    
    %pull across ecodes and spike data
    [nStamps, eventList]  = PL_TrialEvents(plxServer, C.EOTCD, 1); %set a time out to prevent crashes in cases where the plx stream is not on
    
    %don't do any analysis till all the codes have come across
    if isempty(find(eventList == -C.EOTCD, 1))
        drawnow
        continue
    else
        [n, spikeList] = PL_TrialSpikes(plxServer, C.EOTCD, 1);
        d.tParams.codes = -eventList(:,2);
        d.tParams.cdtimes = eventList(:,1) .* (25*10^-6);
        d.tParams.spikeIds = spikeList(:,3);
        d.tParams.spktimes = spikeList(:,1) .* (25*10^-6);%one tick of plexon is 25usec
    end

    %if it's the first trial pull across the relavent header info.
    if find(d.tParams.codes >= 8000, 1)
        dControls('init'); %so that you can clear out an expt by clicking reset states
        setupExperiment();
        setupFig(maxNumNeurons);
    end

    %bail on this loop if the trial was bad, or if there's no header
    if isscalar(find(d.tParams.codes == C.ABORTCD, 1));
        dControls('newTrial')
        continue %bail for bad trials
    elseif ~d.eParams.gotHeader;
        dControls('newTrial')
        continue %bail if you don't have a header
    end
    computeTrialStats(maxNumNeurons);
    updatePlot()    
    dControls('newTrial');
end %main loop

%close the connections with external equipment before exiting
PL_Close(plxServer);

    
    %************************%
    %  NESTED SUBFUNCTIONS   %
    %************************%
    
    
% %
% % Pulls accros the header and unpacks the calibration information
% %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function setupExperiment();
        %roll over the exptParamsHeader and bring across all the content
        for a = 1:size(exptParamsHeader, 2)
            num = dat2num(d.tParams.codes, exptParamsHeader{4,a}, C.VALOFFSET, exptParamsHeader{2,a}, exptParamsHeader{3,a});
            d.eParams = setfield(d.eParams, exptParamsHeader{1,a}, num{1});
        end
        
        %deal with the color directions presented. convert to unit vectors
        colorDirs = reshape(d.eParams.RF_colors, 3,3)'; %colors accross rows
        noColor = find(sum(abs(colorDirs), 2) == 0);
        if ~isempty(noColor)
            colorDirs(noColor, :) = [];
        end
        norms = sqrt(diag(colorDirs * colorDirs'));
        
        %pull out the rest of the important experimental characteristics
        x = 0:255; %the normal range of the gamma look up table
        xx = linspace(0, 255, 2^16); %the desired quantization of the gammaTable
        g1 = reshape(d.eParams.gamma_table, 256, 3);
        d.eParams.gamma_table = [spline(x, g1(:,1), xx)', spline(x, g1(:,2), xx)', spline(x, g1(:,3), xx)'];
        d.eParams.m_mtx = reshape(d.eParams.m_mtx, 3, 3);
        d.eParams.bkgndlms = d.eParams.m_mtx*[d.eParams.bkgnd_r; d.eParams.bkgnd_g; d.eParams.bkgnd_b];
        d.eParams.colorDirs = colorDirs ./ repmat(norms, 1, 3);
        d.eParams.nContrasts = d.eParams.nContrasts+1; %the slave adds the 'zero' condition
        
        %update the flag to signify that we've received the header
        d.eParams.gotHeader = 1;
    end



% %
% % This function fills up the d.stats array. For behavioral data, the
% % relevant arrays are:
% %   d.stat.norms
% %   d.stats.nTrials
% %   d.stats.nTrialsCorrect
% % These arrays are of the demensions {nColors x nSfs}(1:nContrasts). The
% % relevant neural data arrays are:
% %   d.stats.rocArea
% %   d.stats.spikeStat
% % These arrays are of the demnsions {nColors x nSfs x nUnits}(1:nContrasts)
% % and contain the area under the ROC. spikeStat can be either the counts or
% % the f(1) response amplitude.
% %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function computeTrialStats(maxNumNeurons)
       %convert the ecodes
        sptPeriod = dat2num(d.tParams.codes, C.GABORLAMBDACD, C.VALOFFSET, 'float', 0);
        tSfs = (1./sptPeriod{1}) .* d.eParams.pixperdeg;
        cntrst = dat2num(d.tParams.codes, C.CNTRSTLEVCD, C.VALOFFSET, 'int', 0);
        tCntrstLev = cntrst{1};
        tCorrect = dat2num(d.tParams.codes, C.CORRECTCD, C.VALOFFSET, 'int', 0);
        
        %deal with the neural data b/4 the behavioral data. At this point I
        %don't know what color/sfs condition this trial belongs to, so I'll
        %hang onto the rate data until I figure out where to put it. I'm
        %doing this so that I know how many arrays to initialize (which
        %happens after I determine each trial's sf.
        flashOnidx = find(d.tParams.codes == C.STIMONREPORTCD, 1);
        flashOnTime = d.tParams.cdtimes(flashOnidx);
        flashOffidx = find(d.tParams.codes == C.STIMOFFCD, 1);
        flashOffTime = d.tParams.cdtimes(flashOffidx);
        trlTime = flashOffTime-flashOnTime;
        nFrames = dat2num(d.tParams.codes, C.NUMFRAMESCD, C.VALOFFSET, 'long', 0);
        timeFromFrameRate = nFrames{1} ./ d.eParams.frame_rate;
        
        trialNeurons = unique(d.tParams.spikeIds);
        trialNeurons(trialNeurons == 0) = []; %don't count unsorted spikes
        trialNeurons(trialNeurons > maxNumNeurons) = []; % don't consider more than maxNum Neurons
        nNeurons = length(trialNeurons);
        for a = 1:nNeurons;
            l_tSpikes = d.tParams.spikeIds == trialNeurons(a);
            tSpikes = d.tParams.spktimes(l_tSpikes);
            tSpikes(tSpikes <= flashOnTime) = [];
            tSpikes(tSpikes > flashOffTime) = [];
            tRate(a) = length(tSpikes) ./ timeFromFrameRate;
        end
        
        
        %deal with the spatial frequency. I'm dynamically allocating space
        %in arrays, so do some initialization if it's the first trial of a
        %particular sf.
        if ~isempty(find(d.stats.sfs == tSfs,1))
            tSfsIdx = find(d.stats.sfs == tSfs);
        else
            nColors = size(d.eParams.colorDirs, 1);
            d.stats.sfs(end+1) = tSfs;
            tSfsIdx = find(d.stats.sfs == tSfs, 1);
            [d.stats.norms{1:nColors, tSfsIdx}] = deal([]);
            [d.stats.nTrials{1:nColors, tSfsIdx}] = deal(zeros(1, d.eParams.nContrasts));
            [d.stats.nTrialsCorrect{1:nColors, tSfsIdx}] = deal(zeros(1, d.eParams.nContrasts));
            [d.stats.meanRates{1:nColors, tSfsIdx, 1:maxNumNeurons}] = deal(zeros(1, d.eParams.nContrasts));
            [d.stats.unitNumTrials{1:nColors, tSfsIdx, 1:maxNumNeurons}] = deal(zeros(1, d.eParams.nContrasts));
        end
        
        %now deal with the trial statistics
        if tCntrstLev == 1; %'zero' contrast is special
            for a = 1:size(d.eParams.colorDirs, 1);
                for b = 1:length(d.stats.sfs);
                    for c = 1:nNeurons
                        unitIdx = trialNeurons(c);
                        oldMeanRate = d.stats.meanRates{a,b,unitIdx}(tCntrstLev);
                        oldNTrials = d.stats.unitNumTrials{a,b,unitIdx}(tCntrstLev);
                        newMeanRate = ((oldMeanRate .* oldNTrials) + tRate(c)) ./ (oldNTrials+1);
                        d.stats.meanRates{a,b,c}(tCntrstLev) = newMeanRate;
                        d.stats.unitNumTrials{a,b,unitIdx}(tCntrstLev) = oldNTrials+1;
                    end
                    d.stats.norms{a,b}(tCntrstLev) = 10^-5; %something close to zero
                    d.stats.nTrialsCorrect{a,b}(tCntrstLev) = d.stats.nTrialsCorrect{a,b}(tCntrstLev) + tCorrect{1};
                    d.stats.nTrials{a,b}(tCntrstLev) = d.stats.nTrials{a,b}(tCntrstLev) + 1;
                end
            end
        else
            %find the color direction
            RGB = dat2num(d.tParams.codes, [C.ACTRGUNCD C.ACTGGUNCD C.ACTBGUNCD], C.VALOFFSET, 'long', 0);
            rgb = [d.eParams.gamma_table(RGB{1},1); d.eParams.gamma_table(RGB{2}, 2); d.eParams.gamma_table(RGB{3}, 3)];
            lms = d.eParams.m_mtx * rgb;
            LMS = (lms-d.eParams.bkgndlms) ./ d.eParams.bkgndlms;
            tColorDir = LMS./norm(LMS);
            [maxProj, tColorType] = max(abs(d.eParams.colorDirs * tColorDir));
            
            %fill up the trial stats arrays
            d.stats.norms{tColorType, tSfsIdx}(tCntrstLev) = norm(LMS);
            d.stats.nTrials{tColorType, tSfsIdx}(tCntrstLev) = d.stats.nTrials{tColorType, tSfsIdx}(tCntrstLev) +1;
            d.stats.nTrialsCorrect{tColorType, tSfsIdx}(tCntrstLev) = d.stats.nTrialsCorrect{tColorType, tSfsIdx}(tCntrstLev) + tCorrect{1};
        
            %Update the neural data only when the stimulus is in the RF
            xPos = dat2num(d.tParams.codes, C.FLASHPOSXCD, C.VALOFFSET, 'float', 0);
            yPos = dat2num(d.tParams.codes, C.FLASHPOSYCD, C.VALOFFSET, 'float', 0);
            if (d.eParams.rf_x == xPos{1}) && (d.eParams.rf_y == yPos{1})
                for c = 1:nNeurons
                    unitIdx = trialNeurons(c);
                    oldMeanRate = d.stats.meanRates{tColorType, tSfsIdx, unitIdx}(tCntrstLev);
                    oldNTrials = d.stats.unitNumTrials{tColorType, tSfsIdx, unitIdx}(tCntrstLev);
                    newMeanRate = ((oldMeanRate .* oldNTrials) + tRate(c)) ./ (oldNTrials+1);
                    d.stats.meanRates{tColorType, tSfsIdx, unitIdx}(tCntrstLev) = newMeanRate;
                    d.stats.unitNumTrials{tColorType, tSfsIdx, unitIdx}(tCntrstLev) = oldNTrials+1;
                end
            end
        end

        %refresh the contents of the ui controls
        sfsNames = num2str(sort(d.stats.sfs)', 3);
        set(d.figs.ui.sfsMenu, 'string', sfsNames);
        minTrialsPerCond = min(horzcat(d.stats.nTrials{:}));
        maxTrialsPerCond = max(horzcat(d.stats.nTrials{:}));
        set(d.figs.ui.nTrialsText, 'string', sprintf('Trials/Cond \n min: %d \n max: %d', minTrialsPerCond, maxTrialsPerCond));
        
        %timing violations?
        if abs(timeFromFrameRate - trlTime) > (1./d.eParams.frame_rate);
            disp('Timing violation!!');
            disp(timeFromFrameRate - trlTime)
        end
    end



% %
% % Simply sets up the figure according to the number of color directions
% %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function setupFig(maxNumNeurons)
        nColorDirs = size(d.eParams.colorDirs, 1);
        axHeight = 300;
        axWidth = 300;
        axGutter = 50;
        uiControlGutter = 170;
        rSideSlop = 10;
        figHeight = axHeight+100;
        figWidth = ((axWidth+axGutter)*nColorDirs)+uiControlGutter-rSideSlop;
        figLeft = 11+(axWidth.*(3-nColorDirs));
        
        %open a figure for the first run. clicking reset states on rex
        %clears out the current contents of the figure.
        d.figs.figHand = figure(1);
        clf(d.figs.figHand);
        set(d.figs.figHand, 'defaultaxesunits', 'pixels', 'position', [figLeft 220 figWidth figHeight]);

        
        %set up the plotting axes
        for a = 1:nColorDirs;
            ax = axes('position', [(uiControlGutter+((a-1).*(axWidth+axGutter))), 45, axWidth, axHeight]);
            set(ax, 'buttondownfcn', @getCntrstRange, 'UserData', d.eParams.colorDirs(a,:));
            hold on,
            [ax, d.figs.monkHands(a), d.figs.cellHands(a)] = plotyy([5], [0.3], [5], [25], 'semilogx', 'semilogx');
            d.figs.psyFunHands(a) = semilogx(ax(1), 5, 0.5,'go');
            set(d.figs.psyFunHands(a), 'buttondownfcn', @getCntrstRange);
            set(d.figs.monkHands(a), 'color', 'k', 'marker', 'o', 'buttondownfcn', @getCntrstRange);
            set(d.figs.cellHands(a), 'color', 'b', 'marker', 'o', 'buttondownfcn', @getCntrstRange);
            set(ax(2), 'ylim', [0 30], 'Xtick', [], 'YTickMode', 'auto', 'color', 'w', 'ycolor', 'k', 'box', 'off');
            set(ax(1), 'ylim', [0 1], 'Xscale', 'log', 'YTickMode', 'auto', 'color', 'none', 'ycolor', 'k', 'box', 'off');
            hold off
            title(sprintf('Color: [%.2f %.2f %.2f] \n psy alpha: n/a',...
                d.eParams.colorDirs(a,1), d.eParams.colorDirs(a,2), d.eParams.colorDirs(a,3)));
            xlabel('Cone Contrast');
            axes(ax(1)) %makes both y axes visible (???)
            if a==1;
                set(get(ax(1),'Ylabel'),'String','Probability Correct');
                if nColorDirs > 1;
                    set(ax(2), 'yTickLabel', [])
                end
            end
            if (a>1) && (a<nColorDirs)
                set(ax, 'yTickLabel', []); %sets both simultaneously
            end
            if a == nColorDirs
                set(get(ax(2), 'Ylabel'), 'String', 'Rate (imp/sec)');
                if nColorDirs > 1;
                    set(ax(1), 'yTickLabel', []);
                end
            end
        end
        
        %set up the uicontrols
        d.figs.ui.sfsMenu = uicontrol('style', 'popupmenu', 'string', 'none', 'callback', @updatePlot, 'position', [20, figHeight.*0.68, 70, 35]);
        d.figs.ui.sfsMenuText = uicontrol('style', 'text', 'string', 'Spt. Freq', 'position', [18, figHeight.*0.77, 78, 15], 'backgroundColor', get(gcf, 'color'));
        d.figs.ui.unitMenu = uicontrol('style', 'popupmenu', 'string', num2str([1:maxNumNeurons]'), 'callback', @updatePlot, 'position', [20, figHeight.*0.815, 70, 35]);
        d.figs.ui.unitMenuText = uicontrol('style', 'text', 'string', 'Unit No.', 'position', [18, figHeight.*0.91, 78, 15], 'backgroundColor', get(gcf, 'color'));
        d.figs.ui.nTrialsText = uicontrol('style', 'text', 'string', sprintf('Trials/Cond \n none'), 'position', [18, figHeight.*0.53, 78, 42], 'backgroundColor', get(gcf, 'color'));
        d.figs.ui.cntrstRange = uicontrol('style', 'text', 'string', sprintf('   Cntrst Range \n\n L: %.4f \n M: %.4f \n S: %.4f \n\n Scale: %.4f', nan, nan, nan, nan), 'position', [17, figHeight.*0.14, 88, 98], 'horizontalAlignment', 'left');
        d.figs.ui.rangeButton = uicontrol('style', 'pushbutton', 'string' , 'Compute', 'callback', @getCntrstRange, 'position', [29, figHeight.*0.05, 60, 25]);
        drawnow
    end
    
% %
% %
% % This function updates the plot window. It's called with varable input
% % arguments so that the ui control call backs can utilize this
% % function. I don't typically do anything with these inputs b/c by
% % default I query all the UI controls at the beginning of this
% % function. A given plot will only be presented or updated if there are
% % at least 1 trial per contrast for ALL color directions.
% %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function updatePlot(varargin)
        
        sfsUIVal = get(d.figs.ui.sfsMenu, 'value');
        sortSfs = sort(d.stats.sfs);
        sfsToPlot = sortSfs(sfsUIVal);
        sfsIdx = find(d.stats.sfs == sfsToPlot, 1);
        unitIdx = get(d.figs.ui.unitMenu, 'value');
        
        %update all the color directions. every contrast for a givin sf
        %must have more than one trial inorder to plot that sf.
        if all(all(vertcat(d.stats.nTrials{1:end, sfsIdx})));
            for a = 1:size(d.eParams.colorDirs, 1);
                hold on,
                ax = get(d.figs.monkHands(a), 'parent');
                axes(ax);

                %plot the raw data
                performance = d.stats.nTrialsCorrect{a, sfsIdx} ./ d.stats.nTrials{a, sfsIdx};
                set(d.figs.monkHands(a), 'Ydata', performance);
                set(d.figs.monkHands(a), 'Xdata', d.stats.norms{a, sfsIdx});
                set(get(d.figs.monkHands(a), 'parent'), 'ylim', [0.4 1.1])
                low = d.stats.norms{a,sfsIdx}(2).*0.8;
                high = d.stats.norms{a, sfsIdx}(end)*1.2;
                set(get(d.figs.monkHands(a), 'parent'), 'xlim', [low high]);
                set(d.figs.monkHands(a), 'Marker', 'o', 'color', 'k', 'linestyle', 'none')

                %plot the fit to the behavior
                % GDLH added 3/20/09
                err = (performance-.82).^2;
                guesses = [d.stats.norms{a, sfsIdx}(find(err == min(err),1)) 1 1];
                [alpha, beta, gamma] = weibullFit(d.stats.norms{a, sfsIdx}, performance, 'sse', guesses);
                % End of added code
                
                dx = (d.stats.norms{a, sfsIdx}(end) - d.stats.norms{a, sfsIdx}(1)) ./ 100;
                x = [d.stats.norms{a, sfsIdx}(1) : dx : d.stats.norms{a, sfsIdx}(end).*1.1];
                model = gamma - (gamma - 0.5).*exp(-(x./alpha).^beta);
                set(d.figs.psyFunHands(a), 'Ydata', model)
                set(d.figs.psyFunHands(a), 'Xdata', x);
                set(d.figs.psyFunHands(a), 'Marker', 'none', 'linestyle', '-', 'color', 'k')
                
                %plot the contrast response function
                if any(horzcat(d.stats.meanRates{1:size(d.eParams.colorDirs, 1), sfsIdx, unitIdx}));
                    set(d.figs.cellHands(a), 'Ydata', d.stats.meanRates{a, sfsIdx, unitIdx}, 'Xdata', d.stats.norms{a, sfsIdx})
                    set(d.figs.cellHands(a), 'color', 'b', 'linewidth', 1, 'marker', '.')
                    maxRate = max(horzcat(d.stats.meanRates{:,:,:}));
                    set(get(d.figs.cellHands(a), 'parent'), 'YLim', [0, maxRate.*1.1]);
                    set(get(d.figs.cellHands(a), 'parent'), 'xlim', [low high]);
                else
                    %so switching to an 'unavailable unit' blanks the old crf
                    set(d.figs.cellHands(a), 'color', 'w')
                end
                

                %update the title
                title(sprintf('Color: [%.2f %.2f %.2f], sfs: %.2f, unit: %d\n psyFun: %.2f',...
                    d.eParams.colorDirs(a,1), d.eParams.colorDirs(a,2), d.eParams.colorDirs(a,3), sfsToPlot, unitIdx, alpha));
                hold off
                drawnow
            end
        end
        drawnow
    end


% %
% % This function allows the user to click on an axis (in two places) and
% % define a contrast range. The text menue to on the right side of the
% % figure the displays this range in units suitable to use in the Rex user
% % menues
% %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function getCntrstRange(h, ev)
        [desiredContrasts, dummy] = ginput(2);
        colorDir = get(gca, 'userdata'); %colorDir stuffed into user data during setupFig
        LMS = colorDir .* max(desiredContrasts);
        LMS = LMS .* 100; % convert to %b/w 0&100 (which is what the user menu on Rex uses)
        scale = min(desiredContrasts) ./ max(desiredContrasts);
        set(d.figs.ui.cntrstRange, 'string', sprintf('   Cntrst Range \n\n L: %.4f \n M: %.4f \n S: %.4f \n\n Scale: %.4f', LMS(1), LMS(2), LMS(3), scale));
        drawnow
    end



% %
% % Gets called to manipulate the 'd' struct at the beginning of trials and
% % epxerimental runs.
% %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function dControls(process)
        switch process
            case 'newTrial'
                PL_TrialDefine(plxServer, -1000, -C.EOTCD, -1000, -C.EOTCD, 0, 0, [1], 0, 0);
                d.tParams.codes = [];
                d.tParams.cdtimes = [];
                d.tParams.spikeIds = [];
                d.tParams.spktimes = [];
                tCntrstLev = [];
                drawnow
            case 'init' %for the initilaization of the data struct
                d.eParams = struct('gotHeader', [0]);
                d.figs = struct();
                d.stats = struct();
                d.stats.sfs = [];
                d.stats.nTrials = {};
                d.stats.nTrialsCorrect = {};
                d.stats.norms = {};
                d.stats.meanRates = {};
                d.stats.unitNumTrials = {};
        end
    end
    
    
    
end %DTonline



%******************************%
%   NON-NESTED SUBFUNCTIONS    %
%******************************%
function allDone = checkESCkey(allDone)
    [keyIsDown, secs, keyCode] = KbCheck;
    if (keyCode(27)) % esc is key 27
        allDone = 1;
    end
    %for debugging
    if (keyCode (109)) %the minus key on the number pad
        allDone = 2;
    end
end

function stopNow = dealWithUDPMessages(socket, allDone);
    stopNow = allDone;
    
    msgSize = pnet(socket, 'readpacket', 250, 'noblock');
    if ~msgSize
        return
    end
    
    try
        msg = pnet(socket, 'read', msgSize, 'char');
        if strcmpi(msg(1:6), 'return')
            stopNow = 1;
            return
        else
            eval(msg)
        end
    catch
        fprintf('Unknown message %s\n', msg);
        wtf
    end
end


