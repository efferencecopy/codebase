function DTgui

% an interactive tool for doing data analysis on detection data


%open a control pannel. This will specify the different files to analyze
%and keep track of a file history. The control pannel calls a seperate
%function that sets up all the plots.
cntrlPanHand = currentFigOpen('DT GUI Control Pannel');
if isempty(cntrlPanHand)
    cPan.figHand = figure;
else
    cPan.figHand = cntrlPanHand;
end
set(cPan.figHand, 'name', 'DT GUI Control Pannel', 'numbertitle', 'off');
set(cPan.figHand, 'units', 'normalize', 'position', [0.1891    0.7754    0.6047    0.1309]);
uicontrol('style', 'text', 'string', 'Detection Files', 'units', 'normalize', 'position', [.4, .61, .2, .16], 'fontsize', 13, 'backgroundColor', get(gcf, 'color'));
cPan.dtTextHand = uicontrol('style', 'edit', 'string', 'input DT file', 'units', 'normalized', 'position', [.38, .42, .24, .18], 'backgroundColor', [1 1 1], 'fontsize', 12);
uicontrol('style', 'text', 'string', 'Grating Files', 'units', 'normalize', 'position', [.1, .76, .2, .16], 'fontsize', 13, 'backgroundColor', get(gcf, 'color'));
cPan.gtTextHand = uicontrol('style', 'edit', 'string', 'input GT file', 'units', 'normalized', 'max', 4, 'min', 1, 'position', [.08, .25, .24, .5], 'backgroundColor', [1 1 1], 'fontsize', 12);
cPan.applyHand = uicontrol('style', 'pushbutton', 'string', 'Apply', 'callback', {@analyzeData, cPan.figHand}, 'units', 'normalized', 'position', [.7, .63 .2, .2], 'fontsize', 11, 'interruptible', 'off', 'busyaction', 'cancel');
cPan.resetHand = uicontrol('style', 'pushbutton', 'string', 'Reset', 'callback', @resetTextFields, 'units', 'normalized', 'position', [0.7 .3 .2 .2], 'fontsize', 11, 'interruptible', 'off', 'busyaction', 'cancel');
set(cPan.figHand, 'userdata', cPan); 

    %
    %   RESETTING THE TEXT FIELDS
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function resetTextFields(h, ev)
        set(cPan.dtTextHand, 'string', 'input DT file')
        set(cPan.gtTextHand, 'string', 'input GT file')
    end
    
end %DTgui

%***********************************************************
%
%           NON NESTED SUBFUNCTIONS
%
%***********************************************************



%
%   ANALYZE DATA (ROUTER FUNCTION)
%
%   the main analysis script is a non-nested function and calls other
%   non-nested functions to carry out each analysis. GT and DT data are
%   stored in the userdata field of each figure. The grating figure has a
%   sturcture named 'g' that contains these data. The DT structure is
%   called 'd'.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function analyzeData(h, ev, cntrlPanHand)
    cPan = get(cntrlPanHand, 'userdata');
    set(cPan.applyHand, 'string', 'Finding Data Files...', 'FontSize', 9,'BackgroundColor',[0.92549 0.913725 0.847059]);
    drawnow;
    
    %open the appropriate files. Be nice and change the user input to all
    %caps in case they forgot.
    DTfname = upper(get(cPan.dtTextHand, 'string'));
    GTfname = upper(get(cPan.gtTextHand, 'string'));
    
    
    %setup the gratings GUI
    if ~strcmpi(GTfname, 'input GT file') && ~isempty(GTfname)
        %iterate over all the available GT files
        for a = 1:size(GTfname,1)
            fPath = findfile(GTfname(a,:));
            if isempty(fPath)
                set(cPan.applyHand, 'string', sprintf('GT file <%d> not found', a), 'FontSize', 11);
                set(cPan.applyHand,'BackgroundColor',[1 0.8 0.8]);
                drawnow
                return
            end
            g.GT{a} = nex2stro(findfile(GTfname(a,:)));
            g.GT{a}.orients = g.GT{a}.trial(:,find(strcmp(g.GT{a}.sum.trialFields(1,:),'orient')));
            g.GT{a}.sfs = g.GT{a}.trial(:,find(strcmp(g.GT{a}.sum.trialFields(1,:),'sf')));
            g.GT{a}.diams = g.GT{a}.trial(:,find(strcmp(g.GT{a}.sum.trialFields(1,:),'diam')));
            g.GT{a}.protocols = g.GT{a}.trial(:,find(strcmp(g.GT{a}.sum.trialFields(1,:),'protocol')));
            framerate = g.GT{a}.sum.exptParams.framerate; 
            nframes = g.GT{a}.trial(:,find(strcmp(g.GT{a}.sum.trialFields(1,:),'nframes')));
            g.GT{a}.stimon_t = g.GT{a}.trial(:,find(strcmp(g.GT{a}.sum.trialFields(1,:),'stim_on')));
            g.GT{a}.stimoff_t= g.GT{a}.trial(:,strcmp(g.GT{a}.sum.trialFields(1,:),'stim_off'));;
            spikeidx = find(strcmp(g.GT{a}.sum.rasterCells(1,:),getSpikenum(g.GT{a}))); 
            
            Lcc = g.GT{a}.trial(:,find(strcmp(g.GT{a}.sum.trialFields(1,:),'lcont')));
            Mcc = g.GT{a}.trial(:,find(strcmp(g.GT{a}.sum.trialFields(1,:),'mcont')));
            Scc = g.GT{a}.trial(:,find(strcmp(g.GT{a}.sum.trialFields(1,:),'scont')));
            
            %standardize the colors contained in the .trial field. The L
            %cones should have a positive sign. S-iso should also be
            %positive.
            LMS = [Lcc, Mcc, Scc];
            l_negL = LMS(:,1) < 0;
            LMS(l_negL,:) = LMS(l_negL,:) .* -1; %deal with all non-siso colors
            l_negSiso = ismember(sign(LMS), [0 0 -1], 'rows');
            LMS(l_negSiso,:) = LMS(l_negSiso,:) .* -1; %deal with the Siso colors
            g.GT{a}.colordirections = LMS;
            
            g.GT{a}.spikerates = [];
            g.GT{a}.baselines = [];
            baseline_t = 0.25; 
            for i = 1:size(g.GT{a}.trial,1)
                spiketimes = g.GT{a}.ras{i,spikeidx};
                nspikes = sum(spiketimes > g.GT{a}.stimon_t(i) & spiketimes < g.GT{a}.stimoff_t(i));
                g.GT{a}.spikerates = [g.GT{a}.spikerates; nspikes./(g.GT{a}.stimoff_t(i)-g.GT{a}.stimon_t(i))];
                nspikes = sum(spiketimes > g.GT{a}.stimon_t(i)-baseline_t & spiketimes < g.GT{a}.stimon_t(i));
                g.GT{a}.baselines = [g.GT{a}.baselines; nspikes./baseline_t];
            end
        end
        
        
        gtGuiHand = currentFigOpen('Grating Explorer');
        if isempty(gtGuiHand)
            g.figHand = figure;
            set(g.figHand, 'name', 'Grating Explorer', 'numbertitle', 'off')
        else
            g.figHand = gtGuiHand;
        end
        figure(g.figHand); %declare this to be the current fig
        analysisList = {'OrientCart'; 'OrientPolar'; 'SptFreq'; 'Diam'; 'Color'; 'CRF'; 'CVvsRate'; 'LFPSTA'; 'LFPSpect'};
        set(g.figHand, 'units', 'normalized', 'position', [0.0148    0.0813    0.4375    0.5250])
        g.listHand = uicontrol('style', 'popupmenu', 'string', analysisList, 'value', 1, 'callback', {@gratingAnalysis, g.figHand}, 'units', 'normalized', 'position', [0.02 0.89 .15 .1]);
        g.plotcolors = {'b', 'r', 'g', 'm', 'c', 'y'}; %for ploting multiple GT files
        set(g.figHand, 'userdata', g);
        gratingAnalysis([], [], g.figHand);
    end
    
    %setup the detection GUI
    if ~strcmpi(DTfname, 'input DT file') && ~isempty(DTfname)
        fPath = findfile(DTfname);
        if isempty(fPath)
            set(cPan.applyHand, 'string', 'DTspot file not found', 'FontSize', 11); %setting it back from 'finding files'
            set(cPan.applyHand,'BackgroundColor',[1 0.8 0.8]);
            drawnow;
            return
        end
        d.DT = dtobj(DTfname);
        [d.DT, d.grating] = stripOutGratingTrials(d.DT); %remove grating catch trials!!
        d.analParams.val.cellNum = 1;
        d.analParams.val.start = 'gabor on'; %gabor onset
        d.analParams.val.end = 'gabor off'; %gabor offset
        d.analParams.val.lowCutoff = 5; %min number of trials to compute cp
        d.analParams.val.meth = 'rate';
        [d.monk, d.cell, d.expt] = DTunpack(d.DT, d.analParams.val); %use the default params for the first analysis
        
        %initialize the dt gui
        dtGuiHand = currentFigOpen('Detection Explorer');
        if isempty(dtGuiHand);
            d.figHand = figure;
            set(d.figHand, 'name', 'Detection Explorer', 'numbertitle', 'off')
        else
            d.figHand = dtGuiHand;
        end
        figure(d.figHand) %declare this the current fig
        colorDir = reshape(d.DT.sum.exptParams.RF_colors, 3, 3)';
        colorDir(sum(abs(colorDir), 2) == 0, :) = [];
        norms = sqrt(sum(colorDir.^2, 2));
        colorDir = colorDir./repmat(norms, 1, 3);
        colorDir = num2str(colorDir, 2);
        plotTypes = {'CRF'; 'Raster'; 'PSTH'; 'CP'; 'CatchTrials'; 'Summary'};
        parseOpt = {'Flash Loc', 'Choice'};
        trigEvent = {'Flash On', 'Frame On', 'Go Sig', 'Sac Onset'};
        set(d.figHand, 'units', 'normalized', 'position', [0.4813    0.0410    0.5086    0.6475]);
        uicontrol('style', 'text', 'string', 'Color Direction', 'fontsize', 8, 'units', 'normalized', 'position', [0.015 .97 .12 .02],'backgroundColor', get(gcf, 'color'));
        d.colorDirHand = uicontrol('style', 'popupmenu', 'string', colorDir, 'value', 1, 'units', 'normalized', 'position', [.015 .867 .2 .1]);
        uicontrol('style', 'text', 'string', 'PlotType', 'fontsize', 8, 'units', 'normalized', 'position', [0.23 0.97 .1 .02], 'backgroundColor', get(gcf, 'color'));
        d.plotTypeHand = uicontrol('style', 'popupmenu', 'string', plotTypes, 'value', 1, 'units', 'normalized', 'position', [0.23 0.836 .12 .13]);
        uicontrol('style', 'text', 'string', 'Parse Opt', 'units', 'normalized', 'position', [0.37 0.97 .1 .02], 'backgroundColor', get(gcf, 'color'))
        d.parsOptHand = uicontrol('style', 'popupmenu', 'string', parseOpt, 'value', 1, 'units', 'normalized', 'position', [0.37 0.847 .12 .12]);
        uicontrol('style', 'text', 'string', 'Trigger Event', 'units', 'normalized', 'position', [0.5 0.97 .1 .02], 'backgroundColor', get(gcf, 'color'));
        d.trigEventHand = uicontrol('style', 'popupmenu', 'string', trigEvent, 'value', 1, 'units', 'normalized', 'position', [0.5 0.847 .12 .12]);
        uicontrol('style', 'text', 'string', 'Pre-Time', 'units', 'normalized', 'position', [0.64 0.97 .1 .02], 'backgroundColor', get(gcf, 'color'));
        d.preTimeHand = uicontrol('style', 'edit', 'string', num2str(200), 'units', 'normalized', 'position', [0.64 0.932 .1 .03]);
        uicontrol('style', 'text', 'string', 'Post-Time', 'units', 'normalized', 'position', [0.76 0.97 .1 .02], 'backgroundColor', get(gcf, 'color'));
        d.postTimeHand = uicontrol('style', 'edit', 'string', num2str(800), 'units', 'normalized', 'position', [0.76 0.932 .1 .03]);
        d.applyHand = uicontrol('style', 'pushbutton', 'string', 'Apply', 'callback', {@detectionAnalysis, d.figHand}, 'units', 'normalized', 'position', [0.88 0.93 .1 .04], 'fontsize', 9, 'interruptible', 'off', 'busyaction', 'cancel');
        set(d.figHand, 'userdata', d);
        detectionAnalysis([], [], d.figHand);
    end

    set(cPan.applyHand, 'string', 'Apply', 'FontSize', 11) %setting it back from 'finding files'
    drawnow
end

%
%   GRATING ANALYSIS
% A router for subsequent GT analyses. This just calls analysis modules as
% non-nested subfunctions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gratingAnalysis(h, ev, gtFigHand);
    g = get(gtFigHand, 'userdata');
    
    %the uicontrols from the CRF window linger. force them to go away
    if isfield(g,'crfHand')
        set(cell2mat(struct2cell(g.crfHand)), 'visible', 'off')
    end
    
    analysisOptions = get(g.listHand, 'string');
    analysisIdx = get(g.listHand, 'value');
    analysisToRun = analysisOptions{analysisIdx};
    try
    feval([analysisToRun, '_GTanalysisFxn'], gtFigHand);
    catch
        wtf %for debugging
    end
end


function OrientCart_GTanalysisFxn(gtFigHand)
    g = get(gtFigHand, 'userdata');
    subplot(1,1,1); %reset the figure window
    figure(g.figHand)
    cla reset,
    hold on,
    for a = 1:length(g.GT);
        orienttrials = find(g.GT{a}.protocols == 1);
        %bail if there are no orient trials
        if isempty(orienttrials)
            continue
        end
        protocolswitches = orienttrials(find(diff(orienttrials) > 1)+1);
        starts = [find(orienttrials == 1,1,'first'); protocolswitches];
        stops = [];
        for i = 1:length(starts)
            if (i == length(starts))
                stops = [stops; max(orienttrials(orienttrials > starts(i)))];
            else
                stops = [stops; max(orienttrials(orienttrials > starts(i) & orienttrials < starts(i+1)))];
            end
        end

        for i = 1:length(starts)
            trlidxs = [starts(i):stops(i)];
            [trlidxs, p1params] = gratingCheckTrials(g.GT{a}, 'orient', trlidxs);
            x = g.GT{a}.orients(trlidxs);
            y = g.GT{a}.spikerates(trlidxs);
            Ltmp = x == min(x);
            y = [y; y(Ltmp)];
            x = [x; x(Ltmp)+2*pi];

            %pp = csape(x,y,'periodic');
            xx = linspace(0,2*pi,100);
            %fit = ppval(pp,xx);
            subplot(1,length(starts),i);
            hold on;
            plot(180/pi*x,y, [g.plotcolors{a},'.']);
            %plot(180/pi*xx,fit, g.plotcolors{a});
            sf = unique(g.GT{a}.sfs(trlidxs));  % For axis label
            title(['SF: ',num2str(p1params.sf),' cyc/deg']);
            xlabel('orientation (deg)');
            ylabel('response (sp/sec)');
            set(gca,'YLim',[0 max(g.GT{a}.spikerates)]);
            mu = []; sem = [];
            for j = unique(x)'
                mu(find(j == unique(x)))= mean(y(x==j));
                sem(find(j == unique(x))) = std(y(x==j))/sqrt(sum(x==j));
            end
            errorbar(unique(x)*180/pi,mu,sem, g.plotcolors{a});
            plot([0 max(x)*180/pi], repmat(mean(g.GT{a}.baselines),1,2),[g.plotcolors{a},':']);
            cd = unique(g.GT{a}.colordirections(trlidxs,:),'rows');
            t = text(.8*max(x),0.1*max(mu),['color direction: ',num2str(p1params.color*100,'%2.1d, %2.1d, %2.1d')]);
            set(t, 'color', g.plotcolors{a});
            axis tight
            drawnow
        end

    end
end
        
function OrientPolar_GTanalysisFxn(gtFigHand)
    g = get(gtFigHand, 'userdata');
    figure(g.figHand)
    subplot(1,1,1);
    cla reset,
    mu = {};
    sem = {};
    uniqX = {};
    for a = 1:length(g.GT)
        orienttrials = find(g.GT{a}.protocols == 1);
        protocolswitches = orienttrials(find(diff(orienttrials) > 1)+1);
        starts = [find(orienttrials == 1,1,'first'); protocolswitches];
        stops = [];
        for i = 1:length(starts)
            if (i == length(starts))
                stops = [stops; max(orienttrials(orienttrials > starts(i)))];
            else
                stops = [stops; max(orienttrials(orienttrials > starts(i) & orienttrials < starts(i+1)))];
            end
        end
        
        for i = 1:length(starts)
            trlidxs = [starts(i):stops(i)];
            trlidxs = gratingCheckTrials(g.GT{a}, 'orient', trlidxs);
            x = g.GT{a}.orients(trlidxs);
            y = g.GT{a}.spikerates(trlidxs);
            Ltmp = x == min(x);
            y = [y; y(Ltmp)];
            x = [x; x(Ltmp)+2*pi];
            pp = csape(x,y,'periodic');
            xx = linspace(0,2*pi,100);
            fit = ppval(pp,xx);
            for j = unique(x)'
                mu{a,i}(find(j == unique(x)))= mean(y(x==j));
                sem{a,i}(find(j == unique(x))) = std(y(x==j))/sqrt(sum(x==j));
            end
            uniqX{a,i} = unique(x)';
        end
        
        maxEcc(a) = 0;
        for i = 1:length(starts)
            ecc = max(mu{a,i} + sem{a,i});
            maxEcc(a) = max([maxEcc(a) ecc]);
        end
    end
    
    %plotting
    [val, exptWithMax] = max(maxEcc);
    plotOrder = 0:length(g.GT);
    plotOrder(plotOrder == exptWithMax) = [];
    plotOrder(1) = exptWithMax;
    for a = 1:length(g.GT)
        exptToPlot = plotOrder(a);
        for i = 1:size(mu, 2);
            if a > size(uniqX,1) %some times there's no orient data.
                continue
            end
            polar(uniqX{exptToPlot, i},mu{exptToPlot, i}+sem{exptToPlot, i}, [g.plotcolors{exptToPlot}, ':']);
            hold on,
            h = polar(uniqX{exptToPlot, i},mu{exptToPlot, i}, g.plotcolors{exptToPlot});
            set(h,'LineWidth',2);
            polar(uniqX{exptToPlot, i},mu{exptToPlot, i}-sem{exptToPlot, i}, [g.plotcolors{exptToPlot}, ':']);
            %polar(linspace(0,2*pi,100),repmat(mean(g.GT{a}.baselines),1,100), g.plotcolors{a});
        end
    end

end


function SptFreq_GTanalysisFxn(gtFigHand)
    g = get(gtFigHand, 'userdata');
    subplot(1,1,1); %reset the figure window
    cla reset,
    hold on,
    axes
    for a = 1:length(g.GT)
        sftrials = find(g.GT{a}.protocols == 2);
        protocolswitches = sftrials(find(diff(sftrials) > 1)+1);
        starts = [min(sftrials); protocolswitches];
        stops = [];
        for i = 1:length(starts)
            if (i == length(starts))
                stops = [stops; max(sftrials(sftrials > starts(i)))];
            else
                stops = [stops; max(sftrials(sftrials > starts(i) & sftrials < starts(i+1)))];
            end
        end
        
        for i = 1:length(starts)
            trlidxs = [starts(i):stops(i)];
            [trlidxs, p2params] = gratingCheckTrials(g.GT{a}, 'sptFreq', trlidxs);
            x = g.GT{a}.sfs(trlidxs);
            y = g.GT{a}.spikerates(trlidxs);
            pp = csape(x,y,'variational');
            xx = linspace(min(g.GT{a}.sfs),max(g.GT{a}.sfs),100);
            fit = ppval(pp,xx);
            subplot(1,length(starts),i);
            hold on;
            plot(x,y,[g.plotcolors{a},'.']);
            plot(xx,fit,[g.plotcolors{a}, '-']);
            orient = unique(g.GT{a}.orients(trlidxs)); % for axis labeling
            title(['orientation: ',num2str(p2params.orient*180/pi),' deg']);
            set(gca,'XScale','log');
            xlabel('spatial frequency (cyc/deg)');
            ylabel('response (sp/sec)');
            set(gca,'YLim',[0 max(g.GT{a}.spikerates)]);
            mu = []; sem = [];
            for j = unique(x)'
                mu(find(j == unique(x)))= mean(y(x==j));
                sem(find(j == unique(x))) = std(y(x==j))/sqrt(sum(x==j));
            end
        end
        
        if ~isempty(starts)
            errorbar(unique(x),mu,sem, [g.plotcolors{a}, '-']);
            plot([min(x) max(x)], repmat(mean(g.GT{a}.baselines),1,2),[g.plotcolors{a}, ':']);
            axis tight
        end
    end
end

function Diam_GTanalysisFxn(gtFigHand)
    g = get(gtFigHand, 'userdata');
    figure(g.figHand);
    subplot(1,1,1)
    cla reset,
    hold on,
    for a = 1:length(g.GT)
        if (any(g.GT{a}.protocols == 5))
            Lprotocol = g.GT{a}.protocols == 5;
            uniquesizes = unique(g.GT{a}.diams(Lprotocol));
        else
            continue % nothing to do
        end
        mu = []; sem = [];
        for i = 1:length(uniquesizes)
            L = Lprotocol & g.GT{a}.diams == uniquesizes(i);
            tList = find(L);
            tList = gratingCheckTrials(g.GT{a}, 'diam', tList);
            L = zeros(size(L));
            L(tList) = 1; %just a hack to use gratingCheckTrials
            L = logical(L);
            mu = [mu; mean(g.GT{a}.spikerates(L))];
            sem = [sem; sqrt(var(g.GT{a}.spikerates(L))/sum(L))];
        end

        errorbar(uniquesizes, mu, sem, g.plotcolors{a});
        xlabel('Aperture diameter (deg)');
        ylabel('sp/sec');
    end
end


function Color_GTanalysisFxn(gtFigHand)
    g = get(gtFigHand, 'userdata');
    subplot(1,1,1); %reset the figure window
    cla reset,
    hold on,
    for a = 1:length(g.GT)
        Lprotocol = g.GT{a}.protocols == 4;
        if (any(Lprotocol))
            uniquecolordirs = unique(g.GT{a}.colordirections(Lprotocol,:),'rows');
            data = [];
            for i = 1:size(uniquecolordirs,1)
                L = Lprotocol &...
                    g.GT{a}.colordirections(:,1) == uniquecolordirs(i,1) & ...
                    g.GT{a}.colordirections(:,2) == uniquecolordirs(i,2) & ...
                    g.GT{a}.colordirections(:,3) == uniquecolordirs(i,3);
                tList = find(L);
                tList = gratingCheckTrials(g.GT{a}, 'color', tList);
                L = zeros(size(L));
                L(tList) = 1; %just a hack to use gratingCheckTrials
                L = logical(L);
                data = [data; mean(g.GT{a}.spikerates(L)) sqrt(var(g.GT{a}.spikerates(L))/sum(L))];
                nTrials(i) = sum(L);
            end
            % Plotting
            basis = [1 1 0; 1 -1 0; 0 0 1];
            DKLs = (basis*uniquecolordirs')';
            coordinates = sign(DKLs);
            peakfr = max(data(:,1));
            labels = {'L+M','L-M';'S','L-M';'S','L+M'};
            for i = 1:3
                subplot(1,3,i);
                polar(0,peakfr, g.plotcolors{a});
                hold on;
                polar(linspace(0,2*pi,100), repmat(mean(g.GT{a}.baselines),1,100), g.plotcolors{a});
                if a == 1
                    text(0,2*peakfr,labels{i,1},'HorizontalAlignment','center');
                    text(.7*peakfr,.15*peakfr,labels{i,2});
                end
            end

            for i = 1:size(uniquecolordirs,1)
                subplot(1,3,1);
                if (coordinates(i,3) == 0)
                    theta = atan2(coordinates(i,1),coordinates(i,2));
                    polarerrbar(theta, data(i,1), data(i,2), g.plotcolors{a});
                end
                subplot(1,3,2);
                if (coordinates(i,1) == 0)
                    theta = atan2(coordinates(i,3),coordinates(i,2));
                    polarerrbar(theta, data(i,1), data(i,2), g.plotcolors{a});
                end
                subplot(1,3,3);
                if (coordinates(i,2) == 0)
                    theta = atan2(coordinates(i,3),coordinates(i,1));
                    polarerrbar(theta, data(i,1), data(i,2), g.plotcolors{a});
                end
            end
            for i = 1:3
                subplot(1,3,i);
                polar(linspace(0,2*pi,100),repmat(mean(g.GT{a}.baselines),1,100), [g.plotcolors{a}, ':']);
            end
            
            
            if length(g.GT) == 1;
                %charlie's addition for determining the pref color and the CSI
                [dummy, lumIdx] = ismember([1 1 0], sign(uniquecolordirs), 'rows');
                lumResp = data(lumIdx, 1);
                isoLumBasis = [1, -1, 0; 0, 0, 1; -0.0636, 0.0636, 0.45; -0.0636, 0.0636, -0.45];
                isoLumNorms = sqrt(sum(isoLumBasis.^2,2));
                isoLumUnits = isoLumBasis ./ repmat(isoLumNorms, 1, 3);
                uniqueColorNorms = sqrt(sum(uniquecolordirs.^2,2));
                uniqueColorUnits = uniquecolordirs ./ repmat(uniqueColorNorms, 1, 3);
                for a = 1:size(isoLumBasis, 1);
                    clrIdx = abs(isoLumUnits(a,:) * uniqueColorUnits') > (1-eps);
                    isoLumResp(a) = data(clrIdx,1);
                end
                [maxIsoLumResp, maxIdx] = max(isoLumResp);
                CSI = maxIsoLumResp./lumResp;
                prefColorIdx = abs(isoLumUnits(maxIdx,:) * uniqueColorUnits') > (1-eps);
                prefColor = uniquecolordirs(prefColorIdx, :);

                % Greg's linear model stuff
                % This should really be a weighted regression, but for now just
                % fitting the mean responses.
                signmat = 2*fullfact(repmat(2,size(uniquecolordirs,1),1))-3;
                normedcolordirs = mkbasis(uniquecolordirs')';
                normedresponses = data(:,1)./sqrt(sum(transpose(uniquecolordirs.^2)))';
                err = Inf;
                for i = 1:size(signmat,1)
                    b = regress(normedresponses,normedcolordirs.*repmat(signmat(i,:)',1,3));
                    pred = (normedcolordirs.*repmat(signmat(i,:)',1,3))*b;
                    if (norm(pred-normedresponses) < err)
                        linmodelprefcolor = b;
                        err = norm(pred-normedresponses);
                    end
                end
                linmodelprefcolor = mkbasis(linmodelprefcolor);
                if (sum(sign(linmodelprefcolor)) <= -2)
                    linmodelprefcolor = -linmodelprefcolor;
                end
                textax = axes('position',[.25 .2, .5, .1]);
                set(textax,'DefaultTextHorizontalAlignment','center')
                text(.5,.9, sprintf('Pref Color: [%s] \n CSI: %.2f', num2str(prefColor, 3), CSI));
                text(.5,.4, sprintf('Cone weights from linear model: [%s]', num2str(linmodelprefcolor', 3)));
                text(.5, .1, sprintf('Min Num Trials: %d', min(nTrials)));
                set(textax,'Visible','off');
            end
        end
    end

    
    % Support function for grating analysis
    function polarerrbar(theta, x, e, color)
        p(1)=polar([theta theta], [x-e x+e], [color, '-']);
        p(2)=polar([theta+pi theta+pi], [x-e x+e], [color, '-']);
        p(3)=polar(theta, x, [color, '.']);
        p(4)=polar(theta+pi, x, [color, '.']);
        set(p, 'markersize', 15, 'linewidth', 3)
    end
end

function CRF_GTanalysisFxn(gtFigHand)
    g = get(gtFigHand, 'userdata');
    
    %pick a color direction to analyze
    g.crfHand.color = uicontrol('style', 'popupmenu', 'string', num2str([1, -1, 0; 0, 0, 1]), 'callback', @plotGTcrf, 'value', 1, 'units', 'normalized', 'position', [0.4 0.87 .12 .12]);
    g.crfHand.text = uicontrol('style', 'text', 'string', 'Color Dir', 'units', 'normalized', 'position', [0.28 0.94 .1 .04]); 
    set(gtFigHand, 'userdata', g); %so that the uicontrols can be deleted later
    plotGTcrf(1, 1); %call the subFxn with some dummy arguments.
    
    function plotGTcrf(ev, h)
        subplot(1,1,1); %reset the figure window
        cla reset,
        menuVal = get(g.crfHand.color, 'value');
        colorDirs = get(g.crfHand.color, 'string');
        colorToPlot = str2num(colorDirs(menuVal,:));
        legendText = {};
        linetype = {'-', '--', ':'};
        hold on,
        for exptNum = 1:length(g.GT);
            LMS = g.GT{exptNum}.colordirections;
            l_color = ismember(sign(LMS), colorToPlot, 'rows');
            l_crf = g.GT{exptNum}.protocols == 7;
%             if sum(l_crf) ==0
%                 continue %nothing to do.
%             end
            tmp = unique(LMS((l_color & l_crf), :), 'rows');
            contrasts = zeros(size(tmp,1)+1, 3);
            contrasts(2:size(tmp,1)+1, :) = tmp; %add the zero contrast cond.
            norms = sqrt(sum(contrasts.^2, 2));
            gtmu = []; gtsem = [];
            for a = 1:size(contrasts,1)
                l_contrasts = softEq(contrasts(a,:), LMS, 5, 'rows');
                if a == 1;
                    l_trials = l_contrasts;
                else
                    l_trials = l_crf & l_color & l_contrasts;
                end
                gtmu(a) = mean(g.GT{exptNum}.spikerates(l_trials));
                gtsem(a) = std(g.GT{exptNum}.spikerates(l_trials)) ./ sqrt(sum(l_trials));
                nGTTrials(exptNum, a) = sum(l_trials);
            end
            
            % set up the plot
            if length(g.GT) < 3;
                if sum(colorToPlot) == 1;
                    lineColor = 'b';
                elseif sum(colorToPlot) == 0;
                    lineColor = 'r';
                end
                linestyle = linetype{exptNum};
            else
                lineColor = g.plotcolors{exptNum};
                linestyle = '-';
            end
            if all(~isnan(gtmu))
                errorbar(norms, gtmu, gtsem, [lineColor, '.', linestyle])
                legendText{exptNum} = sprintf('GT run: %d', exptNum);
            end
        end
         
         %determine if there is any DT data (from the same color dir) to plot.
         dtFigHand = currentFigOpen('Detection Explorer');
         if ~isempty(dtFigHand)
             d = get(dtFigHand, 'userdata');
             clrIdx = ismember(sign(d.expt.standColors), colorToPlot, 'rows');
             if any(clrIdx)
                 dtmu = cellfun(@mean, d.cell.crfIn{clrIdx});
                 dtsigma = cellfun(@std, d.cell.crfIn{clrIdx});
                 nTrials = cellfun(@length, d.cell.crfIn{clrIdx});
                 dtsem = dtsigma./sqrt(nTrials);
                 errorbar(d.expt.norms{clrIdx}, dtmu, dtsem, 'k.:')
                 legendText{end+1} = 'Detection';
             end
         end
         
         %tidy up the plot
         set(gca, 'xscale', 'log')
         axis tight
         hold off
         set(gtFigHand, 'userdata', g); %so that the uicontrols can be deleted later
         hold off
         if ~isempty(legendText)
             legend(legendText, 'location', 'northwest')
             xlabel(sprintf('Min Num Trials = [%s] (one N per GT expt)', num2str(min(nGTTrials, [], 2)')))
         end
    end
end

function LFPSpect_GTanalysisFxn(gtFigHand)
    
end

%
%   DETECTION ANALYSIS
%
% gets evaluated when the apply button on the Detection analysis sub-pannel
% is pressed. Basically just redirects to the appropriate analysis module
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function detectionAnalysis(h, ev, dtFigHand);
    d = get(dtFigHand, 'userdata');
    
    %the uicontrols from the summary window linger. force them to go away
    if isfield(d.analParams,'hand')
        set(cell2mat(struct2cell(d.analParams.hand)), 'visible', 'off')
    end
    
    %call the appropriate sub function based on what type of plot the user
    %wants to compute.
    plotType = get(d.plotTypeHand, 'string');
    plotTypeIdx = get(d.plotTypeHand, 'value');
    plotType = plotType{plotTypeIdx};
    endString = '_analysisFxn';
    feval([plotType, endString], dtFigHand);
end

%
%   CRF ANALYSIS MODULE FOR (DT)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CRF_analysisFxn(dtFigHand)
    d = get(dtFigHand, 'userdata');
    
    %first for in the RF
    whichColorDir = get(d.colorDirHand, 'value');
    statAvgIn = cellfun(@mean, d.cell.crfIn{whichColorDir});
    statStd = cellfun(@std, d.cell.crfIn{whichColorDir});
    statN = cellfun(@length, d.cell.crfIn{whichColorDir});
    stdErrIn = statStd./sqrt(statN);
    
    %now for out of the RF
    statAvgOut = cellfun(@mean, d.cell.crfOut{whichColorDir});
    statStd = cellfun(@std, d.cell.crfOut{whichColorDir});
    statN = cellfun(@length, d.cell.crfOut{whichColorDir});
    stdErrOut = statStd./sqrt(statN);
    
    %plotting
    figure(d.figHand); % assert this as the current fig
    subplot(1, 1, 1)% set things back to just one axis
    hold on,
    set(gca, 'position', [0.15 0.13 .7 .75])
    errorbar(d.expt.norms{whichColorDir}, statAvgIn, stdErrIn, 'b');
    errorbar(d.expt.norms{whichColorDir}, statAvgOut, stdErrOut, 'k');
    set(gca, 'xlim', [d.expt.norms{whichColorDir}(2).*.97 max(d.expt.norms{whichColorDir}).*1.02]);
    set(gca, 'xscale', 'log')
    set(gca, 'xticklabel', (10.^str2num(get(gca, 'xticklabel')).*100))
    set(get(gca, 'ylabel'), 'string', 'Avg Trial Stat', 'fontsize', 12)
    set(get(gca, 'xlabel'), 'string', 'Percent Cone Contrast', 'fontsize', 12)
    set(get(gca, 'children'), 'linewidth', 2)
    legend('In RF', 'Out of RF')
    legend boxoff
    hold off
    
end


%
%   RASTER ANALYSIS MODULE (DT)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Raster_analysisFxn(dtFigHand)
    d = get(dtFigHand, 'userdata');
    subplot(1, 1, 1) %clear out the old data.
    colorDirInd = strmatch('color_dir', [d.DT.sum.trialFields(1,:)]);
    flashXInd = strmatch('flash_x', [d.DT.sum.trialFields(1,:)]);
    flashYInd = strmatch('flash_y', [d.DT.sum.trialFields(1,:)]);
    correctInd = strmatch('correct', [d.DT.sum.trialFields(1,:)]);
    cntrstLevInd = strmatch('cntrst_lev', [d.DT.sum.trialFields(1,:)]);
    
    %unpack the user settings
    colorDir = get(d.colorDirHand, 'value');
    parseOpt = get(d.parsOptHand, 'string');
    parseIdx = get(d.parsOptHand, 'value');
    parseOpt = parseOpt{parseIdx};
    trigOpt = get(d.trigEventHand, 'string');
    trigIdx = get(d.trigEventHand, 'value');
    trigOpt = trigOpt{trigIdx};
    preTime = str2num(get(d.preTimeHand, 'string'))./1000;
    postTime = str2num(get(d.postTimeHand, 'string'))./1000;
    
    %determine which index in the stro.trial field to extract the time zero
    %value. This will be specified by the trigOpt.
    flashOnIdx = strmatch('rep_flash_on', [d.DT.sum.trialFields(1,:)]);
    frameOnIdx = strmatch('frame_on', [d.DT.sum.trialFields(1,:)]);
    targOnIdx = strmatch('targ_on', [d.DT.sum.trialFields(1,:)]);
    choiceTimeIdx = strmatch('choice_time', [d.DT.sum.trialFields(1,:)]);
    numFramesIdx = strmatch('numframes', [d.DT.sum.trialFields(1,:)]);
    switch lower(trigOpt)
        case 'flash on'
            timeZeroIdx = flashOnIdx;
        case 'frame on'
            timeZeroIdx = frameOnIdx;
        case 'go sig'
            timeZeroIdx = targOnIdx;
        case 'sac onset'
            timeZeroIdx = choiceTimeIdx;
    end
             
    % loop twice, one for each parse option.
    maxNumTrials = 0;
    numContrasts = length(unique(d.DT.trial(:,cntrstLevInd)));
    axAddress = reshape(1:numContrasts.*2, 2, numContrasts)';
    subPlotCounter = 0;
    for a = 1:2
        l_color = d.DT.trial(:,colorDirInd) == colorDir;
        switch lower(parseOpt)
            case 'choice'
                rfx = d.DT.sum.exptParams.rf_x;
                rfy = d.DT.sum.exptParams.rf_y;
                inRF = (d.DT.trial(:,flashXInd) == rfx) & (d.DT.trial(:,flashYInd) == rfy);
                corrects = d.DT.trial(:, correctInd);
                T1choices = zeros(size(d.DT.trial, 1), 1);
                T1choices(corrects & inRF) = 1; %in the RF and correct
                T1choices(~corrects & ~inRF) = 1; %incorrect but not in the RF
                if a == 1
                    l_parseOpt = T1choices;
                elseif a == 2
                    l_parseOpt = ~T1choices;
                end
            case 'flash loc'
                rfx = d.DT.sum.exptParams.rf_x;
                rfy = d.DT.sum.exptParams.rf_y;
                inRF = (d.DT.trial(:,flashXInd) == rfx) & (d.DT.trial(:,flashYInd) == rfy);
                if a == 1
                    l_parseOpt = inRF;
                elseif a == 2
                    l_parseOpt = ~inRF;
                end
        end
        
        %loop through the trials (cntrst by cntrst) making raster plots
        for cntrst = 1:numContrasts;
            l_contrast = d.DT.trial(:, cntrstLevInd) == cntrst;
            if cntrst == 1; %treat the 'zero' contrast separately
                tList = (l_contrast & l_parseOpt);
            else
                tList = (l_color & l_contrast & l_parseOpt);
            end
            
            %place the subplot window in a reasonable place. then do the plotting\
            subPlotCounter = subPlotCounter+1;
            spHand(subPlotCounter) = subplot(numContrasts, 2, axAddress(cntrst, a));            
            height = (0.8./numContrasts) .* 0.8;
            currentPos = get(gca, 'position');
            currentPos(end) = height;
            currentPos(2) = currentPos(2) - .02;
            set(gca, 'position', currentPos)
            hold on,
            zeroTimes = mat2cell(d.DT.trial(tList, timeZeroIdx), ones(sum(tList),1),1);
            spikeTimes = cellfun(@(x,y)(x-y), d.DT.ras(tList,d.expt.cellNum), zeroTimes, 'UniformOutput', 0);
            counter = mat2cell([0:sum(tList)-1]', ones(sum(tList), 1), 1);
            cellfun(@(x, y)plot([x, x]', [zeros(1,length(x))+y; [ones(1, length(x)).*0.8 + y]], 'k'), spikeTimes, counter);
            
            %mark the relavant trial events
            plot([d.DT.trial(tList, frameOnIdx) - d.DT.trial(tList, timeZeroIdx)], [0:sum(tList)-1]+0.5, 'r.');%frame on
            plot([d.DT.trial(tList, flashOnIdx) - d.DT.trial(tList, timeZeroIdx)], [0:sum(tList)-1]+0.5, 'c.');%flash onset
            stimOffTime = [d.DT.trial(tList, flashOnIdx) - d.DT.trial(tList, timeZeroIdx)] + [d.DT.trial(tList, numFramesIdx)./d.DT.sum.exptParams.frame_rate];
            plot(stimOffTime, [0:sum(tList)-1], 'b.'); %flash offset
            plot([d.DT.trial(tList, targOnIdx) - d.DT.trial(tList, timeZeroIdx)], [0:sum(tList)-1], 'g.'); %go signal
            plot([d.DT.trial(tList, choiceTimeIdx) - d.DT.trial(tList, timeZeroIdx)], [0:sum(tList)-1], 'm.') %sac time
            hold off,
            
            %keep track of the nTrials and use this to scale all the axes
            maxNumTrials = max(maxNumTrials, sum(tList));
            

        end
    end
    
    %now loop through again, this time standardizing the axes
    subPlotCounter = 0;
    for a = 1:2;
        for cntrst = 1:numContrasts;
            subPlotCounter = subPlotCounter+1;
            set(spHand(subPlotCounter), 'ytick', [], 'xtick', [-preTime, 0, postTime], 'xticklabel', [], 'box', 'off')
            set(spHand(subPlotCounter), 'ylim', [0 maxNumTrials+1], 'xlim', [-preTime postTime]);
            if a == 1;
                set(get(spHand(subPlotCounter), 'ylabel'), 'string', num2str(cntrst-1));
            end
            if cntrst == numContrasts;
                set(spHand(subPlotCounter), 'xticklabel', [-preTime, 0, postTime]);
                set(get(spHand(subPlotCounter), 'xlabel'), 'string', sprintf('Parse Opt %d', a));
            end
        end
    end
end


%
%   PSTH ANALYSIS MODULE (DT)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PSTH_analysisFxn(dtFigHand)
    d = get(dtFigHand, 'userdata');
    subplot(1, 1, 1) %clear out the old data.
    colorDirInd = strmatch('color_dir', [d.DT.sum.trialFields(1,:)]);
    flashXInd = strmatch('flash_x', [d.DT.sum.trialFields(1,:)]);
    flashYInd = strmatch('flash_y', [d.DT.sum.trialFields(1,:)]);
    correctInd = strmatch('correct', [d.DT.sum.trialFields(1,:)]);
    cntrstLevInd = strmatch('cntrst_lev', [d.DT.sum.trialFields(1,:)]);
    
    %unpack the user settings
    colorDir = get(d.colorDirHand, 'value');
    parseOpt = get(d.parsOptHand, 'string');
    parseIdx = get(d.parsOptHand, 'value');
    parseOpt = parseOpt{parseIdx};
    trigOpt = get(d.trigEventHand, 'string');
    trigIdx = get(d.trigEventHand, 'value');
    trigOpt = trigOpt{trigIdx};
    preTime = str2num(get(d.preTimeHand, 'string'))./1000;
    postTime = str2num(get(d.postTimeHand, 'string'))./1000;
    
    %determine which index in the stro.trial field to extract the time zero
    %value. This will be specified by the trigOpt.
    flashOnIdx = strmatch('rep_flash_on', [d.DT.sum.trialFields(1,:)]);
    frameOnIdx = strmatch('frame_on', [d.DT.sum.trialFields(1,:)]);
    targOnIdx = strmatch('targ_on', [d.DT.sum.trialFields(1,:)]);
    choiceTimeIdx = strmatch('choice_time', [d.DT.sum.trialFields(1,:)]);
    numFramesIdx = strmatch('numframes', [d.DT.sum.trialFields(1,:)]);
    switch lower(trigOpt)
        case 'flash on'
            timeZeroIdx = flashOnIdx;
        case 'frame on'
            timeZeroIdx = frameOnIdx;
        case 'go sig'
            timeZeroIdx = targOnIdx;
        case 'sac onset'
            timeZeroIdx = choiceTimeIdx;
    end
             
    % loop twice, one for each parse option.
    maxRate = 0;
    numContrasts = length(unique(d.DT.trial(:,cntrstLevInd)));
    axAddress = reshape(1:numContrasts.*2, 2, numContrasts)';
    for a = 1:2
        l_color = d.DT.trial(:,colorDirInd) == colorDir;
        switch lower(parseOpt)
            case 'choice'
                rfx = d.DT.sum.exptParams.rf_x;
                rfy = d.DT.sum.exptParams.rf_y;
                inRF = (d.DT.trial(:,flashXInd) == rfx) & (d.DT.trial(:,flashYInd) == rfy);
                corrects = d.DT.trial(:, correctInd);
                T1choices = zeros(size(d.DT.trial, 1), 1);
                T1choices(corrects & inRF) = 1; %in the RF and correct
                T1choices(~corrects & ~inRF) = 1; %incorrect but not in the RF
                if a == 1
                    l_parseOpt = T1choices;
                elseif a == 2
                    l_parseOpt = ~T1choices;
                end
            case 'flash loc'
                rfx = d.DT.sum.exptParams.rf_x;
                rfy = d.DT.sum.exptParams.rf_y;
                inRF = (d.DT.trial(:,flashXInd) == rfx) & (d.DT.trial(:,flashYInd) == rfy);
                if a == 1
                    l_parseOpt = inRF;
                elseif a == 2
                    l_parseOpt = ~inRF;
                end
        end
        
        %loop through the trials (cntrst by cntrst) making raster plots
        for cntrst = 1:numContrasts;
            l_contrast = d.DT.trial(:, cntrstLevInd) == cntrst;
            if cntrst == 1; %treat the 'zero' contrast separately
                tList = (l_contrast & l_parseOpt);
            else
                tList = (l_color & l_contrast & l_parseOpt);
            end

            %place the subplot window in a reasonable place. then do the plotting
            subplot(numContrasts, 2, axAddress(cntrst, a));
            hold on,
            zeroTimes = mat2cell(d.DT.trial(tList, timeZeroIdx), ones(sum(tList),1),1);
            spikeTimes = cellfun(@(x,y)(x-y), d.DT.ras(tList,d.expt.cellNum), zeroTimes, 'UniformOutput', 0);
            binSize = 0.025;
            edges = -preTime:binSize:postTime;
            edges = mat2cell(repmat(edges, sum(tList), 1), ones(sum(tList),1), length(edges));
            countsPerBin = cellfun(@(x, y)histc(x, y), spikeTimes, edges, 'uniformoutput', 0);
            countsPerBin = cellfun(@(x)(x(:)'), countsPerBin, 'uniformoutput', 0); %allign for vertcat
            countsPerBin = vertcat(countsPerBin{:});
            avgCountsPerBin = mean(countsPerBin, 1);
            avgRatePerBin = avgCountsPerBin ./ binSize;
            if (~isempty(avgCountsPerBin))
                bar(edges{1}, avgRatePerBin, 1);
            end
            %keep track of the nTrials and use this to scale all the axes
            maxRate = max(max(avgRatePerBin), maxRate);
        end
    end
    
    %now loop through again, this time standardizing the axes
    for a = 1:2;
        for cntrst = 1:numContrasts;
            subplot(numContrasts, 2, axAddress(cntrst, a))
            hold on,
            plot([0, 0], [0, maxRate], 'k', 'linewidth', 2);
            hold off,
            set(gca, 'ytick', [0, round(maxRate./2), maxRate], 'yticklabel', [], 'xtick', [-preTime, 0, postTime], 'xticklabel', [], 'box', 'off')
            ylim([0 maxRate]);
            xlim([-preTime postTime]);
            if a == 1;
                ylabel(num2str(cntrst-1));
            end
            if cntrst == numContrasts;
                set(gca, 'xticklabel', [-preTime, 0, postTime]);
                xlabel(sprintf('Parse Opt %d', a));
                if a == 1;
                    set(gca, 'yticklabel', [0, round(maxRate./2), maxRate])
                end
            end
            %make sure the size of the subplot is o.k
            height = (0.8./numContrasts) .* 0.8;
            currentPos = get(gca, 'position');
            currentPos(end) = height;
            currentPos(2) = currentPos(2) - .02;
            set(gca, 'position', currentPos)
        end
    end
end


%
%   CHOICE PROBABILITY MODULE (DT)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CP_analysisFxn(dtFigHand)
    d = get(dtFigHand, 'userdata');
    colorDir = get(d.colorDirHand, 'value');
    figure(d.figHand) %assert as the current fig
    subplot(1,1,1); %clear out the old axes
    hold on,
    plot((1:length(d.expt.norms{colorDir, 1})), d.cell.cp.in{colorDir, 1}, 'g.', 'markersize', 15) 
    plot((1:length(d.expt.norms{colorDir, 1})), d.cell.cp.out{colorDir, 1}, 'k.', 'markersize', 15)
    plot([0.5 length(d.expt.norms{colorDir, 1})+0.5], [0.5 0.5], 'k:')
    set(gca, 'xtick', [1:length(d.expt.norms{colorDir, 1})], 'xticklabel', round(d.expt.norms{colorDir, 1}.*10000)./100)
    set(gca, 'units', 'normalized', 'position', [0.13 0.49 0.775 0.4])
    ylim([0 1])
    xlim([0 length(d.expt.norms{colorDir, 1})+1])
    xlabel('Percent Cone Contrast')
    ylabel('Choice Probability')
    legend('Stim In RF', 'Stim Out of RF')
    legend boxoff
    hold off
       
    %print out the pooled CP values.
    bkgndColor = get(d.figHand, 'color');
    ax = axes('position', [0.05 0.05 0.9 0.35]);
    set(ax, 'xcolor', bkgndColor, 'ycolor', bkgndColor, 'color', bkgndColor)
    ylim([0 1]) %makes placing text easier
    xlim([0 1])
    txt = text(0.1, 0.8, sprintf('Pooled Choice Probability \n\n* In RF: %.3f \n* Out of RF: %.3f',...
        d.cell.cp.poolConIn.val(colorDir, 1), d.cell.cp.poolConOut.val(colorDir,1)));
end


%
%   SUMMARY SLIDE FOR DETECTION (DT)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Summary_analysisFxn(dtFigHand);
    d = get(dtFigHand, 'userdata');
    figure(d.figHand) %assert this to be the current fig
    subplot(1,1,1);
    bkgndColor = get(d.figHand, 'color');
    set(gca, 'position', [0.05 0.05 0.9 0.85], 'xcolor', bkgndColor, 'ycolor', bkgndColor, 'color', bkgndColor)
    ylim([0 1]) %makes placing text easier
    xlim([0 1]) 
    
    %Num Trials
    cntrstLevInd = strmatch('cntrst_lev', [d.DT.sum.trialFields(1,:)]);
    colorDirInd = strmatch('color_dir', [d.DT.sum.trialFields(1,:)]);
    colorDirs = get(d.colorDirHand, 'string');
    numContrasts = length(unique(d.DT.trial(:,cntrstLevInd)));
    nTrials = nan(numContrasts, size(colorDirs,1));
    for clr = 1:size(colorDirs,1)
        l_color = d.DT.trial(:,colorDirInd) == clr;
        for cntrst = 1:numContrasts
            l_cntrst = d.DT.trial(:,cntrstLevInd) == cntrst;
            if cntrst == 1;
                nTrials(cntrst, clr) = sum(l_cntrst);
            else
                nTrials(cntrst, clr) = sum(l_cntrst&l_color);
            end
        end
    end
    nTrialTxt = text(0.05, .94, sprintf('Num Trials:\n    min: %d\n    max: %d', min(nTrials(:)), max(nTrials(:))));
    
    %Choice bias, measured for Zero contrast trials.
    flashXInd = strmatch('flash_x', [d.DT.sum.trialFields(1,:)]);
    flashYInd = strmatch('flash_y', [d.DT.sum.trialFields(1,:)]);
    correctInd = strmatch('correct', [d.DT.sum.trialFields(1,:)]);
    rfx = d.DT.sum.exptParams.rf_x;
    rfy = d.DT.sum.exptParams.rf_y;
    inRF = (d.DT.trial(:,flashXInd) == rfx) & (d.DT.trial(:,flashYInd) == rfy);
    corrects = d.DT.trial(:, correctInd);
    T1choices = zeros(size(d.DT.trial, 1), 1);
    T1choices(corrects & inRF) = 1; %in the RF and correct
    T1choices(~corrects & ~inRF) = 1; %incorrect but not in the RF
    tList = d.DT.trial(:, cntrstLevInd) == 1;
    ChoicesAtZero(1) = sum(tList & T1choices);
    ChoicesAtZero(2) = sum(tList & ~T1choices);
    [val, prefChoice] = max(ChoicesAtZero);
    p = 1-binocdf(ChoicesAtZero(prefChoice), sum(ChoicesAtZero), 0.5);
    biasTxt = text(0.05, 0.78, sprintf('Bias: \n    T1 choices: %d\n    T2 choices %d\n    p = %.4f', ChoicesAtZero(1), ChoicesAtZero(2), p));
    
    % ANOVA on noise trials vs. outRF trials
    text(0.05, 0.63, 'Comparison of Zero Contrast to Out RF Responses Kruskal-Wallis:');
    for clrDir = 1:length(d.cell.crfOut)
        nSamples = cellfun(@(x)size(x,2), d.cell.crfOut{clrDir}(:)); %nTrials at each contrast
        group = cellfun(@(x)ones(x,1), mat2cell(nSamples, ones(length(nSamples),1)), 'uniformoutput', 0); %making a 'group' vec for anova1
        group = cellfun(@(x,y)(x.*y), group, mat2cell([1:length(group)]', ones(length(group),1)), 'uniformoutput', 0);
        group = vertcat(group{:});
        responses = horzcat(d.cell.crfOut{clrDir}{:})';
        p(clrDir) = kruskalwallis(responses, group, 'off');
        text(0.05, 0.63-0.04*clrDir, sprintf('    [%s]   p = %.4f', colorDirs(clrDir,:), p(clrDir)));
    end

    % Threshold ratios
    text(0.05, 0.46,'Threshold Ratios:'); 
    for clrDir = 1:length(d.cell.crfOut)
        TR = d.cell.alpha(clrDir)./d.monk.alpha(clrDir);
        text(0.05, 0.46-0.04*clrDir, sprintf('    [%s] => %.3f', colorDirs(clrDir,:), TR));
    end
    
    % User input for manipulating the default params for DT analysis
    text(0.05, 0.29, 'Analysis Parameters:');
    text(0.05, 0.25, '    Cell Number:');
    text(0.05, 0.20, '    Start Time(sec):');
    text(0.05, 0.15, '    End Time(sec):');
    text(0.05, 0.10, '    Min Num Trials for CP:');
    text(0.05, 0.05, '    Analysis Method:');
    d.analParams.hand.cellNum = uicontrol('style', 'edit', 'string', num2str(d.analParams.val.cellNum), 'callback', @updateAnalParams, 'units', 'normalized', 'position', [.28, 0.245, .1, .03]);
    d.analParams.hand.start = uicontrol('style', 'edit', 'string', num2str(d.analParams.val.start), 'callback', @updateAnalParams, 'units', 'normalized', 'position', [.28, 0.205, .1, .03]);
    d.analParams.hand.end = uicontrol('style', 'edit', 'string', num2str(d.analParams.val.end), 'callback', @updateAnalParams, 'units', 'normalized', 'position', [.28, 0.165, .1, .03]);
    d.analParams.hand.lowCutoff = uicontrol('style', 'edit', 'string', num2str(d.analParams.val.lowCutoff), 'callback', @updateAnalParams, 'units', 'normalized', 'position', [.35, 0.124, .1, .03]);
    d.analParams.hand.meth = uicontrol('style', 'edit', 'string', num2str(d.analParams.val.meth), 'callback', @updateAnalParams, 'units', 'normalized', 'position', [.30, 0.075, .1, .03]);
    d.analParams.hand.refresh = uicontrol('style', 'pushbutton', 'string', 'Refresh DT Data', 'callback', {@refreshDTdata}, 'units', 'normalized', 'position', [0.65 0.2 .2 .04], 'fontsize', 10, 'interruptible', 'off', 'busyaction', 'cancel');
    set(d.figHand, 'userdata', d);
    
    function updateAnalParams(h, ev)
        d.analParams.val.cellNum = str2num(get(d.analParams.hand.cellNum, 'string'));
        d.analParams.val.start = str2num(get(d.analParams.hand.start, 'string'));
        d.analParams.val.end = str2num(get(d.analParams.hand.end, 'string'));
        d.analParams.val.lowCutoff = str2num(get(d.analParams.hand.lowCutoff, 'string'));
        d.analParams.val.meth = get(d.analParams.hand.meth, 'string');
        set(d.figHand, 'userdata', d);
    end

    function refreshDTdata(h, ev)
        [d.monk, d.cell, d.expt] = DTunpack(d.DT, d.analParams.val); %use the most recent analysis parameters
        set(d.figHand, 'userdata', d)
        detectionAnalysis([], [], d.figHand);
    end
end

%
% GRATING CATCH TRIAL MODULE (DT)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CatchTrials_analysisFxn(dtFigHand);
    d = get(dtFigHand, 'userdata');
    figure(d.figHand) %assert this to be the current fig
    subplot(1,1,1); %clear out what was plotted previously
    for a = 1:length(d.DT.sum.rasterCells)
        wf_list(a) = strcmpi('_wf', d.DT.sum.rasterCells{a}(end-2:end));
    end
    waveFormColumns = find(wf_list);
    if isempty(waveFormColumns)
        bkgndColor = get(d.figHand, 'color');
        set(gca, 'position', [0.05 0.05 0.9 0.85], 'xcolor', bkgndColor, 'ycolor', bkgndColor, 'color', bkgndColor)
        ylim([0 1]) %makes placing text easier
        xlim([0 1])
        if isempty(d.grating.trial) %if there are no catch trials
            t = text(.3, .5, 'No Grating Catch Trials Present');
        else
            t = text(.16, .55, 'Catch Trials Present But No Waveform Data Saved');    
        end
        set(t, 'fontsize', 13);
        return
    end
    wfIdx = waveFormColumns(d.expt.cellNum);
    
    
    %plot all the waveforms from all catch trials alongside the trial mean
    %wave form.
    storeRate = 1./d.DT.sum.waves.storeRates{1};
    wfs = vertcat(d.grating.ras{:,wfIdx});
    maxVal = max(max(wfs));
    minVal = min(min(wfs));
    meanWfsByTrial = cellfun(@(x) mean(x,1), d.grating.ras(:, wfIdx), 'uniformoutput', false);
    meanWfsByTrial = vertcat(meanWfsByTrial{:});
    nSamples = length(meanWfsByTrial(1,:));
    timeVector = [0 : storeRate : storeRate.*(nSamples-1)]; %in seconds
    timeVector = timeVector .* 10^6; %in usec
    subplot(1,2,1)
    hold on,
    plot(timeVector, wfs', 'k')
    plot(timeVector, meanWfsByTrial', 'r', 'linewidth', 2)
    ylim([minVal, maxVal])
    xlim([0, timeVector(end)]);
    hold off
    xlabel('Time (usec)');
    ylabel('volts?!?!?!');
    
    %now plot the mean firing rate in response to the grating (as a
    %function of trial number);
    nFramesInd = strmatch('numframes', [d.DT.sum.trialFields(1,:)]);
    flashOnInd = strmatch('flash_on', [d.DT.sum.trialFields(1,:)]);
    trialDurations = d.grating.trial(:,nFramesInd) ./ d.DT.sum.exptParams.frame_rate;
    flashStartTime = mat2cell(d.grating.trial(:,flashOnInd), ones(length(d.grating.trial(:,flashOnInd)),1));
    flashEndTime = mat2cell(d.grating.trial(:,flashOnInd)+trialDurations, ones(length(d.grating.trial(:,flashOnInd)),1));
    trialSpikes = cellfun(@(ras, start, stop) ((ras>start)&(ras<=stop)),d.grating.ras(:,d.expt.cellNum), flashStartTime, flashEndTime, 'uniformoutput', false);
    trialRates = cellfun(@(spikes, durs) (sum(spikes)./durs), trialSpikes, mat2cell(trialDurations, ones(length(trialDurations),1)));
    subplot(1,2,2)
    flashStartTime = cell2mat(flashStartTime)./60; %in minutes
    plot(flashStartTime, trialRates, 'bo-', 'linewidth', 2)
    xlim([min(flashStartTime)*0.98, max(flashStartTime)*1.02])
    ylim([0, max(trialRates).*1.1])
    xlabel(sprintf('Catch Trials \n (minutes since expt start)'));
    ylabel('Rate (imp / sec)');
end

%
% DETERMINE IF A FIGURE WINDOW CURRENTLY EXISTS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function guiFigHand = currentFigOpen(figName);
    guiFigHand = [];
    figList = get(0, 'children');
    for a = 1:length(figList)
        if strcmpi(figName, get(figList(a), 'name'))
            guiFigHand = figList(a);
        end
    end
end

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    