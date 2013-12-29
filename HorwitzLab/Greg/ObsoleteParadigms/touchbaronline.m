function p = touchbaronline()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Online analysis stuff
% GDLH 10/31/07
% This code analyzes Touchbar.d data online and will hopefully become a template 
% for future online analysis tools 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global stats;  % 'stats' has to be global because the uicontrol callbacks
                % have to be able to manipulate it.  The convention I will
                % follow is that I will pass stats as an argument to any
                % function that doesn't modify it, but I'll use it as a
                % global in any function that does modify it.

% Some initializations
TouchbarCodes;        % Script that defines a bunch of constants
p = InitPStruct(EOTCD); % Structure of the events from each trial
s = InitPlex();     % Link to Plexon datastream
InitStatsStruct();  % Structure that holds the statistics
[sock, Success] = pnetStart(6665);  % UDP communication with REX
SetUpFig();     % Setting up the figure
stopnow = 0;    % Flag to break out of the endless while loop (hit ESC to set to 1)

% The main loop
while (~stopnow)
    stopnow = CheckForESCKey();
    [n, eventList] = PL_GetTS(s);
    if (n > 0)
       p = ProcessEventList(p, eventList);
    end
    
    if (~p.processnowflag)
        continue;
    end
    % Doing the calculations
    corcdidx = find(p.events == CORRECTCD);
    if (isempty(corcdidx))
        error('Couldn''t find a CORRECTCD');
    end
    stats.corvect = [stats.corvect; p.events(corcdidx+1)-VALOFFSET];
    stats.ntotal = stats.ntotal + 1;
    touch_idx = find(p.events == TOUCHCD);
    fpdim_idx = find(p.events == FPDIMCD);
    release_idx = find(p.events == REWCD);
    if (~isempty(touch_idx) & ~isempty(fpdim_idx) & ~isempty(release_idx))
        touch_t = p.times(touch_idx);
        fpdim_t = p.times(fpdim_idx);
        release_t = p.times(release_idx);
        stats.latencies = [stats.latencies; release_t-fpdim_t];
        stats.nondimdurations = [stats.nondimdurations; fpdim_t-touch_t];
    end
    % Doing the plotting
    if (sum(~isnan(stats.latencies)) > 2)
        subplot(1,3,1);
        LatencyHist(stats);
        subplot(1,3,2);
        LatencyScatter(stats);
        subplot(1,3,3);
        PerformanceVsTime(stats);
        drawnow;
    end
    p = CleanUpEvents(p);
end
PL_Close(s);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Statistics that get updated every trial
function InitStatsStruct()
    global stats;

    stats.corvect = [];  % Binary vector for correct and incorrect responses
    stats.latencies = []; % Vector of release latencies (FPdim to bar release)
    stats.nondimdurations = [];
    stats.ntotal = 0;    % Total number of trials
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up the figure for plotting
function SetUpFig()
    figure(1);
    clf;
    set(gcf,'position',[151 369 915 213]);
    h = uicontrol('style','pushbutton','Callback',@ResetCallback, 'string','RESET');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Function that gets executed when you 
    % click the RESET button
    function ResetCallback(hobj, ed)
      SetUpFig;
      drawnow;
      InitStatsStruct;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting functions below...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A histogram of latencies
function LatencyHist(stats)
    hist(stats.latencies,length(stats.latencies)/5+2);
    xlabel('latency (ms)');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A scatterplot of latency as a function of
% the time from the initial bar touch to FP dim.
function LatencyScatter(stats)
    plot(stats.nondimdurations, stats.latencies,'k.');
    xlabel('nondim duration (ms)');
    ylabel('latency (ms)');
    X = [ones(length(stats.nondimdurations),1) stats.nondimdurations];
    b = regress(stats.latencies, X);
    hold on;
    plot([min(stats.nondimdurations) max(stats.nondimdurations)],...
       b(2)*[min(stats.nondimdurations) max(stats.nondimdurations)]+b(1),'b-');
    textstr = sprintf('intercept = %3.0f\nslope = %3.3f',b(1),b(2));
    xlim = get(gca,'Xlim');
    ylim = get(gca,'Ylim');
    text(xlim(1)+.1*range(xlim), ylim(2)-.2*range(ylim), textstr);
    hold off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adaptively smoothed trace of percent correct
% as a function of trial index.
function PerformanceVsTime(stats)
    cla;
    hold on;
    kernel = ones(floor(length(stats.corvect)/5),1);
    kernel = kernel./sum(kernel);
    tmpcorvect = conv(stats.corvect,kernel);
    tmpcorvect(1:length(kernel)-1) = [];
    tmpcorvect(end-(length(kernel)-1):end) = [];
    plot(tmpcorvect*100);
    ylabel('% correct')
    xlabel('trial');
    set(gca,'xlim',[1 length(tmpcorvect)]);
    title(['% correct: ',num2str(sum(stats.corvect)*100/stats.ntotal,3),' (',num2str(sum(stats.corvect)),'/',num2str(stats.ntotal),')']);
end
