function sacstats = getSacData(stro, curvethresh)

% function sacstats = getSacData(stro,[curvethresh]);
%
%     This function will run through stro data structure
% calculate a bunch of saccade statistics.  Results will be 
% returned in the sacstats structure.
%
%     The format of  sacstats is as follows: each row in each
% field of the data structure is a saccade-event.  Corresponding
% rows in each field refer to the same saccade event.
%
%    The fields of the sacstats structure is as follows:
%
%	trialnums = which trial this saccade comes from
%	starttimes = start time (relative to trial start)
%	endtimes = end time (relative to trial start)
%	durations = duration of the movement
%	amplitudes = amplitude of the eye movement
%	directions = direction of the eye movement
%   peakv = the peak velocity of the eye movement
%   pathlength = the piecewise linear path length of the movement
%
%    The algorithm for finding saccades is as follows: the horizontal
% and vertical position traces are smoothed and differentiated.
% saccades are identified as times at which the eyes moved faster
% than THRESHOLD deg/sec.  Provided the movement is of longer duration
% than MINSACDUR and is briefer than MAXSACDUR, the movement is
% considered a saccade.  If two saccades so identified occur within
% MINSACREFRACT of each other they are treated as a single saccade.
%
%   If the optional argument "curvethresh" is given, saccades are omitted
% if their pathlength exceeds their amplitude + curvethresh.  The 
% empirically determined value of '0.3' for curvethresh seems to work
% pretty well.
% 

% Constants for finding saccades
THRESHOLD = 10;   % Degrees/sec -- speed criterion for saccade detection
MINSACDUR = 8;    % Minimum saccade duration in ms
MAXSACDUR = 100;   % Maximum saccade duration in ms
MINSACREFRACT = 40;   % Minimum time between saccades
SIGMA = 4;		% ms/stdev

hepidx = strcmp(stro.sum.rasterCells(1,:),'AD11');
vepidx = strcmp(stro.sum.rasterCells(1,:),'AD12');
starttimeidx = strcmp(stro.sum.rasterCells(1,:),'anlgStartTime');

samplerate = stro.sum.analog.storeRates{1};
kernel = normpdf([-3:1000/(samplerate*SIGMA):3],0,1);
kernel = diff(kernel./sum(kernel));

figure; axes; 
for i = 1:size(stro.ras,1)
    cla; hold on;
    h = stro.ras{i,hepidx}*4096/400;
    v = stro.ras{i,vepidx}*4096/400;
    
    % 4096 A/D levels = 10 V
    % 1 degree = 40 A/D levels (according to REX)
    % (4096 levels/10 V) * (1 degree/40 levels) = 4096/400.
    % The key is that both REX and PLEXON use 12 bit A/D boards configured
    % for +/- 5 V.
    
   	e1t = stro.ras{i,starttimeidx};
    neyesamp = size(h,1);
    x = linspace(e1t,neyesamp/samplerate+e1t,neyesamp);

    plot(x,h,'Color','green');
    plot(x,v,'Color','red');
    % Smoothing the traces
    tmpe1h = conv(h,kernel)*samplerate;    % in degrees/sec
    tmpe1v = conv(v,kernel)*samplerate;    % in degrees/sec

    tmpe1h([1:length(kernel)-1]) = [];
    tmpe1v([1:length(kernel)-1]) = [];
    tmpe1h([end-length(kernel)+2:end]) = [];
    tmpe1v([end-length(kernel)+2:end]) = [];
    if (length(kernel)/2 == round(length(kernel)/2))    % Even kernel
        x = conv(x,[.5 .5]);   % Finding the intermediate point
    end
    x([1:floor(length(kernel)/2)]) = [];
    x([end-floor(length(kernel)/2)+1:end]) = [];

    plot(x,tmpe1h,'Color','green')
    plot(x,tmpe1v,'Color','red')

    speed = sqrt(tmpe1h.^2+tmpe1v.^2);
    plot(x, speed);
    % Because of the logic below, sacbeginidx is the last sample below THRESHOLD.
    % sacendidx is the last sample that's still above threhold.
    sacbeginidx = find(diff(speed > THRESHOLD) == 1);
    sacendidx = find(diff(speed > THRESHOLD) == -1);
    sacbegintime = (x(sacbeginidx)+x(sacbeginidx+1))/2;
    sacendtime = (x(sacendidx)+x(sacendidx+1))/2;

    % Checks I need to make:
    % 0) Eliminate starts without stops (and stops without starts)
    % 1) Eliminate events that are too short
    % 2) Eliminate events that are too long
    % 3) Merge two events if they are close to 
    % each other in time and the sum doesn't exceed the maximum duration
    beginidx = [];
    endidx = [];
    hbeginidx = [];
    hendidx = [];
    maxspeeds = [];
    pathlens = [];
    for j = 1:length(sacbeginidx)
        tmpendtimes = sacendtime-sacbegintime(j);
        tmpendtime = min(tmpendtimes(tmpendtimes > 0));
        if (~isempty(tmpendtime))   % Unmatched saccade start
            tmpendidx = find(tmpendtimes == tmpendtime);
            if (tmpendtime > MINSACDUR/1000 & tmpendtime < MAXSACDUR/1000)
                idx1 = sacbeginidx(j);
                idx2 = sacendidx(tmpendidx);
                beginidx = [beginidx, idx1];
                endidx = [endidx, idx2];
                idx3 = idx1 + ceil((length(kernel)-1)/2);
                idx4 = idx2 + ceil((length(kernel)-1)/2);
                hbeginidx = [hbeginidx, idx3];
                hendidx = [hendidx, idx4];
                pathlens = [pathlens, findpathlength(h(idx3:idx4),v(idx3:idx4))];
                maxspeeds = [maxspeeds, max(speed(sacbeginidx(j):sacendidx(tmpendidx)))];
            end
        end
    end
    % Moving the start and end times slightly to account for the fact that the "diff"
    % operator takes values at i and i+1 and puts the result in i.
    sacbegintime = (x(beginidx)+x(beginidx+1))/2;
    sacendtime = (x(endidx)+x(endidx+1))/2;
    sacbeginh = h(hbeginidx);
    sacendh = h(hendidx);
    sacbeginv = v(hbeginidx);
    sacendv = v(hendidx);
   
    % Combining events that are too close in time to be considered separate
    % events
    if (length(sacbegintime) > 1)
        latencies = sacbegintime(2:end) - sacendtime(1:end-1);
        L = latencies < MINSACREFRACT/1000;  % 1 in position i means saccades i and i+1 are too close
        for j = find(L)    % recursive method
            maxspeeds(j+1) = max([maxspeeds(j); maxspeeds(j+1)]);;
            pathlens(j+1) = pathlens(j)+pathlens(j+1)+...
                sqrt((sacendh(j)-sacbeginh(j+1))^2+(sacendv(j)-sacbeginv(j+1))^2);
        end
        sacbegintime(logical([0 L])) = [];
        sacendtime(logical([L 0])) = [];
        sacbeginh(logical([0 L])) = [];
        sacendh(logical([L 0])) = [];
        sacbeginv(logical([0 L])) = [];
        sacendv(logical([L 0])) = [];
        maxspeeds(logical([L 0])) = [];
        pathlens(logical([L 0])) = [];
    end

    if (~isempty(sacbegintime))
        plot(sacbegintime,THRESHOLD,'g*');
        plot(sacendtime,THRESHOLD,'r*');
    end
    
    % Getting rid of curved "saccades" if that's what the use wants
    if (nargin > 1)
        L = logical(pathlens > sqrt((sacendh'-sacbeginh').^2+(sacendv'-sacbeginv').^2) + curvethresh);
        sacbegintime(L) = [];
        sacendtime(L) = [];
        sacbeginh(L) = [];
        sacendh(L) = [];
        sacbeginv(L) = [];
        sacendv(L) = [];
        maxspeeds(L) = [];
        pathlens(L) = [];
    end
    sacstats.trialnums{i} = i*ones(length(sacbegintime),1);
    sacstats.durations{i} = (sacendtime-sacbegintime)';
    sacstats.starttimes{i} = sacbegintime';
    sacstats.endtimes{i} = sacendtime';
    sacstats.amplitudes{i} = sqrt((sacendh'-sacbeginh').^2+(sacendv'-sacbeginv').^2)';
    sacstats.directions{i} = atan2(sacendv'-sacbeginv',sacendh'-sacbeginh')';
    sacstats.peakv{i} = maxspeeds';
    sacstats.pathlengths{i} = pathlens';
    %pause
end

end
%%%%%%%%%%%%%%%%%%%%%%%%
% Support function for finding the piecewise-linear 
% length of a saccade.
%
   function [d] = findpathlength(x,y)
    %  function [d] = findpathlength(x,y)
    %
    %  Computes the length of the (piecewise linear) path
    %  that connects the x and y coordinates passed in as 
    %  inputs.  It's just the sum of the indididual distances
    %  point 1 to point2, point2 to point3, etc.

    if (length(x) ~= length(y))
        error('input arguments must have the same length.');
    end

    diffxsq = diff(x).^2;
    diffysq = diff(y).^2;
    diffs = diffxsq+diffysq;
    d = sum(sqrt(diffs));
end