% NOTES
%
% ### TO DO
%
% Export the PSCs for each stimulus location (ordered in the same fashion
% as the file names) and the location matrix itself.
%
% Save the params stuff in the MDB. Do some error checking to make certian
% that the stuff in the params field lands in the correct entry on the mdb.
% Also, I shouldn't save the abfobj's becuase they are redundant with the
% raw data, so I should flag this and omit it from the mdb. Otherwise the
% mdb will get huge.
%
% Add the analysis tags to the searchable text fields? (this would be
% overwritten when new files are added or when updates are folded into the
% mdb.
%
% Make a text file. It should document each HVA and neurons that are
% located there. Maybe make a different text file for L2/3 PY cells, L2/3
% FS cells.
%
% For now, the text file should include
% * mouse name
% * cell number
% * stimulation location
% * HVA location
%
%
%
% ### ANALYSIS
%
% Look at norm PSC amp train for each TF. All cells plotted on the same set
% of axes (one axis for each TF, one figure for each HVA)
%
% Look at PPratio P1:P2 as function of TF. All cells get plotted on the
% same set of axes (one axis for each HVA).
%
% Head to head comparison for HVAs: plot the average PPratio vs TF for each
% HVA on the same set of axes.
%
% Some how compare the effect of different stimulus location so that we can
% figure out where to stimulate (instead of stimulating at the soma)
%
% Plot the PN:P1 ratio for each TF, fit an exponential to the data sets
% where there is only monotonic depression. Compare tau's across HVAs
%
% What happens to the magnitued and frequency of spontaneous events
% following the pulse train protocol? For facilitating synapses, is there
% less spontaneous activity? [Although the spontaneously active synapses
% might not be the ones that undergo evoked release]





%% 3D PLOT OF TF vs PP ratio vs STIM LOCATION

fin

% specify which neurons you want to analyze
attribute.area = '.';
attribute.type = 'PY';
attribute.layer = '2/3';
attribute.multiple = 1;

% other parameters
pulsenumber = 2; % which pulse to compare with the first for PP_ratio
[~, pltclr] = hvaPlotColor(attribute.area);
newfig = 1;

% locate the neurons in the mouseDB
[mdbidx, cellnum] = cellListFromXL('Cell_Library', 1, attribute);


% bin cortical distance and TFs. Then setup the data matrix, which will be
% Ndists by Ntfs by Ncells. 'Dist' defines bin centers. So if the distance
% from the objective to cell body is 75um than this data point will fall
% in the first bin (the 0 to 100 bin).
dist = 50:50:600;
binedges = dist-50;
TFs = [5,10,20,40];
popdat_pp = nan(numel(TFs), numel(dist), numel(mdbidx));
popdat_sum = nan(numel(TFs), numel(dist), numel(mdbidx));
pltTF = repmat(TFs(:), 1, size(popdat_pp,2));
pltDist = repmat(binedges, size(popdat_pp,1), 1);


% iterate over the data files and compile the PPratios according to the TF
% and the distance of the objective to the soma.
mdb = initMouseDB(0, true);
if newfig; 
    f_pp = figure; hold on,
    set(f_pp, 'name', 'Paired Pulse Ratio')
    f_sum = figure; hold on,
    set(f_sum, 'name', 'Pulse Sum Value')
end
for neuron = 1:numel(mdbidx);
    
    % grab the data
    celldat = mdb.mice{mdbidx(neuron)}.popAnly{cellnum(neuron)};
    
    % each neuron has several files. Iterate through these and assign the
    % PPratio to the appropriate row and col idx for the popdat matrix
    for fl = 1:numel(celldat.files)
        
        % don't analyze if the TF is < 5;
        if celldat.TF(fl) < 5-eps
            fprintf('TF = %.3f discarded\n', celldat.TF(fl))
            continue
        end
        
        % deal with the PP ratio analysis
        rowIdx = softEq(TFs, celldat.TF(fl), 5);
        assert(~isempty(rowIdx), 'ERROR: could not find TF');
        [~, colIdx] = min(abs(dist - norm(celldat.stimLoc(fl,:))));       
        PP_ratio = celldat.pscamp{fl}(pulsenumber) ./ celldat.pscamp{fl}(1);
        popdat_pp(rowIdx, colIdx, neuron) = PP_ratio;
        
        % now calculate the kooky thing across pulses
        normpsc = celldat.pscamp{fl} ./ celldat.pscamp{fl}(1);
        normpsc = normpsc(1:5); % just take the first 5
        popdat_sum(rowIdx, colIdx, neuron) = sum(normpsc)./5;
        
        
    end
    
    % plot the raw PP ratio data
    figure(f_pp)
    tmp = popdat_pp(:, :, neuron);
    plot3(pltTF(:), pltDist(:), tmp(:), 'o', 'markeredgecolor', pltclr)
    
    % plot the raw PP_sum data
    figure(f_sum)
    tmp = popdat_sum(:, :, neuron);
    plot3(pltTF(:), pltDist(:), tmp(:), 'o', 'markeredgecolor', pltclr)
    
    
end


% now plot the average across cells, in a different color, and connect the
% lines.
avg_pp = exp(nanmean(log(popdat_pp),3));
avg_sum = exp(nanmean(log(popdat_sum),3));
figure(f_pp)

% connect things acquired at the same eccentricity
for a = 1:numel(binedges)
    l_valid = ~isnan(avg_pp(:,a));
    
    figure(f_pp)
    plot3(TFs(l_valid)', repmat(binedges(a), sum(l_valid), 1), avg_pp(l_valid, a)', [pltclr,'-'])

    figure(f_sum)
    plot3(TFs(l_valid)', repmat(binedges(a), sum(l_valid), 1), avg_sum(l_valid, a)', [pltclr,'-'])
end

% connect things acquired at the same TF
for a = 1:numel(TFs)    
    l_valid = ~isnan(avg_pp(a,:));
    
    figure(f_pp)
    plot3(repmat(TFs(a), 1, sum(l_valid)), binedges(l_valid), avg_pp(a, l_valid)', [pltclr,'-'])
    
    figure(f_sum)
    plot3(repmat(TFs(a), 1, sum(l_valid)), binedges(l_valid), avg_sum(a, l_valid)', [pltclr,'-'])
        
end



% tidy up the plots
figure(f_pp)
set(gca,'view', [-74    14], 'zscale', 'log', 'xscale', 'log', 'xtick', TFs)
axis tight
xlabel('Temporal Frequency')
ylabel('Distance')
zlabel('Paired Pulse Ratio')

figure(f_sum)
set(gca,'view', [-74    14], 'zscale', 'log', 'xscale', 'log', 'xtick', TFs)
axis tight
xlabel('Temporal Frequency')
ylabel('Distance')
zlabel('Paired Pulse Ratio')




%% PP VS. DISTANCE FOR A SINGLE TF

fin


% specify which neurons you want to analyze
attribute.area = '.';
attribute.type = 'PY';
attribute.layer = '2/3';
attribute.multiple = 1;

% other parameters
TF = 20;
pulsenumber = 2; % which pulse to compare with the first for PP_ratio

% locate the neurons in the mouseDB
[mdbidx, cellnum] = cellListFromXL('Cell_Library', 1, attribute);

% load the mdb
mdb = initMouseDB(0, true);

N = numel(cellnum);
PPratio = {};
stimLoc = {};
for a = 1:N
    
   % grab the data
    celldat = mdb.mice{mdbidx(a)}.popAnly{cellnum(a)}; 
    
    % look for files were the desired TF was actually tested
    l_TF = softEq(TF, celldat.TF, 5);
    if ~any(l_TF); continue; end
    
    %store the PPratio
    pnum = repmat({pulsenumber}, size(celldat.TF));
    tmp_ppratio = cellfun(@(x,y)  x(y)./x(1), celldat.pscamp, pnum);
    
    % store the stimulus position
    tmp_dist = sqrt(sum(celldat.stimLoc.^2, 2));
    
    
    % weed out the non-intended TFs,average across multiple instances of a
    % particular distance
    tmp_dist = tmp_dist(l_TF);
    tmp_ppratio = tmp_ppratio(l_TF);
    unique_dist = unique(tmp_dist);
    unique_ppratio = [];
    for i = 1:numel(unique_dist);
        l_dist = tmp_dist == unique_dist(i);
        unique_ppratio(i) = mean(tmp_ppratio(l_dist));
    end   
    
    % package the unique distances (which should now be in ascending order
    % of stimulus location)
    stimLoc{a} = unique_dist;
    if isempty(stimLoc{a});
        keyboard
    end
    PPratio{a} = unique_ppratio;
    
    
end

% empty cells means that there was no data at that TF/dist combo. weed
% these out
l_empty = cellfun(@isempty, stimLoc);
stimLoc(l_empty) = [];
PPratio(l_empty) = [];
N = numel(PPratio);


% make a figure of all the data
f = figure; hold on,
set(gcf, 'position', [285     5   787   800])
for a = 1:N
    plot(stimLoc{a}, PPratio{a}, '-k.')
end
hold off

% plot the trend line (bined by distance)
start = 200;
delta = 100;
maxDist = max(cellfun(@max,stimLoc));
maxDist = 100 .* ceil(maxDist./100);
binRedge = start:delta:maxDist
binLedge = binRedge-delta;
binLedge(1) = 0
binned_ppratio = nan(N, size(binLedge,2));
for a = 1:N
    for i = 1:numel(PPratio{a})
        % which bin does the PPratio live in?
        bin_idx = (binLedge<=stimLoc{a}(i)) & (binRedge > stimLoc{a}(i));
        binned_ppratio(a, bin_idx) = PPratio{a}(i);
    end
end

binCenter = [binLedge(1), binLedge(2:end)+(delta/2)]

figure(f), hold on,
avg =  nanmean(binned_ppratio, 1);
plot(binCenter, avg, '-bo', 'linewidth', 4, 'markerfacecolor', 'b', 'markersize', 10);
sem = nanstd(binned_ppratio, 1) ./ sqrt(sum(~isnan(binned_ppratio), 1));
plot([binCenter; binCenter], [avg-sem; avg+sem], '-b', 'linewidth', 4)
set(gca, 'yscale', 'log', 'ytick', [0.01, [.2:.4:1.6]], 'yticklabel', [0.01, [.2:.2:1.2]])
axis tight
xlim([binCenter(1)-40, binCenter(end)+30])

datPts = sum(~isnan(binned_ppratio), 1);
minY = min(cellfun(@min, PPratio)).*0.8;
maxY = max(cellfun(@max, PPratio)).*1.05;
for a = 1:numel(datPts)
    t = text(binCenter(a)-20, minY, ['(',num2str(datPts(a)),')']);
    set(t, 'Fontsize', 25)
end
ylim([minY.*0.75 maxY])

% add the labels
set(gca, 'fontsize', 35, 'linewidth', 2)
xlabel('LED location (um from soma)')
ylabel(sprintf('P%d to P1 ratio', pulsenumber))
title(sprintf('PPratio vs. Dist for %d Hz', TF))



%% PP VS.TF FOR A SINGLE STIMULUS LOCATION



% specify which neurons you want to analyze
attribute.area = 'AL';
attribute.type = 'PY';
attribute.layer = '2/3';
attribute.multiple = 1;

% other parameters
distCrit = 350;
pulsenumber = 2; % which pulse to compare with the first for PP_ratio

% locate the neurons in the mouseDB
[mdbidx, cellnum] = cellListFromXL('Cell_Library', 1, attribute);

% load the mdb
mdb = initMouseDB(0, true);

N = numel(cellnum);
PPratio = {};
TFs = {};
for a = 1:N
    
    % grap the data
    celldat = mdb.mice{mdbidx(a)}.popAnly{cellnum(a)};
    
    % figure out what distances were tested. Grab the indicies to the first
    % location that is greater than the distance criterion
    tmp_dist = sqrt(sum(celldat.stimLoc.^2,2));
    uniqueDists = unique(tmp_dist);
    if ~any(uniqueDists >= distCrit)
        continue
    end
    firstDistIdx = find(uniqueDists > distCrit, 1, 'first');
    l_gtCrit = tmp_dist == uniqueDists(firstDistIdx);
    
    % determine the PPratio, and store it along with the TF. Do this in tmp
    % arrays. Sort the tmp arrays by TF, and then store in a perminant
    % array
    TFs{a} = celldat.TF(l_gtCrit);
    pnum = repmat({pulsenumber}, size(TFs{a}));
    PPratio{a} = cellfun(@(x,y) x(y)./x(1), celldat.pscamp(l_gtCrit), pnum);
    
end


% figure out the unique TFs tested across each file. Reorganize the data
% into a big matrix of [Ncells x UniqueTFs]
bigTFvector = cat(2, TFs{:});
uniqueTFs = unique(round(bigTFvector));
uniqueTFs = [5, 20, 40];
N = numel(TFs);
data = nan(N, numel(uniqueTFs));
for a = 1:N
    for i = 1:numel(TFs{a})
        colIdx = uniqueTFs == round(TFs{a}(i));
        data(a, colIdx) = PPratio{a}(i);
    end
end


% Plot the PPratios as a function of TF. Once for all the data, and again
% for the average +/- SEM
avg = nanmean(log10(data),1);
N = sum(~isnan(data),1);
sem = nanstd(log10(data),[],1)./sqrt(N);
if ~exist('fhand', 'var')
    fhand = figure; hold on,
    set(gca, 'fontsize', 35)
    xlabel('Temporal Frequency')
    ylabel('Log(PP Ratio)')
else
    figure(fhand); hold on,
    set(gca, 'fontsize', 35)
    xlabel('Temporal Frequency')
    ylabel('Log(PP Ratio)')
end

[clr_raw, clr_avg] = hvaPlotColor(attribute.area);
plot(uniqueTFs', log10(data'), '-', 'color', clr_raw)
plot(uniqueTFs, avg, '-o','color', clr_avg,  'linewidth', 4, 'markerfacecolor', clr_avg, 'markersize', 10)
plot([uniqueTFs; uniqueTFs], [avg+sem; avg-sem], '-','color', clr_avg,  'linewidth', 4)
set(gca, 'xscale', 'log', 'linewidth', 2, 'fontsize', 35)
set(gca, 'ytick', log10([0.66 1 1.5]), 'ytickLabel', [0.66 1 1.5])
xlim([4 45])


%% PPratio VS. pulse TIME for all TFs.

% specify which neurons you want to analyze
attribute.area = 'AL';
attribute.type = 'PY';
attribute.layer = '2/3';

% other parameters
distCrit = 350;
defaultTFs = [40];

% locate the neurons in the mouseDB
[mdbidx, cellnum] = cellListFromXL('Cell_Library', 1, attribute);

% load the mdb
mdb = initMouseDB(0, true);


N = numel(cellnum);
PPratio = nan(numel(defaultTFs), 4, N);
for a = 1:N
    
    % grap the data
    celldat = mdb.mice{mdbidx(a)}.popAnly{cellnum(a)};
    
    % figure out what distances were tested. Grab the indicies to the first
    % location that is greater than the distance criterion
    tmp_dist = sqrt(sum(celldat.stimLoc.^2,2));
    uniqueDists = unique(tmp_dist);
    if ~any(uniqueDists >= distCrit)
        continue
    end
    firstDistIdx = find(uniqueDists > distCrit, 1, 'first');
    l_gtCrit = tmp_dist == uniqueDists(firstDistIdx);
    
    % determine the PPratio, and store it along with the TF. Do this in tmp
    % arrays. Sort the tmp arrays by TF, and then store in a perminant
    % array
    TFs = celldat.TF(l_gtCrit);
    tmp_ppratio = cellfun(@(x) x./x(1), celldat.pscamp(l_gtCrit), 'uniformoutput', false);
    
    % oder the PPratios by TF in the big array
    for i = 1:numel(TFs)
        rowIdx = round(TFs(i)) == defaultTFs;
        if ~any(rowIdx)
            continue
        end
        
        PPratio(rowIdx, :, a) = tmp_ppratio{i}(1:size(PPratio,2));
    end
    
end


% plot a surface of the PPratio
if ~exist('fhand', 'var')
    fhand = figure; hold on,
    set(gca, 'fontsize', 35)
%     ylabel('Temporal Frequency')
%     xlabel('Pulse Number')
else
    figure(fhand); hold on,
    set(gca, 'fontsize', 35)
%     ylabel('Temporal Frequency')
%     xlabel('Pulse Number')
end


% Things to plot when more than one default TF is defined
if numel(defaultTFs)>1
    [X,Y] = meshgrid(1:5,1:numel(defaultTFs));
    avg = nanmean(PPratio, 3);
    
    pHand = surf(X, Y, avg);
    [clr_raw, clr_avg] = hvaPlotColor(attribute.area);
    set(pHand, 'Facecolor', clr_avg, 'FaceAlpha', .2, 'EdgeColor', clr_avg, 'linewidth', 2)
    set(gca, 'YTick', [1:numel(defaultTFs)], 'YTickLabel', defaultTFs, 'view', [38 10])
    set(gca, 'linewidth', 2)
end

% Things to plot if only one default TF is defined
if numel(defaultTFs)==1
    avg_log  = nanmean(log10(PPratio), 3);
    sem_log = nanstd(log10(PPratio), [], 3)./sqrt(sum(~isnan(PPratio), 3));
    dT = 1./defaultTFs.*1000; % dt in msec
    pulseTimes = 0:dT:(numel(avg_log)-1).*dT;
    
    [clr_raw, clr_avg] = hvaPlotColor(attribute.area);
    plot(pulseTimes, avg_log, '-o','color', clr_avg,  'linewidth', 4, 'markerfacecolor', clr_avg, 'markersize', 10)
    plot([pulseTimes; pulseTimes], [avg_log-sem_log; avg_log+sem_log], '-','color', clr_avg,  'linewidth', 4)
    set(gca, 'ytick', log10([0.25 0.5 1 2]), 'ytickLabel', [0.25 0.5 1 2])
    ylabel('Log(PP Ratio)')
    xlabel('Pulse Time')
    axis tight
    xlim([-10 pulseTimes(end)+20])
end

%% PPratio VS. PULSE MAGNITUDE AT THE SOMA

% specify which neurons you want to analyze
attribute.area = '.';
attribute.type = 'PY';
attribute.layer = '2/3';
%attribute.multiple = 1;

% a few other paramters
pulsenumber = 2;

% locate the neurons in the mouseDB
[mdbidx, cellnum] = cellListFromXL('Cell_Library', 1, attribute);

% load the mdb
mdb = initMouseDB(0, true);

PPratio = [];
pulseMag = [];
for a = 1:numel(mdbidx)
    
    % grab the data
    celldat = mdb.mice{mdbidx(a)}.popAnly{cellnum(a)};
    
    % only consider things at the soma
    tmpDist = sqrt(sum(celldat.stimLoc.^2, 2));
    l_soma = find(round(tmpDist) == 0, 1, 'first');
    
    %only consider 20 Hz.
    if round(celldat.TF(l_soma)) ~= 20;
        continue
    end
    
    % store the data
    PPratio(end+1) = celldat.pscamp{l_soma}(pulsenumber) ./ celldat.pscamp{l_soma}(1);
    pulseMag(end+1) = celldat.pscamp{l_soma}(1);
    
    
end

figure
plot(pulseMag, PPratio, 'ko')
xlabel('Pulse Magnitude')
ylabel('P2 to P1 ratio')



