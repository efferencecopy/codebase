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





%% TF vs PP ratio vs STIM LOCATION

fin

% specify which neurons you want to analyze
%attribute.area = 'PM';
attribute.type = 'PY';
attribute.layer = '2/3';
%attribute.multiple = 1;

% other parameters
pulsenumber = 4; % which pulse to compare with the first for PP_ratio
pltclr = 'k';
newfig = 1;

% locate the neurons in the mouseDB
[mdbidx, cellnum] = cellListFromXL('Cell_Library', 1, attribute);


% bin cortical distance and TFs. Then setup the data matrix, which will be
% Ndists by Ntfs by Ncells. 'Dist' defines bin centers. So if the distance
% from the objective to cell body is 75um than this data point will fall
% in the first bin (the 0 to 100 bin).
dist = 50:100:650;
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








