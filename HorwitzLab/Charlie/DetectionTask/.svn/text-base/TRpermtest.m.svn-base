function [p, permstats, Fobs] = TRpermtest()
    % works on the log(TRs)

    % prepare the batch data for further analysis
    global cardVsIntBatchPath blpBatchPath
    preprocessDTbatchData


    % identify the neuons that meed the inclusion criteria
    l_valid = ~(commonExclusions | (out.errors(:, neuroThresh1Ind) & out.errors(:, neuroThresh2Ind)));
    validExpts = find(l_valid);
    
    %cycle through the valid neurons and assign the TRs to their native color
    %directions
    colors = {'L-M', 'S-iso', 'SwM', 'SwL'}; %defines the mapping b/w colum and color dir
    dataMap = nan(sum(l_valid), 4);
    for a = 1:size(dataMap,1);
        ex = validExpts(a);
        exColors = out.dat(ex).expt.standColors;
        tmpTRs = log(rawTRs(ex,:));
        tmpTRs(isnan(tmpTRs)) = -inf;

        %assign the TRs
        if any(ismember(sign(exColors), [1 -1 0], 'rows'))
            dataMap(a,1) = tmpTRs(ismember(sign(exColors), [1 -1 0], 'rows'));
        elseif any(ismember(sign(exColors), [0 0 1], 'rows'))
            dataMap(a,2) = tmpTRs(ismember(sign(exColors), [0 0 1], 'rows'));
        end

        if any(ismember(sign(exColors), [1 -1 -1], 'rows'))
            dataMap(a,3) = tmpTRs(ismember(sign(exColors), [1 -1 -1], 'rows'));
        elseif any(ismember(sign(exColors), [1 -1 1], 'rows'))
            dataMap(a,4) = tmpTRs(ismember(sign(exColors), [1 -1 1], 'rows'));
        end

    end

    % calculate the observed F-statistic, but convert the -infs (used as
    % place holders) into nans (which will be ignored by the test)
    tmpData = dataMap;
    tmpData(tmpData<0) = nan;
    K = size(tmpData,2);
    grand_N = sum(~isnan(tmpData(:)));
    Fobs = fstat(tmpData, K, grand_N);


    % set things up for the permutations
    dataMap = ~isnan(dataMap); % a map of which conditions were tested
    filtRawTRs = log(rawTRs(l_valid,:));
    nPerms = 5000;
    permstats = nan(nPerms, 1);
    for a = 1:nPerms

        %shuffle the data
        tmpTRs = filtRawTRs(randperm(sum(l_valid))',:);
        tmpTRs(isnan(tmpTRs)) = -inf; %necessary for book keeping
        tmpData = nan(size(dataMap));
        for i = 1:size(tmpData,1)
            rowIdx = randperm(2); %assign the "shuffeled" TRs to randomized positions in the real data map
            tmpData(i, find(dataMap(i,:),1, 'first')) = tmpTRs(i,rowIdx(1));
            tmpData(i, find(dataMap(i,:),1, 'last')) = tmpTRs(i,rowIdx(2));
        end

        % calculate the F-statistic, but convert the -infs back to nans
        tmpData(tmpData<0) = NaN;
        permstats(a) = fstat(tmpData, K, grand_N);

    end
    
    %plot some stuff
    dx = 0.01;
    x = 0:dx:6;
    y = fpdf(x, K-1, grand_N-K);

    edges = linspace(0, max(permstats).*1.05, 30);
    counts = histc(permstats, edges);
    counts = counts ./ nPerms;

    figure, hold on,
    bar(edges, counts, 'type', 'histc')
    plot(x,y, 'b')
    plot([Fobs, Fobs], [0, max(counts)], 'r')

    p = sum(permstats>Fobs) ./ numel(permstats);
end

function F = fstat(tmpData, K, grand_N)
    group_Ns = sum(~isnan(tmpData), 1);
    grandMean = nanmean(tmpData(:));

    %determine the between groups mean
    groupMeans = nanmean(tmpData,1);
    bwGroups = sum(group_Ns .* ((groupMeans - grandMean).^2)) ./ (K-1);

    %determine the within groups mean
    withinGroups = (bsxfun(@minus,tmpData, groupMeans).^2); %squared diffs from group means
    withinGroups(isnan(withinGroups)) = 0; %changing the NaNs to Zeros won't affect the sums
    withinGroups = sum(withinGroups,1); %w/in groups SSE
    withinGroups = sum(withinGroups); % sum of w/in groups SSE
    withinGroups = withinGroups ./ (grand_N - K);

    F = bwGroups ./ withinGroups;
end


