function out = iafGetAvgMovie(allTrialMovies)

% allTrialMovies should be a cell array of movies. One movie for each
% trial. And each movie should be Height x Width x Nframes


Ntrials = numel(allTrialMovies);
movieSize = cellfun(@size, allTrialMovies, 'uniformoutput', false);
movieSize = unique(cat(1, movieSize{:}), 'rows');
assert(size(movieSize, 1) == 1, 'ERROR, movie size is not consistent across trials')

% calculate the "sum" portion of the mean
runningSum = zeros (movieSize(1), movieSize(2), movieSize(3));
NvalidTrials = 0;
for i_trl = 1:Ntrials
    if any(isnan(allTrialMovies{i_trl}(:)))
        continue % don't add in any nans
    end
    runningSum = runningSum + allTrialMovies{i_trl};
    NvalidTrials = NvalidTrials + 1;
end

% divide by N to finish off the mean
out = runningSum ./ NvalidTrials;