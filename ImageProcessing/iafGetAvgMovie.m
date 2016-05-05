function out = iafGetAvgMovie(allTrialMovies)

% allTrialMovies should be a cell array of movies. One movie for each
% trial. And each movie should be Height x Width x Nframes


movieSize = cellfun(@size, allTrialMovies, 'uniformoutput', false);
movieSize = unique(cat(1, movieSize{:}), 'rows');
assert(size(movieSize, 1) == 1, 'ERROR, movie size is not consistent across trials')

allTrialMovies = cellfun(@(x) x(:), allTrialMovies, 'uniformoutput', false);

allFrames = cat(2, allTrialMovies{:}); % all the data
allFrames(abs(allFrames)==inf) = NaN; % remove the infs
allFrames = nanmean(allFrames, 2); % average across trials, but ignore NaNs
out = reshape(allFrames, movieSize);

