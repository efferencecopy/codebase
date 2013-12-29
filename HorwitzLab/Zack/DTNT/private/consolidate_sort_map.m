% This function converts a spatial frequency (SF) -> filename map from size Px2
% to Mx2, where P and M are the number of non-unique and unique SFs,
% respectively. Since experimenters may probe SFs in any order, we consolidate
% all like-SF files in the map so the code can recreate the color direction
% search in one swoop.

function sf_map = consolidate_sort_map(sf_map, listed_sfs)
[unique_vals,sorted_pos] = unique(listed_sfs, 'first');
sf_poss = arrayfun(@(x) {find(listed_sfs == x)}, unique_vals);

% append the all like-SF files into the first SF entry of the map
for curr_pos = sf_poss(:)'
    if numel(curr_pos{1}) > 1
        dest = curr_pos{1}(1);
        src = curr_pos{1}(2:end);
        appended_dest = [sf_map(dest,2); sf_map(src,2)];
        sf_map(dest,2) = {vertcat(appended_dest{:})};
    end
end
sf_map = sf_map(sorted_pos,:);
