function smooth_p1_amps = test_vca_make_smooth_p1(fake_p1_amps, N)

initmean = nanmean(fake_p1_amps(1:5));
endmean = nanmean(fake_p1_amps(end-4:end));
tmp = [ones(1,N).*initmean , fake_p1_amps(:)', ones(1,N+2).*endmean];
smooth_p1_amps = nan(size(tmp));
for i_swp = 1:(numel(smooth_p1_amps)-N)
    smooth_p1_amps(i_swp) = nanmean(tmp(i_swp:(i_swp+(N-1))));
end

assert(rem(N, 2)==0, 'ERROR: N saps in box car must be even')
half_width = N./2;
smooth_p1_amps(1:half_width)=[];
smooth_p1_amps(numel(fake_p1_amps)+1:end) = [];