function cnt_LUM = luminanceContrast(cnt_LMS, Mmtx, bkgnd_rgb, monspect, cnt_RGB)

% load in Vlambda and calculate the LUM due to bkgnd
load('T_vos1978_Y');
Vlambda = T_vos1978_Y;
lum_bkgnd = bkgnd_rgb * monspect' * Vlambda(:);

if ~exist('cnt_RGB', 'var') || isempty(cnt_RGB)
    
    bkgnd_lms = Mmtx * bkgnd_rgb(:);
    lms = (cnt_LMS(:) + 1) .* bkgnd_lms(:);
    rgb_stim = Mmtx \ lms;
    
else 
    
    rgb_stim = (cnt_RGB(:)+1) .* bkgnd_rgb(:);
    
end


lum_stim = rgb_stim' * monspect' * Vlambda(:);
cnt_LUM = (lum_stim - lum_bkgnd) ./ lum_bkgnd;
