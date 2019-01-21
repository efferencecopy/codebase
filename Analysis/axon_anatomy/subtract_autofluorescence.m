function dat = subtract_autofluorescence(dat, all_hvas)


mouse_names = fieldnames(dat);
plt_clr = 'k'; % so that the user is agnostic to area
for i_mouse = 1:numel(mouse_names)
    mouse = mouse_names{i_mouse};
    
    
    
    hf = figure;
    hf.Position = [16         254        1416         468];
    ax1 = subplot(1,3,1); hold on,
    raw_profiles = struct();
    for i_hva = 1:numel(all_hvas)
        
        % initialize the (temporary) structure to collect the raw
        % z-profiles
        raw_profiles.(all_hvas{i_hva}).xx_common = [];
        raw_profiles.(all_hvas{i_hva}).yy_avg = [];
        
        % skip hvas for which there is no data
        has_no_data = isempty(dat.(mouse).(all_hvas{i_hva}).images);
        if has_no_data; continue; end
        
        % compute yy_profile average across different lengthed profiles
        [xx_common, yy_avg] = compute_avg_profile(dat, mouse, all_hvas{i_hva});
        raw_profiles.(all_hvas{i_hva}).xx_common = xx_common;
        raw_profiles.(all_hvas{i_hva}).yy_avg = yy_avg;
        
        plot(xx_common, yy_avg, '-', 'color', plt_clr, 'linewidth', 2)
        
    end
    hf.Name = mouse;
    ylabel('Raw Fluroescence');
    axis tight
    
    % compute the corrected z-profiles
    baselined_profiles = struct();
    for i_hva = 1:numel(all_hvas)
        % initialize the output data structure
        baselined_profiles.(all_hvas{i_hva}).xx_common = [];
        baselined_profiles.(all_hvas{i_hva}).yy_baselined = [];
        
        % skip hvas for which there is no data
        has_no_data = isempty(dat.(mouse).(all_hvas{i_hva}).images);
        if has_no_data; continue; end
        
        % plot the raw image as a reference
        subplot(1,3,2)
        img_raw = dat.(mouse).(all_hvas{i_hva}).images{2};
        image(max(img_raw, [], 3));
        colormap((gray(256)))
        
        % grab the raw data
        xx_hva = raw_profiles.(all_hvas{i_hva}).xx_common;
        yy_hva = raw_profiles.(all_hvas{i_hva}).yy_avg;
        
        ax2 = subplot(1,3,3);
        cla(ax2)
        plot(xx_hva, yy_hva,  '-', 'color', plt_clr, 'linewidth', 2)
        ylabel('Baseline subtracted Fluroescence');
        axis tight
        
        % ask the user to add a line
        [x_bkgnd, y_bkgnd] = ginput(2);
        coeff = [x_bkgnd(:), ones(numel(x_bkgnd),1)] \ y_bkgnd(:);
        plot(xx_common, coeff(1).*xx_common+coeff(2), '--k');
        
        % computed the baselined profile
        yy_baseline = coeff(1) .* xx_hva + coeff(2);
        subtracted_profile = yy_hva - yy_baseline;
        
        % omit stuff that's negative
        l_neg = subtracted_profile < 0;
        subtracted_profile(l_neg) = 0;
        
        % store the new values
        baselined_profiles.(all_hvas{i_hva}).xx_common = xx_hva;
        baselined_profiles.(all_hvas{i_hva}).yy_baselined = subtracted_profile;
        
    end
    
    % include these data in the dat struct for persistance
    dat.(mouse).avg_profiles = raw_profiles;
    dat.(mouse).baselined_profiles = baselined_profiles;
    dat.(mouse).baseline_coeffs = coeff;
    
end