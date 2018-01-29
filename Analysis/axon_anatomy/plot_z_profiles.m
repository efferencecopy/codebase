function plot_z_profiles(dat, all_hvas)

hf = figure;
hf.Units = 'normalized';
hf.Position = [0.55 0.29111 0.58542 0.36];

ax1 = subplot(1,2,1); hold on,
n_z_planes = [];
all_mice = fieldnames(dat);
for i_mouse = 1:numel(all_mice)
    mouse = all_mice{i_mouse};
    for i_hva = 1:numel(all_hvas)
        hva = all_hvas{i_hva};
        if isempty(dat.(mouse).(hva).images)
            continue
        end
        
        for i_img = 1:numel(dat.(mouse).(hva).z_profile)
            yy = dat.(mouse).(hva).z_profile{i_img};
            yy = yy ./ max(yy);
            xx = linspace(-numel(yy)/2, numel(yy)/2, numel(yy));
            plot(xx, yy)
        end
        n_z_planes = cat(1, n_z_planes, dat.(mouse).(hva).n_zplanes{:});
    end
end
ax1.Title.String = 'Z profile above threshold';
ax1.XLabel.String = 'Z-planes from center';
ax1.YLabel.String = 'Normalized intensity';
ax1.YLim = [0.48 1.005];


ax2 = subplot(1,2,2);
histogram(n_z_planes);
ax2.XLabel.String = 'Number of focal planes above threshold';
ax2.YLabel.String = 'Number of images';