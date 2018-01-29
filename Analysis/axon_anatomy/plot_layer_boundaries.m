function fig_hand = plot_layer_boundaries(layers, y_val, bar_width, plt_clr, fig_hand)

if ~exist('fig_hand', 'var') || isempty(fig_hand)
    fig_hand = figure; hold on,
else
    figure(fig_hand); hold on,
end

if ~exist('plt_clr', 'var') || isempty(plt_clr)
    plt_clr = 'k';
end

% each field in "layers" should be a distribution of values. Each layer
% boundry gets a distribution (eg., L23-L4 boundary), and the distribution
% is composed of the values from the different microscope images.
yy_l1_l23 = cellfun(@mean, layers.L1_L23);
yy_l23_l4 = cellfun(@mean, layers.L23_L4);
yy_l4_l5 = cellfun(@mean, layers.L4_L5);
yy_l5_l6 = cellfun(@mean, layers.L5_L6);
yy_axons = cellfun(@mean, layers.bottom);

% make the stacked bar figure
layer_bounds = [mean(yy_l1_l23), mean(yy_l23_l4), mean(yy_l4_l5), mean(yy_l5_l6), mean(yy_axons)];
y = diff(layer_bounds); % barh plots the cumsum of "Y" not the explicit locations of each boundary
x = [y_val, y_val-bar_width];
hb = barh(x, [y ; nan(size(y))], 0.8, 'stacked');
set(hb, 'FaceAlpha', 0, 'LineWidth', 1, 'EdgeColor', plt_clr);
bl = get(hb, 'BaseLine');
cellfun(@(x) set(x, 'Color', 'none'), bl);


% add the horizontal error bars
bounds_half_sigma = [std(yy_l23_l4), std(yy_l4_l5), std(yy_l5_l6), std(yy_axons)] ./ 2;
l_pos = cumsum(y) - bounds_half_sigma;
r_pos = cumsum(y) + bounds_half_sigma;
err_lines = [l_pos; r_pos];
hold on,
plot(err_lines, ones(size(err_lines)).*y_val, '-', 'color', plt_clr, 'linewidth', 2)

