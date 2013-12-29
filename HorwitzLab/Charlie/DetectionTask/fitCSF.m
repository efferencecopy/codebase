function fitParams = fitCSF(colors, sfs, data, order, fitType)
%
% POLYNOMIAL FITS OF SPATIAL CONTRAST SENSITIVITY FUNCTIONS
%
%   EXAMPLE: fitParams = fitCSF(colors, sfs, data, [order], [fitType]);
%
% By default, order is set to 2. fitType defaults to fitting all the data,
% but can be set to fit just the mean sensitivity. fitType can be 'all', or
% 'mean'
%
% CAH

if ~exist('order', 'var')
    order = 2;
end
if ~exist('fitType', 'var')
    fitType = 'all';
end

fitParams.hand = figure;
plotOrder = [1,     -1,       0;...
             0,      0,       1;...
             1,      1,       1;...
             0.1386, -0.1386, -0.9806;...
             0.1386, -0.1386, 0.9806];
plotOrder = plotOrder ./ repmat(sqrt(sum(plotOrder.^2, 2)), 1, 3);
colorOrder = {'r', 'b', 'k', 'g', 'm'};
sensitivity = 1./(data./100); %divide by 100 to put percentage cone contrast b/w 0&1
avgSen = nanmean(sensitivity, 3);
SEM = nanstd(sensitivity, 0, 3)./sqrt(sum(~isnan(data),3));
for a = 1:size(plotOrder, 1);
    idxToPlot = softEq(plotOrder(a,:), colors, 4, 'rows');
    if any(idxToPlot)
        hold on,
        h.csf = [];
        h.csf = errorbar(sfs, avgSen(idxToPlot, :), SEM(idxToPlot, :),  colorOrder{a});
        set(h.csf, 'marker', '.', 'linestyle', 'none')

        % fit the spatial contrast sensitivity functions with a
        % polynomial regression.
        sfsToFit = ~isnan(avgSen(idxToPlot, :)); %some of the sfs might not be represented and will screw up the fitting
        if (sum(sfsToFit) >= 3)
            switch lower(fitType)
                case 'all'
                    sfsIdxs = find(sfsToFit);
                    tmp_data = [];
                    tmp_sfs = [];
                    for sf = 1:length(sfsIdxs);
                        tmp = sensitivity(idxToPlot, sfsIdxs(sf), :);
                        tmp(isnan(tmp)) = [];
                        tmp = permute(tmp, [1,3,2]);
                        tmp_sfs = [tmp_sfs, repmat(sfs(sf), 1, length(tmp))];
                        tmp_data = [tmp_data, tmp];
                    end
                    polyCoeff = polyfit(tmp_sfs, tmp_data, order);
                case 'mean'
                    polyCoeff = polyfit(sfs(sfsToFit), avgSen(idxToPlot, sfsToFit), order);
            end
            x = [min(sfs)*.8 : 0.001 : max(sfs)*1.2];
            model = polyval(polyCoeff,x);
            model(model<0) = nan; %so that neg values don't throw a warning during logplot
            plot(x, model, colorOrder{a}, 'linewidth', 2);

            %compile the fits into a single structure array
            fitParams.colors = plotOrder;
            fitParams.polyCoeff(a, 1:3) = polyCoeff;
        end
        hold off  
    end
end
maxy = max(max(max(sensitivity))) .* 1.2;
miny = min(min(min(sensitivity))) .* 0.95;
set(gca, 'Xscale', 'log', 'Yscale', 'log', 'Xlim', [0.2 3.5462], 'Ylim', [miny, maxy]);
