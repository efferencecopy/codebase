function [r,p] = prettycorr(X,varnames,type)
% function [r,p] = prettycorr(X,<varnames>,<type>)
%
% Makes an interactive image of a correlation matrix.
% Clock on a cell to see an small scatter plot.
% GDLH 9/2/13

k = size(X,2);
if (nargin < 3)
    type = 'Pearson';
end
if (~any(strcmp(type,{'Pearson','Kendall','Spearman'})))
    error('type (3rd arg) must be "Pearson" (the default), "Kendall", or "Spearman"');
end
if (nargin == 1 )
    varnames = cell(k,1);
end

if (length(varnames) ~= k)
    error('mismatch in # of variables and  labels');
end
[r,p] = corr(X,'type',type);

Makeplot(r,p);

    function Makeplot(r,p)
        figure; axes; hold on;
        him = image((r+1)*32);  colormap(gray(64));
        axis ij;
        for i = 1:k
            for j = 1:k
                if (i == j)
                    continue
                end
                
                h = text(i,j-.2,['r = ',num2str(r(i,j),2)],'HorizontalAlignment','center','ButtonDownFcn',@dispscatter);
                if (r(i,j)<0)
                    set(h,'color',[1 1 1]);
                end
                if (p(i,j) < 0.05)
                    set(h,'color',[1 0 0]);
                end
                
                h = text(i,j+.2,['p = ',num2str(p(i,j),2)],'HorizontalAlignment','center','ButtonDownFcn',@dispscatter);
                if (r(i,j)<0)
                    set(h,'color',[1 1 1]);
                end
                if (p(i,j) < 0.05)
                    set(h,'color',[1 0 0]);
                end
            end
        end
        set(gca,'XTick',1:k,'XTickLabel',varnames);
        set(gca,'YTick',1:k,'YTickLabel',varnames);
        axis tight
        set(him,'ButtonDownFcn',@dispscatter);
    end

    % Support function for prettycorr
    function dispscatter(h,ev)
        whichpt = get(gca,'CurrentPoint');
        whichpt = round(whichpt(1,[1 2]));
        whichpt = min([whichpt; k k]);
        figure('position',[ 78   517   239   211]);
        plot(X(:,whichpt(1)),X(:,whichpt(2)),'k.');
        xlabel(varnames(whichpt(1)));
        ylabel(varnames(whichpt(2)));
        lsline
    end
end