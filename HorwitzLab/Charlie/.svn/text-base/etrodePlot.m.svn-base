function etrodePlot(txtFile)
%
% EXAMPLE: etrodePlot('textFile.txt')
%
% Plots the electrode penetration locations from a tab delimited text
% file. The contrast of each point on the plot indicates the number of
% penetrations.
%
% CAH 08/08


if ~strcmpi(txtFile(end-3:end), '.txt')
    error('not a valid text file')
end

fid = fopen(txtFile);

data = [];
if fid > 0;
    while (1)
        line = fgetl(fid);
        if ~ischar(line)
            break
        end
        if (numel(line) > 0) && ~strcmpi(line(1), '%')
            data = [data;str2num(line)];
        end
    end
else
    fclose(fid);
    error('unable to open file')
end
fclose(fid);

uniquePos = unique([data(:,2), data(:,3)], 'rows');
for a = 1:size(uniquePos, 1)
    nVisits(a) = sum([uniquePos(a,1) == data(:,2)] & [uniquePos(a,2) == data(:,3)]);
end

%normalize the intensities for plotting
normNVisits = 1 - [nVisits ./ max(nVisits)]; %remember that [1 1 1] is white so subtract these vals from 1

%now plot
figure
hold on,
plot([5 -5 -5 5 5], [5 5 -5 -5 5], 'k');
for a = 1:size(uniquePos, 1)
scatter(uniquePos(a,1), uniquePos(a,2), 'MarkerFaceColor', repmat(normNVisits(a), 1, 3), 'MarkerEdgeColor', repmat(normNVisits(a), 1, 3))
end
hold off
xlim([-5.5 5.5])
ylim([-5.5 5.5])
set(gca, 'Xtick', linspace(-5, 5, 11), 'XtickLabel', linspace(-5, 5, 11));
cMap = (repmat(linspace(max(normNVisits), min(normNVisits), 64)', 1, 3));
colormap(cMap);
h = colorbar;
l = get(h, 'YtickLabel');
l = linspace(min(nVisits), max(nVisits), length(l));
set(h, 'YtickLabel', l);
title('Position of Electrode Penetrations')
ylabel('(-)Posterior / (+)Anterior')
xlabel('(-)Medial / (+)Lateral')


%adding a comment and committing from charlie's laptop to check svn






