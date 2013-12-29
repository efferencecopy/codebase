% A SCRIPT THAT PLOTS DETECTION THRESHOLDS OVER TIME. THIS IS IMPORTANT IN
% DETERMINING WHEN AN OBSERVER'S PSYCHOMETRIC PERFORMANCE ASYMPTOTES.


%unpack the .nex files
TOLERANCE = 0.5; %perf must be 0.82 +/- this value to be included

fNames = fnamesFromTxt2('N:\NexFiles\Charlie\Apollo\text files\questTraining.txt');

for a = 1:length(fNames)
    DT = dtobj(fNames{a}{1});
    [d(a).thresh, d(a).colors, d(a).sfs] = DTunpack(DT, 'mode', 10, TOLERANCE);
    close(gcf)
    
    %Standardize the representation of colors by converting to unit vecs
    norms = sqrt(sum(d(a).colors.^2, 2));
    d(a).colors = d(a).colors ./ repmat(norms, 1, 3);
    
    %transform the threshold data into decimal percents (i.e., between 0&1)
    d(a).thresh = d(a).thresh ./ 100;
end

%aggregate the data into a matrix of dimensions nColors x nSfs x nExpts.
%Start by setting up the index matricies "sfs" and "colors". Next, add the
%data to a giant matrix.
colors = [];
sfs = [];
for expt = 1:length(d)
    for s = 1:length(d(expt).sfs);
       if find(sfs == d(expt).sfs(s),1)
           continue
       else
           sfs(end+1) = d(expt).sfs(s);
       end
    end
end

for expt = 1:length(d)
    for c = 1:size(d(expt).colors,1);
        l_color = softEq(d(expt).colors(c,:), colors, 5, 'rows');
        if any(l_color)
            continue
        else
            colors(end+1,:) = d(expt).colors(c,:);
        end
    end
end

%Here's the meat of it. Now I'll add the data.
data = nan(size(colors, 1), length(sfs), length(d));
for expt = 1:length(d)
    for c = 1:size(d(expt).colors, 1)
        rowInd = softEq(d(expt).colors(c,:), colors, 5, 'rows');
        for s = 1:length(d(expt).sfs)
            colInd = [sfs == d(expt).sfs(s)];
            data(rowInd, colInd, expt) = d(expt).thresh(c, s);
        end
    end
end


%plotting.
colSigns = sign(colors);
lineSty = {'k.-', 'k.:', 'k.--', 'ks-', 'ko:', 'kv-', 'kd:', 'k+-'};

for c = 1:size(colSigns,1);
    figure, hold on,
    for a = 1:length(sfs)
        colorIdx = softEq(colSigns(c,:), colSigns, [], 'rows');
        sensitivities = 1./data(colorIdx, a, :);
        sensitivities = permute(sensitivities, [1, 3, 2]);
        plot(sensitivities(~isnan(sensitivities)), lineSty{a})
    end
    title(sprintf('Color: [%s]', num2str(colSigns(c,:))));
    legend(num2str(sfs(:)))
    ylabel('Sensitivity')
end


