%% open a data file

%locate the batch file
if ispc
    cd ('N:\NexFiles\Charlie')
elseif ismac
    cd('/Users/charliehass/LabStuff/NexFiles/Charlie')
end
[fname, fpath] = uigetfile('*.mat', 'Pick a .mat file');
load([fpath, fname]);

%compile some basic things
[l_SwM, l_SwL, l_S, l_LvM, l_other] = deal(logical(zeros(length(out.dat),1)));
for a = 1:length(out.dat)
    
    %the preferred color in the isolum plane
    l_isolum = (out.dat(a).gratings.color.colors * [1;1;0]) == 0;
    isoColors = out.dat(a).gratings.color.colors(l_isolum, :);
    isoLumMax = max(out.dat(a).gratings.color.colresp(l_isolum, 1));
    l_prefColor = out.dat(a).gratings.color.colresp(l_isolum, 1) == isoLumMax;
    if sum(l_prefColor) == 1
        prefIsolum(a,:) = isoColors(l_prefColor, :);
    else
        prefIsolum(a,:) = nan(1,3);
    end
    
    %color sensitivity index
    lumIdx = ismember(sign(out.dat(a).gratings.color.colors), [1 1 0], 'rows');
    lumRate = out.dat(a).gratings.color.colresp(lumIdx, 1);
    rawCSIs(a,1) = isoLumMax ./ lumRate;
    
    %modulation ratio
    rawModRatio(a,1) = out.dat(a).gratings.modulationratio;
    
    %spatial frequency
    rawSF(a,1) = out.dat(a).gratings.sf.prefSF;
    
    %determine which pref color goes with which neuron
    if all(sign(prefIsolum(a,:)) == [0 0 1])
        l_S(a) = 1;
    elseif all(sign(prefIsolum(a,:)) == [1 -1 0])
        l_LvM(a) = 1;
    elseif all(sign(prefIsolum(a,:)) == [1 -1 1])
        l_SwL(a) = 1;
    elseif all(sign(prefIsolum(a,:)) == [1 -1 -1])
        l_SwM(a) = 1;
    else
        l_other(a) = 1;
    end
end


%%  SUMMARY OF PREFERED COLORS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
lvmAxis = [1 -1 0] ./ norm([1 -1 0]);
SAxis = [0 0 .2];
x = prefIsolum * lvmAxis';
y = prefIsolum * SAxis';
thetas = atan(y./x);
bins = linspace(-pi/4,pi/2,4);
rose([thetas;thetas+pi], [bins(:); bins(:)+pi]);
set(get(gca, 'children'), 'linewidth', 2)
title('Prefered Color')
nTies = sum(isnan(prefIsolum(:,1)));
xlabel(sprintf('N = %d expts \n %d tie(s)', size(prefIsolum,1), nTies));



%% cone weights ?!?!? although I don't really know how the prefcolor is calculated
rawWeights = [];
for a = 1:length(out.dat)
    tmp = out.dat(a).gratings.color.prefcolor;
    tmp = tmp ./ sum(abs(tmp)); %normalize the weights
    %make sure that the point all fall in the upper part of the diamond
    if sum(sign(tmp(1:2))) == 0 %i.e., L-M or -L+M
        if tmp(1) > 0;
            tmp = tmp.*-1;
        end
    elseif abs(sum(sign(tmp(1:2)))) == 2 %i.e., L+M or -L-M
        if tmp(1) < 0;
            tmp = tmp.*-1;
        end
    end
    rawWeights(a,:) = tmp; 
end

%plot all the data
negS = rawWeights(:,3) < 0;
figure, hold on,
plot(rawWeights(~negS,1), rawWeights(~negS,2), 'bo')
plot(rawWeights(negS,1), rawWeights(negS,2), 'bo', 'markerfacecolor', 'b')
plot([0, 1], [1, 0], 'k')
plot([0, -1], [1, 0], 'k')
plot([0 0], [1 0], 'k')
plot([1 -1], [0 0], 'k')
plot([-.5 .5], [.5 .5], 'k+', 'markersize', 15)
axis equal
title('Normalized Cone Weights')
hold off



%now according to pref isolum color
figure, hold on,
plot(rawWeights(l_S,1), rawWeights(l_S,2), 'ko', 'markerfacecolor', 'b');
plot(rawWeights(l_LvM,1), rawWeights(l_LvM,2), 'ko', 'markerfacecolor', 'r');
plot(rawWeights(l_SwL,1), rawWeights(l_SwL,2), 'ko', 'markerfacecolor', 'm');
plot(rawWeights(l_SwM,1), rawWeights(l_SwM,2), 'ko', 'markerfacecolor', 'g');
plot([0, 1], [1, 0], 'k')
plot([0, -1], [1, 0], 'k')
plot([0 0], [1 0], 'k')
plot([1 -1], [0 0], 'k')
plot([-.5 .5], [.5 .5], 'k+', 'markersize', 15)
axis equal
hold off
title('Pref Isolum Color')

%now according to CSI
l_CO = rawCSIs(:) >= 2;%color only
l_CL = (rawCSIs(:) >=0.5) & (rawCSIs(:) < 2); %color luminance
l_LO = rawCSIs(:) < 0.5;
figure, hold on,
plot(rawWeights(l_CO,1), rawWeights(l_CO,2), 'ko', 'markerfacecolor', 'r');
plot(rawWeights(l_CL,1), rawWeights(l_CL,2), 'ko', 'markerfacecolor', 'y');
plot(rawWeights(l_LO,1), rawWeights(l_LO,2), 'ko', 'markerfacecolor', 'k');
plot([0, 1], [1, 0], 'k')
plot([0, -1], [1, 0], 'k')
plot([0 0], [1 0], 'k')
plot([1 -1], [0 0], 'k')
plot([-.5 .5], [.5 .5], 'k+', 'markersize', 15)
axis equal
hold off
title('Pref Isolum Color')

%% Pref Color vs. CSI & SF

%bar chart depicting the percentage of CSI types for each pref isolum color
figure, hold on,
bar(1,sum([l_CO & l_S])./sum(l_S),0.9, 'r')
bar(2,sum([l_CL & l_S])./sum(l_S),0.9, 'y')
bar(3,sum([l_LO & l_S])./sum(l_S),0.9, 'k')

bar(5,sum([l_CO & l_SwL])./sum(l_SwL),0.9, 'r')
bar(6,sum([l_CL & l_SwL])./sum(l_SwL),0.9, 'y')
bar(7,sum([l_LO & l_SwL])./sum(l_SwL),0.9, 'k')

bar(9,sum([l_CO & l_SwM])./sum(l_SwM),0.9, 'r')
bar(10,sum([l_CL & l_SwM])./sum(l_SwM),0.9, 'y')
bar(11,sum([l_LO & l_SwM])./sum(l_SwM),0.9, 'k')

bar(13,sum([l_CO & l_LvM])./sum(l_LvM),0.9, 'r')
bar(14,sum([l_CL & l_LvM])./sum(l_LvM),0.9, 'y')
bar(15,sum([l_LO & l_LvM])./sum(l_LvM),0.9, 'k')

bar(17,sum(l_CO)./length(rawCSIs),0.9, 'r')
bar(18,sum(l_CL)./length(rawCSIs),0.9, 'y')
bar(19,sum(l_LO)./length(rawCSIs),0.9, 'k')
legend('Color Only', 'Color-Lum', 'Lum-Only')
set(gca, 'xtick', [2 6 10 14 18], 'XTickLabel', {'S-iso', 'SwL', 'SwM', 'L-M', 'All Cells'})
ylim([0,1])
title('Cell Type by Prefered Color')



% bar chart depicting the percentage of color types for each CSI type
figure, hold on,
bar(1,sum([l_CO & l_S])./sum(l_CO),0.9, 'b')
bar(2,sum([l_CO & l_SwL])./sum(l_CO),0.9, 'm')
bar(3,sum([l_CO & l_SwM])./sum(l_CO),0.9, 'g')
bar(4,sum([l_CO & l_LvM])./sum(l_CO),0.9, 'r')

bar(6,sum([l_CL & l_S])./sum(l_CL),0.9, 'b')
bar(7,sum([l_CL & l_SwL])./sum(l_CL),0.9, 'm')
bar(8,sum([l_CL & l_SwM])./sum(l_CL),0.9, 'g')
bar(9,sum([l_CL & l_LvM])./sum(l_CL),0.9, 'r')

bar(11,sum([l_LO & l_S])./sum(l_LO),0.9, 'b')
bar(12,sum([l_LO & l_SwL])./sum(l_LO),0.9, 'm')
bar(13,sum([l_LO & l_SwM])./sum(l_LO),0.9, 'g')
bar(14,sum([l_LO & l_LvM])./sum(l_LO),0.9, 'r')
legend('S-iso', 'SwL', 'SwM', 'L-M')
set(gca, 'xtick', [2.5 7.5 12.5], 'XTickLabel', {'Color Only', 'Color Lum', 'Lum Only'})
ylim([0,1])
title('Pref Color by CSI Type')


%% MODULATION RATIOS

% (1) first a simple histograms
figure
hist(rawModRatio, 35)
xlabel('Modulation Ratio')


% (2) modratio vs. CSI
l_nanModRatio = isnan(rawModRatio);
filtCSIs = rawCSIs(~l_nanModRatio);
l_CO = filtCSIs >=2;
l_CL = (filtCSIs >= 0.5) & (filtCSIs < 2);
l_LO = filtCSIs < 0.5;
valForInf = max(filtCSIs(~isinf(filtCSIs))) .* 1.05;
filtCSIs(isinf(filtCSIs)) = valForInf;
filtModRatio = rawModRatio(~l_nanModRatio);
[rho, p] = corr(filtCSIs(:), filtModRatio(:));
figure
hold on,
plot(filtCSIs(l_CO), filtModRatio(l_CO), 'ko', 'markerfacecolor', 'r')
plot(filtCSIs(l_CL), filtModRatio(l_CL), 'ko', 'markerfacecolor', 'Y')
plot(filtCSIs(l_LO), filtModRatio(l_LO), 'ko', 'markerfacecolor', 'k')
set(gca, 'xscale', 'log', 'yscale', 'log')
title(sprintf('spearmans rho=%.3f, p=%.3f', rho, p));
xlabel('Color Sensitivity Index')
ylabel('Modulation Ratio')

% (3) Mod Ratio for each color... but only consider CO anc CL cells
l_CO = rawCSIs >=2;
l_CL = (rawCSIs >= 0.5) & (rawCSIs < 2);
l_LO = rawCSIs < 0.5;
data = [rawModRatio(l_LvM & (l_CO|l_CL));rawModRatio(l_SwM & (l_CO|l_CL));rawModRatio(l_SwL & (l_CO|l_CL));rawModRatio(l_S & (l_CO|l_CL))];
figure
group = char(repmat('LvM', sum(l_LvM & (l_CO|l_CL)), 1), repmat('SwM', sum(l_SwM & (l_CO|l_CL)), 1), repmat('SwL', sum(l_SwL & (l_CO|l_CL)), 1), repmat('S-iso', sum(l_S & (l_CO|l_CL)), 1));
boxplot(data, group, 'notch', 'on')
ylabel('Modulation Ratio')


figure, hold on,
plot(ones(sum(l_LvM)), rawModRatio(l_LvM), 'ro')
plot(ones(sum(l_SwM))*2, rawModRatio(l_SwM), 'go')
plot(ones(sum(l_SwL))*3, rawModRatio(l_SwL), 'mo')
plot(ones(sum(l_S))*4, rawModRatio(l_S), 'bo')









%% Spatial Frequency

% (1) histo of SFs
figure
hist(rawSF, 25)
xlabel('Preferred Spatial Frequency')
ylabel('Counts')


% (2) Anova on SF by Pref Color and Box Plot
anovaData = [rawSF(l_LvM); rawSF(l_SwM); rawSF(l_SwL); rawSF(l_S)];
anovaGroup = [repmat('LvM ', sum(l_LvM),1); repmat('SwM ', sum(l_SwM),1); repmat('SwL ', sum(l_SwL),1); repmat('Siso', sum(l_S),1)];
p = anova1(anovaData, anovaGroup);
ylabel('Spatial Frequency')


% (3) SF broken down by color (plotting all the data)
figure
hold on,
plot(ones(sum(l_LvM)), rawSF(l_LvM), 'ro')
plot(ones(sum(l_SwM)).*2, rawSF(l_SwM), 'go')
plot(ones(sum(l_SwL)).*3, rawSF(l_SwL), 'mo')
plot(ones(sum(l_S)).*4, rawSF(l_S), 'bo')
plot(1, median(rawSF(l_LvM)), 'k+', 'markersize', 10, 'linewidth', 4)
plot(2, median(rawSF(l_SwM)), 'k+', 'markersize', 10, 'linewidth', 4)
plot(3, median(rawSF(l_SwL)), 'k+', 'markersize', 10, 'linewidth', 4)
plot(4, median(rawSF(l_S)), 'k+', 'markersize', 10, 'linewidth', 4)
xlim([0 5])
ylabel('Spatial Frequency')
xlabel('Pref Color')


% (4) SF as a function of CSI
figure
semilogy(rawSF, rawCSIs, 'ko')
ylabel('Color Sensitivity Index')
xlabel('Pref Spatial Frequency')







