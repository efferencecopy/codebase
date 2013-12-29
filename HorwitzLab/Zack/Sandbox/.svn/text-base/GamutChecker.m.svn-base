%% generate color dirs
% define the color directions.... just putting points on a sphere for now
nColors = 2000;
tmp = ceil(sqrt(nColors));
az = linspace(0, 2*pi, tmp);
el = linspace(0, pi/2, tmp);
inds = fullfact([tmp, tmp]);
dirs_pol = [az(inds(:,1))', el(inds(:,2))'];
[x,y,z] = sph2cart(dirs_pol(:,1), dirs_pol(:,2), ones(size(dirs_pol,1),1));
%z = z.*3; % scale the s-cone axis
colorDirs = [x,y,z];
norms = sqrt(sum(colorDirs.^2,2));
colorDirs = bsxfun(@rdivide, colorDirs, norms);
colorDirs(abs(colorDirs)<1e-14) = 0;
colorDirs = [colorDirs; 1 -1 0; -1 1 0; 0 0 -1; 0 0 1];
colorDirs = unique(colorDirs, 'rows'); %remove the duplicates
nColors = size(colorDirs,1);

lvsm_idxs = find(abs(colorDirs(:,1)) == 1 & abs(colorDirs(:,2)) == 1 & colorDirs(:,3) == 0);
siso_idxs = find(colorDirs(:,1) == 0 & colorDirs(:,2) == 0 & abs(colorDirs(:,3)) == 1);

%%
load Dell4BitsCal.mat
s = load('T_cones_smj');
fns = fieldnames(s);
fundamentals = s.(fns{1});
fundWavelengthSpacing = s.(fns{2});
calData = cals{end};
ptb.bkgndRGB = round(255 .* calData.bgColor); %bkgnd voltages discretized b/w 0&255
ptb.bkgndrgb = [calData.gammaTable(ptb.bkgndRGB(1)+1, 1), calData.gammaTable(ptb.bkgndRGB(2)+1, 2), ...
    calData.gammaTable(ptb.bkgndRGB(3)+1, 3)]; %add one to create and index b/w 1 & 256

ptb.maxdac = 2^16-1;
ptb.gammaTable = calData.gammaTable;
ptb.monSpd = SplineSpd(calData.S_device, calData.P_device, fundWavelengthSpacing);
ptb.M = fundamentals * ptb.monSpd;
ptb.invM = inv(ptb.M);
ptb.bkgndlms = ptb.M * ptb.bkgndrgb';
ptb.invGamma = InvertGammaTable(calData.gammaInput, ptb.gammaTable, ptb.maxdac+1);

for c = 1:size(colorDirs, 1);
    guess = 0;
    incSize = 0.1; %1 tenth of a %CC
    while(1)
        tmpLMS = colorDirs(c,:) .* guess;
        
        %peak of the gabor:
        peakrgb = ptb.invM*((1+tmpLMS./100)' .* ptb.bkgndlms);
        peakRGB = round(ptb.maxdac .* peakrgb) + 1;
        peakInGamut = all(peakRGB < (ptb.maxdac+1)) && all(peakRGB > 0);
        
        %trough of the gabor:
        troughrgb = ptb.invM*((1-tmpLMS./100)' .* ptb.bkgndlms);
        troughRGB = round(ptb.maxdac .* troughrgb) + 1;
        troughInGamut = all(troughRGB < (ptb.maxdac+1)) && all(troughRGB > 0);
        
        if (peakInGamut && troughInGamut)
            guess = guess + incSize;
        else
            maxContrast(c) = guess - incSize;
            break
        end
    end
end

% the _diffs stuff is because WhiteNoise uses units of difference from background.
% gamut_in_diffs = colorDirs .* bsxfun(@times, ptb.bkgndlms, maxContrast/100)';
gamut_in_cc = bsxfun(@times,colorDirs,maxContrast'/100);

fprintf('%s:\n', cals{end}.describe.monitor);
fprintf('\tL-M gamut length: %g\n', sqrt(sum(diff(gamut_in_cc(lvsm_idxs,:)).^2)));
fprintf('\tS-iso gamut length: %g\n', sqrt(sum(diff(gamut_in_cc(siso_idxs,:)).^2)));

% figure; hold on;
% ourcolors = 'br';
% for i = 0:1
%     plot3((2*i-1)*gamut_in_diffs(: ,1),(2*i-1)*gamut_in_diffs(: ,2),(2*i-1)*gamut_in_diffs(: ,3),[ourcolors(i+1) '.']);
%     plot3((2*i-1)*gamut_in_diffs(: ,1),(2*i-1)*gamut_in_diffs(: ,2),(2*i-1)*gamut_in_diffs(: ,3),[ourcolors(i+1) '.']);
% end
% title('diffs')
% xlabel('L'); ylabel('M'); zlabel('S');

% figure; hold on;
% ourcolors = 'br';
% for i = 0:1
%     plot3((2*i-1)*gamut_in_cc(: ,1),(2*i-1)*gamut_in_cc(: ,2),(2*i-1)*gamut_in_cc(: ,3),[ourcolors(i+1) '.']);
%     plot3((2*i-1)*gamut_in_cc(: ,1),(2*i-1)*gamut_in_cc(: ,2),(2*i-1)*gamut_in_cc(: ,3),[ourcolors(i+1) '.']);
% end
% 
% title('cc')
% xlabel('L'); ylabel('M'); zlabel('S');
