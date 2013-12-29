
% Test of spatial homogeniety
% For now just painting the screen all one color to look for gross
% violations of homogeneiety.  Someday incorporate PR650 measurements and a
% target that moves among 9 locations.

NPIX = 100;
NREPS = 2;
NLOCS = 3;
NWAVELENGTHS = 81;
RGBs = [255 0 0; 0 255 0; 0 0 255; 255 255 255];

% Turning 
input('Turn on PR650 and then hit return now');
PR650init(1,0);

% Open the window and set up drawing rects
[window,screenRect] = Screen(0,'OpenWindow',[0 0 0]);
oldClut = Screen('ReadNormalizedGammaTable', window);
Screen('LoadNormalizedGammaTable',window,linspace(0,1,256)'*ones(1,3));

width = screenRect(3);
height = screenRect(4);
bkgndtex = Screen('MakeTexture', window, zeros(height,width,3)); 
Screen('DrawTexture', window, bkgndtex, [],[], [], 0);
xs = round(linspace(0,width-NPIX,NLOCS));
ys = round(linspace(0,height-NPIX,NLOCS));
tmp = fullfact([length(xs) length(ys)]);
rects = [xs(tmp(:,1))' ys(tmp(:,2))'];
rects = [rects, rects+repmat(NPIX,size(tmp,1),2)];

for i = 1:size(rects,1)
    data{i}.rect = rects(i,:);
    data{i}.spds = [];
end

% Main loop
for i = 1:NREPS
   for j = randperm(size(rects,1))
      spds = zeros(NWAVELENGTHS, size(RGBs,1));
      for k = randperm(size(RGBs,1))
        Screen('FillRect', window, RGBs(k,:), rects(j,:));
        Screen('Flip', window);
        if (all(spds(:) == 0))
            input('');
        end
        spds(:,k) = PR650measspd;
      end
      data{j}.spds = cat(3,data{j}.spds,spds);
   end
end

Screen('LoadNormalizedGammaTable', window, oldClut);
Screen('Close', window);

% Data analysis
% Change this so that we're looking at the position of each light in the
% CIE chromaticity diagram?  Look at the ratios?

xys= [];
whiteindex = find(all(RGBs == 255,2));
load '/Applications/Psychtoolbox/PsychColorimetricData/PsychColorimetricMatFiles/T_xyz1964.mat'
figure; axes; hold on;
xlabel('x'); ylabel('y'); title('White chromaticity');
for i = 1:size(data,2)
    tmp = data{i}.spds;
    white = squeeze(tmp(:,whiteindex,:));
    XYZ = T_xyz1964*white;
    x = XYZ(1,:)./sum(XYZ);
    y = XYZ(2,:)./sum(XYZ);
    plot(x,y);
    text(x(1),y(1),num2str(i));
end


% Preprocessing the data
integrals = zeros(NLOCS, NLOCS, size(RGBs,1));
for i = 1:size(data,2)
    x = data{i}.rect(1);
    y = data{i}.rect(2);
    spds = data{i}.spds;
    integral = [];
    for j = 1:size(spds,2)
        integrals(find(xs == x), find(ys == y), j) = mean(sum(squeeze(spds(:,j,:))));  
    end
end
% Doing the plotting
figure; colormap gray;
whiteindex = find(all(RGBs == 255,2));
for i = 1:size(integrals,3)
    normalizedgunvals = squeeze(integrals(:,:,i))./squeeze(integrals(:,:,whiteindex))
    subplot(2,size(integrals,3),i);
    imagesc(normalizedgunvals);
    axis image;
    set(gca,'XTick',[],'YTick',[]);
    contrast = (max(normalizedgunvals(:)) - min(normalizedgunvals(:)))./min(normalizedgunvals(:))
    cv = std(normalizedgunvals(:))./mean(normalizedgunvals(:));
    title(num2str(contrast));
end

