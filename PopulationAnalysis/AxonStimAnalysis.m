%% EB_031014_B Cell 1

fin


file_DCsteps = '2014_03_28_0000';
photo = 'cell_1_tdTomato_5x';
photoPath = findfile(photo, [GL_DATPATH, filesep, 'EB_031014_B'], '.jpg')

expts = {'2014_03_28_0009', 'distal pos 1 FS closed', [-179 -48];...
         '2014_03_28_0011', 'distal pos 1 FS open', [-179 -48];...
         '2014_03_28_0013', 'distal pos 2', [-288 -118];...
         '2014_03_28_0016', 'Soma FS closed', [0 0];...
         '2014_03_28_0019', 'distal pos 3', [-178 -169];...
         '2014_03_28_0020', 'distal pos 3 half moon', [-178 -169];...
         '2014_03_28_0021', 'distal pos 4 half moon', [-153 -147];...
         '2014_03_28_0022', 'distal pos 5 half moon', [-113 -111]};

Vhold = -85;
validCh = 'HS2_';
pixperum = 131.746 ./ 329.6195;



%% EB_031014_B Cell 2

fin


file_DCsteps = '2014_03_28_0023';
photo = 'cell_2_tdTomato_5x_470nm';
photoPath = findfile(photo, [GL_DATPATH, filesep, 'EB_031014_B'], '.jpg')

expts = {'2014_03_28_0026', 'soma FS open', [0 0];...
         '2014_03_28_0027', 'soma FS closed', [0 0];...
         '2014_03_28_0029', 'Lateral site 1', [215 -51];...
         '2014_03_28_0030', 'Lateral site 2', [176 117];...
         '2014_03_28_0032', 'Lateral site 3', [299 409];...
         '2014_03_28_0036', 'Lateral site 4', [-70 291];...
         '2014_03_28_0037', 'HalfMoon Incorrect', [-70 291];...
         '2014_03_28_0039', 'HalfMoon Correct', [-70 291]};

Vhold = -85;
validCh = 'HS1_';
pixperum = 131.746 ./ 329.6195;



%% EB_031014_E Cell 3: Single pulses

% check to make sure none of the files analyzed use a broad range of LED
% volatages

fin


file_DCsteps = '2014_03_25_0002';
photo = 'cell_3_tdTomato';
photoPath = findfile(photo, [GL_DATPATH, filesep, 'EB_031014_E'], '.jpg')

expts = {'2014_03_25_0007', 'Cortex 1', [-257 201];...
         '2014_03_25_0008', 'Cortex 2', [-357 291];...
         '2014_03_25_0009', 'Axons 1', [-678 384];...
         '2014_03_25_0010', 'Axons 2', [-833 48];...
         '2014_03_25_0011', 'Cortex 3 Slow', [-302 29];...
         '2014_03_25_0012', 'Cortex 3 Fast', [-302 29];...
         '2014_03_25_0014', 'Cortex 3 Half Moon', [-302 29];...
         '2014_03_25_0015', 'Cortex 3 No Half Moon', [-302 29]};

Vhold = -85;
validCh = 'HS2_';
pixperum = 131.746 ./ 329.6195;


%% EB_031014_E Cell 3: Pulse trains

% check to make sure none of the files analyzed use a broad range of LED
% volatages

fin


file_DCsteps = '2014_03_25_0002';
photo = 'cell_3_tdTomato';
photoPath = findfile(photo, [GL_DATPATH, filesep, 'EB_031014_E'], '.jpg');

expts = {'2014_03_25_0016', 'Cortex 3 Half Moon', [-302 29];...
         '2014_03_25_0018', 'Soma FS Open', [0 0]};

Vhold = -85;
validCh = 'HS2_';
pixperum = 131.746 ./ 329.6195;
%% Run the analysis


% colors for various plots
map = colormap('jet'); close;
clrIdx = round(linspace(1,size(map,1), size(expts,1)));


% plot the DC steps routine
ax = abfobj(file_DCsteps);
idx_Vm = eval(['ax.idx.',validCh,'Vm']);
idx_Iclamp = eval(['ax.idx.',validCh,'Iclamp']);
figure
plot(ax.tt, permute(ax.dat(:, idx_Vm,:), [1,3,2]))

% show an image of the slice
img = imread(photoPath);
figure
imshow(img);
centPos = round(ginput(1));
stimPoints = cell2mat([expts(:,3)]);
stimPoints = round(stimPoints .* pixperum); %now in pix
stimPoints = bsxfun(@plus, stimPoints, centPos); % pix relative to neuron
%stimPoints = fliplr(stimPoints); % swap x,y and r,c notation
%stimPoints(:,1) = size(img, 1) - stimPoints(:,1);
hold on,
for a = 1:size(stimPoints,1)
    plot(stimPoints(a,1), stimPoints(a,2), 'o', 'markeredgecolor', map(clrIdx(a),:), 'markerfacecolor', map(clrIdx(a),:))
end
drawnow



% iterate over the whole cell recordings. check to make sure it's at the
% correct holding potential, then save the average
avgCurrent = [];
avgCmd = [];
access = [];
for a = 1:size(expts,1)
    
    disp(a)
    
    %load the data
    ax = abfobj(expts{a,1});
    idx_Im = eval(['ax.idx.',validCh,'Im']);
    idx_Vclamp = eval(['ax.idx.',validCh,'Vclamp']);
    idx_SecVm = eval(['ax.idx.',validCh,'secVm']);
    
    % check the Vhold
    actVhold = ax.dat(1:100, idx_SecVm, :);
    actVhold = mean(mean(actVhold, 3));
    if (actVhold-Vhold) > 1; error('Vhold is wrong'); end
    
    % pull out the data
    raw = mean(ax.dat(:,idx_Im,:),3);
    raw = raw - mean(raw(1:100));
    avgCurrent(:,a) = raw;

    % get the series resistanc
    Ra = ax.getRaRin('linear');
    access = cat(1, access, Ra(:));
end


figure, hold on
for a = 1:size(expts,1)
    plot(ax.tt,  avgCurrent(:,a), 'color', map(clrIdx(a),:), 'linewidth', 2)
end
leg = expts(:,2);
legend(leg')
xlabel('Time (sec)')
ylabel('Baseline subtracted current (pA)')


figure
plot(access)
ylim([0 max(access).*1.05])
xlabel('Sweep Number')
ylabel('Series Resistance (mOhm)')






