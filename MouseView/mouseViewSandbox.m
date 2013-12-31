%% experimenting with unpacking an abf file

fin

fileName = '13d20004';
mouse = 'EMX_3';
datPath = [GL_DATPATH, mouse, filesep, 'Physiology'];
filePath = findfile(fileName, datPath, '.abf');


% open the file

[dat, ~, h] = abfload(filePath);
h.nExperimentType

%
% plot all the sweeps, or the one sweep if it's gap free.
%
%%%%%%%%%%

% define the time base
sampRate = 1/(h.si.*10^(-6));
N = size(dat,1);
tt = [0:N-1]./sampRate;

sweep = 1;
figure
nplts = numel(h.recChNames);
for a = 1:nplts;
    subplot(nplts, 1, a)
    plot(tt, dat(:,a,sweep))
    ylabel(h.recChNames{a})
    if a == 1
        xlabel('Time (ms)')
    end
end

% seems like ABF load is applying "parameter" related header info based off
% of only one recorded channel. Is it possible to output the V/I clamp
% value for each analog input channel?

