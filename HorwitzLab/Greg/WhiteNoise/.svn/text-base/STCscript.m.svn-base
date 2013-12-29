% A front end for STCGUi.m
% For loading one or two files.
%filename = {'K052708001'};  % Blue/Yellow
%filename = {'K102809004'};  % Double-opponent
%filename = {'K071709009'};    % Red-green with luminance PC1?
%filename = {'K021709003','K021709005'};    % Complex
%filename = {'K021609004','K021609005','K021609008'};
%filename = {'K010809002','K010809004','K010809008'};  % Blue-yellow
%filename = {'K031209001.2'};  % Blue-yellow (DOuble opponent)
%filename = {'K042309006','K042309008'};    % Color opponent
filename = {'S051410001.1'};    % Color opponent

 stro = {};
 for i = 1:length(filename)
     tmpstro = nex2stro(findfile(char(filename{i})));
     if (isempty(stro))
         stro = tmpstro;
     else
         stro = strocat(stro, tmpstro);
     end
 end
i = 1;
spikenum = 1;
noisetype = 'GUN';

% For loading a whole mess of files
%[fnames, spikenum] = fnamesFromTxt2('C:\NO BACKUP\NexFiles\nexfilelists\Greg\BYcandidates.txt');

data = {};
for i = 1:length(fnames)
    i
    stro=nex2stro(findfile(fnames{i}));
    spikeidx = find(strcmp(stro.sum.rasterCells(1,:),spikenum(i)));
    
    framerate = stro.sum.exptParams.framerate;
    nstixperside = stro.sum.exptParams.nstixperside;
    ntrials = length(stro.sum.absTrialNum);
    stimonidx = find(strcmp(stro.sum.trialFields(1,:),'stim_on'));
    stimoffidx = find(strcmp(stro.sum.trialFields(1,:),'all_off'));
    nframesidx = find(strcmp(stro.sum.trialFields(1,:),'num_frames'));
    noisetypeidx = find(strcmp(stro.sum.trialFields(1,:),'noise_type'));
    sigmaidxs = strmatch('sigma',stro.sum.trialFields(1,:));
    
%     hepidx = find(strcmp(stro.sum.rasterCells(1,:),'AD11'));
%     vepidx = find(strcmp(stro.sum.rasterCells(1,:),'AD12'));
%     anlgStartTimeidx = find(strcmp(stro.sum.rasterCells(1,:),'anlgStartTime'));
%     eyestart_t = [stro.ras{:,anlgStartTimeidx}]';
%     eyesampperiod = 1/stro.sum.analog.storeRates{1};
    
    gammaTable = stro.sum.exptParams.gamma_table;
    gammaTable = reshape(gammaTable, length(gammaTable)/3, 3);
    gammaTable1 = interp1(linspace(0,255,256),gammaTable,linspace(0,255,65536), 'spline');
    invgamma = InvertGamma(gammaTable, 0);
    
    % Reconstructing the M matrix
    fundamentals = stro.sum.exptParams.fundamentals;
    fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
    mon_spd = stro.sum.exptParams.mon_spd;
    mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
    mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
    M = fundamentals'*mon_spd;
    
    % Getting the background rgb/lms
    ridx = find(strcmp(stro.sum.trialFields(1,:),'bkgnd_r'));
    gidx = find(strcmp(stro.sum.trialFields(1,:),'bkgnd_g'));
    bidx = find(strcmp(stro.sum.trialFields(1,:),'bkgnd_b'));
    bkgndRGB = [mode(stro.trial(:,ridx)), mode(stro.trial(:,gidx)), mode(stro.trial(:,bidx))];
    bkgndrgb = [gammaTable(bkgndRGB(1)+1,1); gammaTable(bkgndRGB(2)+1,2); gammaTable(bkgndRGB(3)+1,3)];
    bkgndlms = M*bkgndrgb;
    
    % Stripping out everything but gun noise trials.
    maxT = 9;
    Lgunnoise = stro.trial(:,noisetypeidx) == 1;
    Lconenoise = stro.trial(:,noisetypeidx) == 2;
    if (strcmp(noisetype, 'CONE'))
        stro.ras(Lgunnoise,:) = [];
        stro.trial(Lgunnoise,:) = [];
    elseif (strcmp(noisetype,'GUN'))
        stro.ras(Lconenoise,:) = [];
        stro.trial(Lconenoise,:) = [];
    end
    
    out = getWhtnsStats(stro,maxT,'STCOVfull', {nstixperside^2, 3, 1, maxT}, stro.sum.rasterCells{spikenum(i)});
    STA = out{1};
    STC = out{2};
    STV = diag(STC);
    nspikes = out{3};

    data{i}.STA = STA;
    data{i}.STC = STC;
    data{i}.maxT = maxT;
    data{i}.nstixperside = nstixperside;
    if (exist('fnames'))
        data{i}.filename = char(fnames{i});
    else
        data{i}.filename = char(filename{1});
    end
end

save STCGUIdata data

%%
i = 1; STCGUI(data{i}.STA, data{i}.STC, data{i}.nstixperside, data{i}.maxT, data{i}.filename)
