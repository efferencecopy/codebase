% Figures for a microsaccade paper
% 
% Section 1: Some summary statistics of the fixational saccades (amplitude,
% direction, main sequence)
%
% Section 2: Probability of correct given that a saccade was made at time
% 't' as a function of 't'.
%
% Section 3: Distribution of saccade directions and amplitudes
% as a function of correct and incorrect at different times in the trial.
%
% Section 4: Barplots showing effects of color direction on saccadic
% suppression. OUTDATED
%
% Section 5: Detection threshold performance summary.
% Section 5.1: Bitmap of a gabor.
%
% Section 6: Barplots showing effects of color direction on saccadic
% suppression summarized as changes in percent correct and threshold.
%
% Section 6.1 Line plots showing effects of color direction on saccadic
% suppression summarized as changes in percent correct and threshold
%
% Section 6.2 Bar plots a la section 6, but now with standard error bars.
%
% Section 6.3 Bar plots a la section 6, but now 2-D
%
% Section 7: STAs and saccade-triggered PSTHs for a few example neurons
%
% Section 8: A comparison of mean saccade-triggered PSTHs across a
% population of color-opponent and nonopponent cells
%
% Section 8.1: Additional panel showing microsaccade-triggered PSTHs
% for L-M DTspot data.
%
% Section 9: Percent correct as a function of saccade amplitude (following
% the suggestion of reviewer 1).
%
% Section 10: Percent correct as a function of time separately for large
% and small amplitude saccades.
%
% Section 11: Charlie's code experimenting with different ways of
% normalizing post-saccade spike density functions.
%
% Section 12: An animated gif showing fixational eye movements
%
%%
% Section 1
% A few summary plots
% Largely take from Section 1 of DTEManalysis.m

% Setting up figure
figure('Units','inches','PaperPosition',[0 0 8.5 11],'Position',[0 0 8.5 11],'PaperOrientation','Portrait','DefaultAxesUnits','inches');
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultTextFontName','Helvetica','DefaultTextFontSize',12);
colordef white;

AMPLITUDEBINS = linspace(0,1,20);
AMPLITUDEBINWIDTH = AMPLITUDEBINS(2)-AMPLITUDEBINS(1);
AXESSIZE = 1.7;
filenamelist = {'N:\NexFiles\nexfilelists\Greg\DTEM\KaliQuestLrg.txt',...
    'N:\NexFiles\nexfilelists\Greg\DTEM\SednaQuestLrg.txt'};

for a = 1:length(filenamelist)
    amplitudes = [];
    directions = [];
    peakv = [];
    filenames = fnamesFromTxt(filenamelist{a});
    for b = 1:size(filenames,1)
        stro = nex2stro(findfile(filenames(b,:)));
        stro = DTfilterquesttrials(stro,'PaperDefaults');
        if isempty(stro.sum.analog.sigid) || isempty(stro.trial)
            continue;
        end
        ntrials = size(stro.trial,1);
        stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_on'));
        stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_off'));
        correct = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'correct'));
        sacstats = getSacData(stro);
        close;
        
        for i = 1:ntrials
            st = sacstats.starttimes{i};
            Lsac = (sacstats.starttimes{i} > stimon_t(i)) & (sacstats.endtimes{i} < stimoff_t(i));
            if any(Lsac)
                if any(sacstats.amplitudes{i}(Lsac) > 2)
                    keyboard
                end
                amplitudes = [amplitudes; sacstats.amplitudes{i}(Lsac)];
                directions = [directions; sacstats.directions{i}(Lsac)];
                peakv = [peakv; sacstats.peakv{i}(Lsac)];
            end
        end
    end
    
    % Amplitude histogram
    axes('Position',[1+3*(a-1) 5.5 AXESSIZE AXESSIZE]);
    [n,x] = histc(amplitudes,AMPLITUDEBINS);
    h = bar(AMPLITUDEBINS,n,'histc');
    set(h,'FaceColor','black','EdgeColor',[1 1 1])
    set(gca,'XLim',[AMPLITUDEBINS(1) AMPLITUDEBINS(end)])
    set(gca,'TickDir','out')
    set(gca,'XTick',[0 .25 .5 .75 1],'XTickLabel',[0 .25 .5 .75 1]);
 
    xlabel('Amplitude (deg)','FontSize',12); ylabel('Count','FontSize',12);
    
    % Main sequence
    axes('Position',[1+3*(a-1) 3 AXESSIZE AXESSIZE]);
    [n,x] = hist2([amplitudes,peakv],[50 50; 0 10; 1 60]);
    mean(amplitudes)
    imagesc(flipud(n)); 
    colormap(1-gray.^.5); 
    colormap(1-gray); 
    axis xy;
    hold on;
 %   plot([1 42],[1 50],'k:');

    % Need these numbers for the equation in DTEManalysis section 9
    x{1}([1 42])
    x{2}([1 50])
    
%     % Sanity check: omitting 'curved' saccades
%     a = amplitudes;
%     v = peakv;
%     v_hat= (125-20)/.8367*a+20;
%     L = v < v_hat;
%     [n,x] = hist2([amplitudes(L),peakv(L)],[50 50; 0 20; 1 125]);
%     imagesc(flipud(n)); colormap(1-gray.^.5); axis xy;
%     hold on;
%     plot([1 42],[1 50],'k:');
%     % End of sanity check.
    
    countingwin = [x{1}(1) x{1}(end) x{2}(1) x{2}(end)];
    m = (x{1}(end)-x{1}(1))./diff(get(gca,'XLim')+[.5 -.5]);
    b = x{1}(1)-m;  % assuming first bin is at '1'
    set(gca,'XTick',([0 .25 .5 .75 1]-b)./m,'XTickLabel',[0 .25 .5 .75 1]);
    m = (x{2}(end)-x{2}(1))./diff(get(gca,'YLim')+[.5 -.5]);
    b = x{2}(1)-m;  % assuming first bin is at '1'
    bins = 25:25:150;
    set(gca,'YTick',(bins-b)./m,'YTickLabel',bins);
    xlabel('Amplitude (deg)','FontSize',12); ylabel('Peak speed (deg/sec)','FontSize',12);
    
    % Directions
    axes('Position',[1+3*(a-1) .5 AXESSIZE AXESSIZE]);
    [rho, theta]= hist(directions,20);
    if(a == 1)
        polar(0,2000,'w.');
    else
        polar(0,4000,'w.'); 
    end
    hold on;
    h = polar([theta theta(1)],[rho rho(1)],'k-');
    set(h,'Linewidth',2);
    
    % Monkey label
    axes('Position',[1+3*(a-1) 7.5 AXESSIZE 1],'Box','off','Visible','off');
    text(0.5,0,['Monkey ',filenames(1,1)],'HorizontalAlignment','center','FontSize',16)
end

%%
% Section 2
% Probability of correct | saccade at time t
% Peristimulus saccade time histogram

figure('Units','inches','PaperPosition',[0 0 8.5 11],'Position',[0 0 8.5 11],'PaperOrientation','Portrait','DefaultAxesUnits','inches');
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultTextFontName','Helvetica','DefaultTextFontSize',12);
colordef white;

AXESSIZE = 1.7;
filenamelist = {'N:\NexFiles\nexfilelists\Greg\DTEM\KaliQuestLrg.txt',...
    'N:\NexFiles\nexfilelists\Greg\DTEM\SednaQuestLrg.txt'};
binwidth = .05; %sec
bins = [-.2:binwidth:.767];
axishandles = [];
ntrialstot = [0 0 0 0]; % Kali, Sedna, Kali no frames, Sedna no frames
for a = 1:length(filenamelist)
    x = zeros(length(bins),1);
    n = zeros(length(bins),1);
    filenames = fnamesFromTxt(filenamelist{a});
    for b = 1:size(filenames,1)
        stro = nex2stro(findfile(filenames(b,:)));
        stro = DTfilterquesttrials(stro,'PaperDefaults');

        if isempty(stro.sum.analog.sigid) || isempty(stro.trial)
            continue;
        end
        
        ntrials = size(stro.trial,1);
        stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_on'));
        stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_off'));
        targon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targ_on'));
        correct = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'correct'));
        
        sacstats = getSacData(stro);
        close;
        
        % Looping over trials
        for k = 1:length(bins)
            for j = 1:ntrials
                dirs = sacstats.directions{j};
                amps = sacstats.amplitudes{j};
                st = sacstats.starttimes{j};
                L = logical(sacstats.starttimes{j} >= targon_t(j));
                st(L) = [];    % Omitting saccades that begin after targets appear
                amps(L) = [];  % Omitting saccades that begin after targets appear
                dirs(L) = [];  % Omitting saccades that begin after targets appear
                Lsac = (st > stimon_t(j)+bins(k)-binwidth/2) & (st < stimon_t(j)+bins(k)+binwidth/2);
                if (any(Lsac))
                    if (correct(j))
                        x(k) = x(k)+1;
                    end
                    n(k) = n(k)+1;
                end
            end
        end
        ntrialstot(a) = ntrialstot(a)+ntrials;
    end
    
    % Saccade frequency v time
    axishandles = [axishandles axes('Position',[1+3*(a-1) 5.5 AXESSIZE AXESSIZE])];
    plot(bins, n./(binwidth*ntrialstot(a)),'k-','Linewidth',2)
    ylabel('Rate (saccades/sec)','FontSize',12);
    set(gca,'Xlim',[bins(1) bins(end)+eps],'Box','off');
     
    % P(correct) v time
    axes('Position',[1+3*(a-1) 3 AXESSIZE AXESSIZE]);
    hold on;
    p = sum(x,2)./sum(n,2);
    se = sqrt((p.*(1-p))./sum(n,2));
    h = patch([bins fliplr(bins)],[p+se; flipud(p-se)]','b');
    set(h,'FaceColor',[.5 .5 .5]);
    plot(bins, sum(x,2)./sum(n,2),'k-','Linewidth',2);
    plot([0 .167 .5 .666],[0.65 .75 .75 0.65],'k:');
    ylabel('P(cor | saccade at t)','FontSize',12);
    xlabel('Time (s)','FontSize',12);
    set(gca,'Xlim',[bins(1) bins(end)+eps])
    if (a == 1)
        set(gca,'YLim',[.64 .851]);
    else
       set(gca,'YLim',[.715 .831]);
    end
    
    % Monkey label
    axes('Position',[1+3*(a-1) 7.5 AXESSIZE 1],'Box','off','Visible','off');
    text(0.5,0,['Monkey ',filenames(1,1)],'HorizontalAlignment','center','FontSize',16);
end

% Adding a dashed line for microsaccade frequency in the "noframes" condition

filenamelist = {'N:\NexFiles\nexfilelists\Greg\DTEM\KaliNoFrames.txt',...
    'N:\NexFiles\nexfilelists\Greg\DTEM\SednaNoFrames.txt'};
for a = 1:length(filenamelist)
    filenames = fnamesFromTxt(filenamelist{a});
    n = zeros(length(bins),1);
    for b = 1:size(filenames,1)
        stro = nex2stro(findfile(filenames(b,:)));
        if isempty(stro.sum.analog.sigid) || isempty(stro.trial)
            continue;
        end
        
        ntrials = size(stro.trial,1);
        stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_on'));
        stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_off'));
        targon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targ_on'));
        correct = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'correct'));
        
        sacstats = getSacData(stro);
        close;
        
        % Looping over trials
        for k = 1:length(bins)
            for j = 1:ntrials
                dirs = sacstats.directions{j};
                amps = sacstats.amplitudes{j};
                st = sacstats.starttimes{j};
                L = logical(sacstats.starttimes{j} >= targon_t(j));
                st(L) = [];    % Omitting saccades that begin after targets appear
                amps(L) = [];  % Omitting saccades that begin after targets appear
                dirs(L) = [];  % Omitting saccades that begin after targets appear
                Lsac = (st > stimon_t(j)+bins(k)-binwidth/2) & (st < stimon_t(j)+bins(k)+binwidth/2);
                if (any(Lsac))
                    n(k) = n(k)+1;
                end
            end
        end
        ntrialstot(a+2) = ntrialstot(a+2)+ntrials;
    end
    
    % Saccade frequency v time
    axes(axishandles(a)); hold on;
    rate =n/(binwidth*ntrialstot(a+2));
    h = plot(bins, rate,'k:','Linewidth',1)
end

%%
% Section 3
% Distribution of saccade endpoints at different times in the trial
% preceding corrects/incorrects and T1/T2 choices.
ACHROMONLY = 0;
figure('Units','inches','PaperPosition',[0 0 8.5 11],'Position',[0 0 8.5 11],'PaperOrientation','Portrait','DefaultAxesUnits','inches');
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultTextFontName','Helvetica','DefaultTextFontSize',12);
colordef white;

AXESSIZE = 0.5;
AXESMARGIN = 0.1;
filenamelist = {'N:\NexFiles\nexfilelists\Greg\DTEM\SednaQuestLrg.txt',...
    'N:\NexFiles\nexfilelists\Greg\DTEM\KaliQuestLrg.txt'};

t_binwidth = .15; %sec
s_binwidth = .1; %deg
t_bins = [-.2:t_binwidth:.767];
s_bins = [-.5:s_binwidth:.5];
niter = 10000;
minnumbertrials = 10;
maxlim = 0;
alpha = 0.05;
for a = 1:length(filenamelist)
    cordata = zeros(length(s_bins),length(s_bins),length(t_bins));
    incdata = zeros(length(s_bins),length(s_bins),length(t_bins));
    T1data = zeros(length(s_bins),length(s_bins),length(t_bins));
    T2data = zeros(length(s_bins),length(s_bins),length(t_bins));
    
    filenames = fnamesFromTxt(filenamelist{a});
    for b = 1:size(filenames,1)
        stro = nex2stro(findfile(filenames(b,:)));
        stro = DTfilterquesttrials(stro,'PaperDefaults');

        if isempty(stro.sum.analog.sigid) || isempty(stro.trial)
            continue;
        end
        
        ntrials = size(stro.trial,1);
        
        if (ACHROMONLY)
            lms = stro.sum.exptParams.RF_colors;
            lms = reshape(lms,[3,size(lms,1)/3]);
            lms(:,all(lms == 0)) = [];
            lms = mkbasis(lms);
            whichcol = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'color_dir'));
        end
        
        stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_on'));
        stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_off'));
        targon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targ_on'));
        correct = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'correct'));
        flashside = -1*sign(stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_x')));
        flashside = max(flashside, 0); % 1 is to the left
        choice = (correct & flashside) | (~correct & ~flashside);
        sacstats = getSacData(stro);
        close;
        
        % Looping over trials
        for k = 1:length(t_bins)
            for j = 1:ntrials
                
                if (ACHROMONLY)
                    coloridx = find(abs(lms(:,whichcol(j))'*[.5774 .5774 .5774]') > .99);  % Avoiding roundoff errors
                    if (coloridx ~= 1)
                        continue;
                    end
                end
                
                dirs = sacstats.directions{j};
                amps = sacstats.amplitudes{j};
                st = sacstats.starttimes{j};
                L = logical(sacstats.endtimes{j} >= targon_t(j));
                st(L) = [];    % Omitting saccades that terminate after targets appear
                amps(L) = [];  % Omitting saccades that terminate after targets appear
                dirs(L) = [];  % Omitting saccades that terminate after targets appear
                Lsac = (st > stimon_t(j)+t_bins(k)-t_binwidth/2) & (st < stimon_t(j)+t_bins(k)+t_binwidth/2);
                if (any(Lsac))
                    [x,y] = pol2cart(dirs(Lsac),amps(Lsac));
                    xidxs = find(histc(x,s_bins));
                    yidxs = find(fliplr(histc(y,s_bins)));
                    if (correct(j))
                        cordata(yidxs,xidxs,repmat(k,sum(Lsac),1)) =  cordata(yidxs,xidxs,repmat(k,sum(Lsac),1)) + 1;
                    else
                        incdata(yidxs,xidxs,repmat(k,sum(Lsac),1)) =  incdata(yidxs,xidxs,repmat(k,sum(Lsac),1)) + 1;
                    end
                    if (choice(j))
                        T1data(yidxs,xidxs,repmat(k,sum(Lsac),1)) =  T1data(yidxs,xidxs,repmat(k,sum(Lsac),1)) + 1;
                    else
                        T2data(yidxs,xidxs,repmat(k,sum(Lsac),1)) =  T2data(yidxs,xidxs,repmat(k,sum(Lsac),1)) + 1;
                    end
                end
            end
        end
    end
    ntimepoints = length(t_bins);
    
    % For each pixel P(corr)
    data = cordata./(cordata+incdata);
    data(cordata+incdata < minnumbertrials) = 2;
    
    for frame = 1:length(t_bins) % Probability correct
        % Resampling test (assuming an empirical multinomial distribution
        % for endpoints).
        tmp = zeros(niter,1);
        ncor = sum(sum(cordata(:,:,frame)));
        ninc = sum(sum(incdata(:,:,frame)));
        empdistn = (cordata(:,:,frame)+incdata(:,:,frame))./(ninc+ncor);
        % Bug workaround from Mathworks site
        % (http://www.mathworks.com/support/bugreports/644205)
        e = min([0 cumsum(empdistn(:))'],1);
        e(end) = 1;
        empdistn = reshape(diff(e),size(empdistn));
        
        for iter = 1:niter
            if (iter == 1)
                resampcor = cordata(:,:,frame);
                resampinc = incdata(:,:,frame);
            else
                resampcor = reshape(mnrnd(ncor,empdistn(:)),size(empdistn));
                resampinc = reshape(mnrnd(ninc,empdistn(:)),size(empdistn));
            end
            resampcor = resampcor./sum(sum(resampcor));
            resampinc = resampinc./sum(sum(resampinc));
            mncor = [sum([length(s_bins):-1:1]*resampcor) sum([1:length(s_bins)]*resampcor')];
            mninc = [sum([length(s_bins):-1:1]*resampinc) sum([1:length(s_bins)]*resampinc')];
            tmp(iter) = norm(mncor-mninc);
        end
        p = sum(tmp > tmp(1))/niter;
        
        axes('Position',[(AXESSIZE+AXESMARGIN)*frame+2 2.75+3*(a-1) AXESSIZE AXESSIZE],'XTick',[],'YTick',[],'Box','off');
        imagesc(data(:,:,frame),[0 1]);
        axis square; hold on;
        
       % if (frame == length(t_bins))
            % Plotting distance between endpoint centroids
            empcordist = cordata(:,:,frame)./sum(sum(cordata(:,:,frame)));
            empincdist = incdata(:,:,frame)./sum(sum(incdata(:,:,frame)));
            mncor = [sum([1:length(s_bins)]*empcordist') sum(fliplr([1:length(s_bins)])*empcordist)]
            mninc = [sum([1:length(s_bins)]*empincdist') sum(fliplr([1:length(s_bins)])*empincdist)]
            overallmn = mean([mncor; mninc]);
            MNSCALEFACTOR = 10;
            mncor = MNSCALEFACTOR*(mncor-overallmn) + overallmn;
            mninc = MNSCALEFACTOR*(mninc-overallmn) + overallmn;
            mncor(2) = length(s_bins)-mncor(2)+1;
            mninc(2) = length(s_bins)-mninc(2)+1;
            
            plot([mncor(1) mninc(1)],[mncor(2) mninc(2)],'k-');
            plot(mncor(1),mncor(2),'wo','MarkerFaceColor','white','MarkerSize',5);
            plot(mninc(1),mninc(2),'ko','MarkerFaceColor','black','MarkerSize',5);
       % end
        
        title(num2str(t_bins(frame)));
        set(gca,'XTick',[],'YTick',[]);
        
        if ( p < alpha)
            text(length(s_bins)/2,length(s_bins)+3,'*');
        end
     end
    axes('Position',[2 2.75+3*(a-1) AXESSIZE AXESSIZE],'XTick',[],'YTick',[],'Box','off','Visible','off');
    text(0.3,0.5,'Cor/Inc','HorizontalAlignment','center','FontSize',12);
   
    % marginal histogram requested by Reviewer 1
  %  axes('Position',[(AXESSIZE+AXESMARGIN)*(frame+1)+2-AXESMARGIN 2.75+3*(a-1) AXESSIZE AXESSIZE],'XTick',[],'YTick',[],'Box','off');
  %  data(cordata+incdata < minnumbertrials) = nan;
  %  marginalmean = nanmean(data(:,:,frame),2);
  %  h = barh(fliplr(s_bins),marginalmean,'FaceColor','black');
  %  set(gca,'Ylim',[min(s_bins)-s_binwidth/2 max(s_bins)+s_binwidth/2],'Ytick',[],'Color','none');
  %  set(gca,'Xlim',[min(marginalmean)*.9 max(marginalmean)*1.1],'Xtick',[.7 .9],'Box','off');
    
    % Now, T1 vs T2
    data = T1data./(T1data+T2data);
    data(T1data+T2data < minnumbertrials) = 2;

    for frame = 1:ntimepoints  % T1 vs T2 choices
        % Resampling test (assuming an empirical multinomial distribution
        % for endpoints).
        tmp = zeros(niter,1);
        nT1 = sum(sum(T1data(:,:,frame)));
        nT2 = sum(sum(T2data(:,:,frame)));
        empdistn = (T1data(:,:,frame)+T2data(:,:,frame))./(nT1+nT2);
        % Bug workaround from Mathworks site
        % (http://www.mathworks.com/support/bugreports/644205)
        e = min([0 cumsum(empdistn(:))'],1);
        e(end) = 1;
        empdistn = reshape(diff(e),size(empdistn));
        
        for iter = 1:niter
            if (iter == 1)
                resampT1 = T1data(:,:,frame);
                resampT2 = T2data(:,:,frame);
            else
                resampT1 = reshape(mnrnd(nT1,empdistn(:)),size(empdistn));
                resampT2 = reshape(mnrnd(nT2,empdistn(:)),size(empdistn));
            end
            resampT1 = resampT1./sum(sum(resampT1));
            resampT2 = resampT2./sum(sum(resampT2));
            mncor = [sum([length(s_bins):-1:1]*resampT1) sum([1:length(s_bins)]*resampT1')];
            mninc = [sum([length(s_bins):-1:1]*resampT2) sum([1:length(s_bins)]*resampT2')];
            tmp(iter) = norm(mncor-mninc);
        end
        p = sum(tmp > tmp(1))/niter;
        
        axes('Position',[(AXESSIZE+AXESMARGIN)*frame+2 2+3*(a-1) AXESSIZE AXESSIZE],'XTick',[],'YTick',[],'Box','off');
        imagesc(data(:,:,frame),[0 1]);
        axis square;
        set(gca,'XTick',[],'YTick',[]);
        if (p < alpha)
            text(length(s_bins)/2,length(s_bins)+3,'*');
        end
    end
    axes('Position',[2 2+3*(a-1) AXESSIZE AXESSIZE],'XTick',[],'YTick',[],'Box','off','Visible','off');
    text(0.3,0.5,'T1/T2','HorizontalAlignment','center','FontSize',12);

    % Monkey name
    axes('Position',[(AXESSIZE+AXESMARGIN)+2 3.75+3*(a-1) ntimepoints*(AXESSIZE+AXESMARGIN) 1],'Box','off','Visible','off');
    text(0.5,0,['Monkey ',filenames(1,1)],'HorizontalAlignment','center','FontSize',16)
end


% Color scale bar
axes('Position',[7 4 .25 1.25]);
imagesc(flipud(linspace(0, 1,500)'),[0 1]);
set(gca,'XTick',[],'YTick',[1 250 500],'YTickLabel',[1 .5 0]);
ylabel('Probability','Fontsize',12);
set(gca,'YAxisLocation','right');
cmap = colormap(jet);
hs = get(gcf,'Children');
for i = 1:length(hs)
   set(hs(i),'Clim',[0 1.1]);
end
colormap([cmap; .5 .5 .5]);


%%
% Section 4
% Barplots showing effects of saccades on probability correct marginalized
% on color and spatial frequency.

% filenamelist = {'N:\NexFiles\nexfilelists\Greg\DTEM\SednaQuestLrg.txt',...
%     'N:\NexFiles\nexfilelists\Greg\DTEM\KaliQuestLrg.txt'};
% 
% figure('Units','inches','PaperPosition',[0 0 8.5 11],'Position',[0 0 8.5 11],'PaperOrientation','Portrait','DefaultAxesUnits','inches');
% set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
% set(gcf,'DefaultTextFontName','Helvetica','DefaultTextFontSize',12);
% colordef white;
% AXESSIZE = 1.5;
% 
% offsets = [.167 .5];   % relative to stim on, when to look for saccades
% colordirs = mkbasis([1 1 1; 0 0 1; 1 -1 0]');
% periods = [206 58 16];
% alpha = 0.05;
% for a = 1:length(filenamelist)
%     x_sac = zeros(length(colordirs), length(periods));
%     n_sac = zeros(length(colordirs), length(periods));
%     x_nosac = zeros(length(colordirs), length(periods));
%     n_nosac = zeros(length(colordirs), length(periods));
%     
%     filenames = fnamesFromTxt(filenamelist{a});
%     for b = 1:size(filenames,1)
%         stro = nex2stro(findfile(filenames(b,:)));
%         stro = DTfilterquesttrials(stro,'PaperDefaults');
% 
%         if isempty(stro.sum.analog.sigid) || isempty(stro.trial)
%             continue;
%         end
%         
%         ntrials = size(stro.trial,1);
%         lms = stro.sum.exptParams.RF_colors;
%         lms = reshape(lms,[3,size(lms,1)/3]);
%         lms(:,all(lms == 0)) = [];
%         lms = mkbasis(lms);
%         whichcol = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'color_dir'));
%         stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_on'));
%         stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_off'));
%         targon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targ_on'));
%         spatialPeriods = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'gabor_lambda'));
%         correct = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'correct'));
%         sacstats = getSacData(stro);
%         close;
%         
%         % Looping over trials
%         for j = 1:ntrials
%             coloridx = find(abs(lms(:,whichcol(j))'*colordirs) > .99);  % Avoiding roundoff errors
%             if (isempty(coloridx))
%                 disp('Got here 1');
%                 continue;
%             end
%             periodidx = find(spatialPeriods(j) == periods);
%             if (isempty(periodidx))
%                 disp('Got here 2');
%                 %  keyboard
%                 continue;
%             end
%             st = sacstats.starttimes{j};
%             L = logical(sacstats.endtimes{j} >= targon_t(j));
%             st(L) = [];    % Omitting saccades that terminate after targets appear
%             Lsac = (st > stimon_t(j)+offsets(1)) & (st < stimon_t(j)+offsets(2));
%             if (any(Lsac))
%                 if (correct(j))
%                     x_sac(coloridx, periodidx) = x_sac(coloridx,periodidx)+1;
%                 end
%                 n_sac(coloridx,periodidx) = n_sac(coloridx,periodidx)+1;
%             end
%             if (~any(Lsac))
%                 if (correct(j))
%                     x_nosac(coloridx, periodidx) = x_nosac(coloridx,periodidx)+1;
%                 end
%                 n_nosac(coloridx,periodidx) = n_nosac(coloridx,periodidx)+1;
%             end
%         end
%     end
%     
%     % color analysis
%     p1 = sum(x_sac,2)./sum(n_sac,2);
%     p2 = sum(x_nosac,2)./sum(n_nosac,2);
%     p =  (sum(x_sac,2)+sum(x_nosac,2))./(sum(x_sac,2)+sum(x_nosac,2));
%     se = sqrt(p1.*(1-p1)./sum(n_sac,2) + p2.*(1-p2)./sum(n_nosac,2));
%     axes('Position',[1 2.5+3*(a-1) AXESSIZE AXESSIZE]);
%     hold on;
%     bar(p2-p1,'black');
%     for i = 1:3
%         plot([i i],[-se(i) se(i)]+(p2(i)-p1(i)),'k-','Linewidth',2);
%         [h,p] = equalproptest([sum(x_sac(i,:)) sum(x_nosac(i,:))],[sum(n_sac(i,:)) sum(n_nosac(i,:))],alpha);
%         if (h)
%             plot(i,(p2(i)-p1(i))+se(i)+.01,'k*');
%         end
%     end
%     set(gca,'XTick',[1 2 3],'XTickLabel',{'Ach','S','L-M'});
%     set(gca,'Ylim',[-.02 .15]);
%     ylabel('Relative suppression');
%     
%     % Joint analysis
%     sfs = stro.sum.exptParams.pixperdeg ./ periods;
%     axes('Position',[3.5 2.5+3*(a-1) AXESSIZE AXESSIZE]);
%     hold on;
%     p = (x_nosac./n_nosac)-(x_sac./n_sac);
%     bar3(flipud(p));
%     h = get(gca,'Children');
%     set(h(1),'EdgeColor',[1 1 1]);
%     set(gca,'XTick',[1 2 3],'XTickLabel',num2str(sfs',2));
%     set(gca,'YTick',[1 2 3],'YTickLabel',fliplr({'Ach','S','L-M'}));
%     set(gca,'Zlim',[.9*min(p(:)) 1.1*max(p(:))]);
%     set(gca,'Zlim',[-.02 .15]);
%     set(gca,'View',[-53 26]);
%     colormap(flipud(gray));
%     zlabel('Relative suppression');
%     
%     % SF analysis  
%     p1 = sum(x_sac,1)./sum(n_sac,1);
%     p2 = sum(x_nosac,1)./sum(n_nosac,1);
%     p =  (sum(x_sac,1)+sum(x_nosac,1))./(sum(x_sac,1)+sum(x_nosac,1));
%     se = sqrt(p1.*(1-p1)./sum(n_sac,1) + p2.*(1-p2)./sum(n_nosac,1));
%     axes('Position',[6 2.5+3*(a-1) AXESSIZE AXESSIZE]);
%     hold on;
%     bar(p2-p1,'black');
%     for i = 1:3
%         plot([i i],[-se(i) se(i)]+fliplr((p2(i)-p1(i))),'k-','Linewidth',2);
%         [h,p] = equalproptest([sum(x_sac(:,i)) sum(x_nosac(:,i))],[sum(n_sac(:,i)) sum(n_nosac(:,i))],alpha);
%         if (h)
%             plot(i,(p2(i)-p1(i))+se(i)+.01,'k*');
%         end
%     end
%     set(gca,'XTick',[1 2 3],'XTickLabel',num2str(sfs',2));
%     set(gca,'Ylim',[-.02 .15]);
%     ylabel('Relative suppression');
%    
%     
%     % Monkey label
%     axes('Position',[3.5 4.5+3*(a-1) AXESSIZE 1],'Box','off','Visible','off');
%     text(0.5,0,['Monkey ',filenames(1,1)],'HorizontalAlignment','center','FontSize',16)
% end

%%
% Section 5
% Detection performance summary

filenamelist = {'N:\NexFiles\nexfilelists\Greg\DTEM\SednaQuestLrg.txt',...
    'N:\NexFiles\nexfilelists\Greg\DTEM\KaliQuestLrg.txt'};

figure('Units','inches','PaperPosition',[0 0 8.5 11],'Position',[0 0 8.5 11],'PaperOrientation','Portrait','DefaultAxesUnits','inches');
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultTextFontName','Helvetica','DefaultTextFontSize',12);
colordef white;
AXESSIZE = 1.5;
colordirs = mkbasis([1 1 1; 0 0 1; 1 -1 0]');
periods = [206 58 16];
axeshandles = [];
for a = 1:length(filenamelist)
    filenames = fnamesFromTxt(filenamelist{a});
    thresholds = nan*ones(length(colordirs), length(periods), size(filenames,1));
    for b = 1:size(filenames,1)     
        stro = nex2stro(findfile(filenames(b,:)));
  %      stro = DTfilterquesttrials(stro,'PaperDefaults');  Don't filter
  %      trials for this analysis!

        if isempty(stro.sum.analog.sigid) || isempty(stro.trial)
            continue;
        end
        
        [tmpthresholds, tmpcolorDirs, tmpsfs] = DTquestUnpackGH(stro, 'mode');
        %Color is on the rows, sf is on the columns of tmpthresholds.
        sfs = stro.sum.exptParams.pixperdeg ./ periods;
        tmpcolorDirs = mkbasis(tmpcolorDirs')';
        tmpperiods = stro.sum.exptParams.pixperdeg./tmpsfs;
        for i = 1:size(colordirs,1)
            for j = 1:length(periods)
                coloridx = find(colordirs(:,i)'*mkbasis(tmpcolorDirs') > .99);
                periodidx = find(periods(j) == round(tmpperiods));
                if (~isempty(coloridx) & ~isempty(periodidx))
                    thresholds(i,j,b) = tmpthresholds(coloridx,periodidx)/100;
                end
            end
        end
    end
    
    axeshandles(a) = axes('Position',[1+3*(a-1) 2 2 2]); hold on;
    plotcolors = [0 0 0; 0 0 1; 1 0 0]
    set(gca,'ColorOrder',plotcolors)
    mn = nanmean(1./thresholds,3);
    se = nanstd(1./thresholds,[],3)./sqrt(sum(~isnan(thresholds),3));
    plot(sfs,mn','Linewidth',2)
    set(gca,'XScale','log','YScale','log')
    title(['Monkey ',filenames(1,1)],'HorizontalAlignment','center','FontSize',16)
    for i = 1:3
        for j = 1:3
            h = plot([sfs(j) sfs(j)],[mn(i,j)+se(i,j) mn(i,j)-se(i,j)]);
            set(h,'Linewidth',2,'Color',plotcolors(i,:));
        end
    end
end
for a = 1:2
    set(axeshandles(a),'XLim',[.9*min(sfs) 1.1*max(sfs)]);
    set(axeshandles(a),'YLim',[9 500]);
    set(axeshandles(a),'Xticklabel',1,'yticklabel',[10 100]);
    xlabel('SF (cyc/deg)','FontSize',12); ylabel('Sensitivity','FontSize',12);
end

%%
% Section 5.1
% Bitmap of a Gabor

npix = 200;
imagegain = .5;
y = cos(linspace(0,6*pi,npix))';
gaussvect = normpdf(linspace(-npix/2,npix/2,npix),0,npix/4);
gaussmat = gaussvect'*gaussvect;
mat = gaussmat.*repmat(y,1,npix);
mat = mat./max(abs(mat(:)))/2*imagegain+.5;
figure;
image(mat*225)
set(gca,'Xtick',[],'YTick',[],'Box','off');
axis image;
colormap(gray(255));

%%
% Section 6
% Microsaccadic suppression summary.  Both "relative suppression index"
% (i.e. change in proportion correct) and something having to do with
% thresholds.
cd ('C:\Matlab Code\Analysis\Sandbox\Greg\DTemstuff');
filenamelist = {'N:\NexFiles\nexfilelists\Greg\DTEM\SednaQuestLrg.txt',...
    'N:\NexFiles\nexfilelists\Greg\DTEM\KaliQuestLrg.txt'};

figure('Units','inches','PaperPosition',[0 0 8.5 11],'Position',[0 0 8.5 11],'PaperOrientation','Portrait','DefaultAxesUnits','inches');
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultTextFontName','Helvetica','DefaultTextFontSize',12);
colordef white;
AXESSIZE = 1.5;

offsets = [.167 .5];   % relative to stim on, when to look for saccades
colordirs = mkbasis([1 1 1; 1 -1 0]');
periods = [206 58 16];
alpha = 0.05;
for monkidx = 1:length(filenamelist)
    x_sac = zeros(size(colordirs,2), length(periods));
    n_sac = zeros(size(colordirs,2), length(periods));
    x_nosac = zeros(size(colordirs,2), length(periods));
    n_nosac = zeros(size(colordirs,2), length(periods));
    data = cell(length(periods), size(colordirs,2));
    
    filenames = fnamesFromTxt(filenamelist{monkidx});
    for b = 1:size(filenames,1)
        stro = nex2stro(findfile(filenames(b,:)));
        stro = DTfilterquesttrials(stro,'PaperDefaults');

        if isempty(stro.sum.analog.sigid) || isempty(stro.trial)
            continue;
        end
        
        ntrials = size(stro.trial,1);
        lms = stro.sum.exptParams.RF_colors;
        lms = reshape(lms,[3,size(lms,1)/3]);
        lms(:,all(lms == 0)) = [];
        lms = mkbasis(lms);
        whichcol = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'color_dir'));
        stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_on'));
        stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_off'));
        targon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targ_on'));
        spatialPeriods = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'gabor_lambda'));
        correct = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'correct'));
        flashR = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_R'));
        flashG = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_G'));
        flashB = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_B'));
        bkgndrgb = [stro.sum.exptParams.bkgnd_r, stro.sum.exptParams.bkgnd_g, stro.sum.exptParams.bkgnd_b];
        M = reshape(stro.sum.exptParams.m_mtx, 3, 3);
        bkgndlms = M * bkgndrgb';
        x = 0:255; %the normal range of the gamma look up table
        xx = linspace(0, 255, 2^16); %the desired quantization of the gammaTable
        g1 = reshape(stro.sum.exptParams.gamma_table, 256, 3);
        gammaTable = [spline(x, g1(:,1), xx)', spline(x, g1(:,2), xx)', spline(x, g1(:,3), xx)'];
        
        sacstats = getSacData(stro);
        close;
        
        % Looping over trials
        % And getting data for precent correct analysis
        saccadeoccurred = nan*ones(ntrials,1);
        for j = 1:ntrials
            coloridx = find(abs(lms(:,whichcol(j))'*colordirs) > .99);  % Avoiding roundoff errors
            if (isempty(coloridx))
                disp('Got here 1');
                continue;
            end
            periodidx = find(spatialPeriods(j) == periods);
            if (isempty(periodidx))
                disp('Got here 2');
                %  keyboard
                continue;
            end
            st = sacstats.starttimes{j};
            L = logical(sacstats.endtimes{j} >= targon_t(j));
            st(L) = [];    % Omitting saccades that terminate after targets appear
            Lsac = (st > stimon_t(j)+offsets(1)) & (st < stimon_t(j)+offsets(2));
            if (any(Lsac))
                saccadeoccurred(j) = 1;
                if (correct(j))
                    x_sac(coloridx, periodidx) = x_sac(coloridx,periodidx)+1;
                end
                n_sac(coloridx,periodidx) = n_sac(coloridx,periodidx)+1;
            end
            if (~any(Lsac))
                saccadeoccurred(j) = 0;
                if (correct(j))
                    x_nosac(coloridx, periodidx) = x_nosac(coloridx,periodidx)+1;
                end
                n_nosac(coloridx,periodidx) = n_nosac(coloridx,periodidx)+1;
            end
        end
        
        % Stuff for getting data for threshold analysis
        for a = 1:size(colordirs,2)
            for b = 1:length(periods)
                coloridx = find(abs(lms'*colordirs(:,a)) > .99);  % Avoiding roundoff errors
                if (isempty(coloridx))
                    continue;
                end
                L = whichcol == coloridx;
                L = L & spatialPeriods == periods(b);
                L = L & stro.trial(:,strcmp(stro.sum.trialFields(1,:),'numframes')) > 0;
                % Above, need to eliminate 0 frame trials.  What are these?
                RGB =[flashR(L), flashG(L), flashB(L)];
                rgb = [gammaTable(RGB(:,1)+1, 1), gammaTable(RGB(:,2)+1, 2), gammaTable(RGB(:,3)+1, 3)];
                cc = ((M*rgb')'-repmat(bkgndlms',sum(L),1)) ./ repmat(bkgndlms',sum(L),1);
                ccnorms = sqrt(sum(cc.^2,2));
                if (any(ccnorms > 1.5))
                    keyboard
                end
                tmp = [saccadeoccurred(L), correct(L), ccnorms];
                data{b,a} = cat(1,data{b,a},tmp);
            end
        end
    end
    
    % Joint analysis
    sfs = stro.sum.exptParams.pixperdeg ./ periods;
    axes('Position',[5 2.5+3*(monkidx-1) AXESSIZE AXESSIZE]);
    hold on;
    relsup = flipud((x_nosac./n_nosac)-(x_sac./n_sac));
    bar3(relsup);
    h = get(gca,'Children');
    set(h(1),'EdgeColor',[1 1 1]);
    set(gca,'XTick',[1 2 3]-.5,'XTickLabel',num2str(sfs',2));
    set(gca,'YTick',[1 2]-.5,'YTickLabel',fliplr({'Ach','L-M'}));
    set(gca,'Zlim',[.9*min(relsup(:)) 1.1*max(relsup(:))]);
    set(gca,'Ylim',[-1 3],'Xlim',[0 3]);
    set(gca,'View',[-28 24]);
    colormap(flipud(gray));
    zlabel('\Delta Percent Correct');
    for a = 1:size(colordirs,2)
        for b = 1:length(periods)
            h = equalproptest([x_sac(a,b) x_nosac(a,b)],[n_sac(a,b) n_nosac(a,b)],0.05);
            if (h)
                whichcol = -a+3;
                h = plot3(b,whichcol,relsup(whichcol,b),'k*');
                if (b == length(periods))
                    set(h,'MarkerEdgeColor','white');
                end
            end
        end
    end

    % Thresholds
    % Doing a significance test using the asymptotic -2*log likelihood thing
    out = [];
    for a = 1:size(colordirs,2)
        for b = 1:length(periods)
            tmp = data{b,a};
            Lsac = logical(tmp(:,1));
            [fittedparams, success(1)] = weibullFitGH(tmp(:,3), tmp(:,2), 'sse', [mean(tmp(:,3)), 1]);
            [fittedparams_full, success(2),ses] = weibullFitGH([tmp(:,3), Lsac], tmp(:,2), 'mle', [fittedparams(1),fittedparams(2)]);
            [fittedparams_nested, success(2)] = weibullFitGH(tmp(:,3), tmp(:,2), 'mle', [fittedparams(1),fittedparams(2)]);
            pred_nested = 1-.5*exp(-(tmp(:,3)./fittedparams_nested(1)).^fittedparams_nested(2));
            pred_full = 1-.5*exp(-(tmp(:,3)./(fittedparams_full(1)+tmp(:,1)*fittedparams_full(4))).^fittedparams_full(2));
            dev_nested = sum(log(sum([tmp(:,2).*pred_nested ~tmp(:,2).*(1-pred_nested)],2)));
            dev_full = sum(log(sum([tmp(:,2).*pred_full ~tmp(:,2).*(1-pred_full)],2)));
            teststat = -2*(dev_nested-dev_full);
            p(a,b) = 1-chi2cdf(teststat,1);
            out(b,a,1) = fittedparams_full(1)+fittedparams_full(4);  % Threshold with saccade would be high
            out(b,a,2) = fittedparams_full(1);  % Threshold with saccade would be high
            out(b,a,3) = ses(3); % Standard error of beta(1) (aka "delta" in weibullFitGH.m)
        end
    end
    p
    out(:,:,1)-out(:,:,2) % estimate of beta(1)
    out(:,:,3) % se
    
    % Now in "log units"
    raweffectsize = flipud(log10(out(:,:,1)./out(:,:,2))')
    effectsize = min(raweffectsize,0.2);  % upper bound for plotting
    axes('Position',[2 2.5+3*(monkidx-1) AXESSIZE AXESSIZE]);
    hold on;
    bar3(effectsize);
    h = get(gca,'Children');
    set(h(1),'EdgeColor',[1 1 1]);
    set(gca,'XTick',[1 2 3]-.5,'XTickLabel',num2str(sfs',2));
    set(gca,'YTick',[1 2]-.5,'YTickLabel',fliplr({'Ach','L-M'}));
    set(gca,'Zlim',[.9*min(effectsize(:)) 1.1*max(effectsize(:))]);
    if (any(effectsize(:) ~= raweffectsize(:)))
        set(gca,'Ztick',[0 0.05 .1 .15 .2])
        set(gca,'Zticklabel',[0 0.05 .1 .15 round(max(raweffectsize(:))*10)/10]);
    end
    set(gca,'View',[-28 24]);
    colormap(flipud(gray));
    zlabel('\Delta Threshold (log units)');

    for a = 1:size(colordirs,2)
        for b = 1:length(periods)
            if (p(a,b) < 0.05)
                whichcol = -a+3;
                h = plot3(b,whichcol,effectsize(whichcol,b),'k*');
                if (b == length(periods))
                    set(h,'MarkerEdgeColor','white');
                end
            end
        end
    end
    set(gca,'Ylim',[-1 3],'Xlim',[0 3]);
    
    % Monkey label
    axes('Position',[3.5 4.5+3*(monkidx-1) AXESSIZE 1],'Box','off','Visible','off');
    text(0.5,0,['Monkey ',filenames(1,1)],'HorizontalAlignment','center','FontSize',16)
end

%%
% Section 6.1
% Microsaccadic suppression summary.  Change in percent correct and 
% change in detection threshold.  As a bar plot.

cd ('C:\Matlab Code\Analysis\Sandbox\Greg\DTemstuff');
filenamelist = {'N:\NexFiles\nexfilelists\Greg\DTEM\SednaQuestLrg.txt',...
    'N:\NexFiles\nexfilelists\Greg\DTEM\KaliQuestLrg.txt'};

figure('Units','inches','PaperPosition',[0 0 8.5 11],'Position',[0 0 8.5 11],'PaperOrientation','Portrait','DefaultAxesUnits','inches');
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultTextFontName','Helvetica','DefaultTextFontSize',12);
set(gcf,'DefaultAxesColorOrder',[0 0 0; 1 0 0]);
colordef white;
AXESSIZE = 1.5;

offsets = [.167 .5];   % relative to stim on, when to look for saccades
colordirs = mkbasis([1 1 1; 1 -1 0]');
periods = [206 58 16];
alpha = 0.05;
for monkidx = 1:length(filenamelist)
    x_sac = zeros(size(colordirs,2), length(periods));
    n_sac = zeros(size(colordirs,2), length(periods));
    x_nosac = zeros(size(colordirs,2), length(periods));
    n_nosac = zeros(size(colordirs,2), length(periods));
    data = cell(length(periods), size(colordirs,2));
    
    filenames = fnamesFromTxt(filenamelist{monkidx});
    for b = 1:size(filenames,1)
        stro = nex2stro(findfile(filenames(b,:)));
        stro = DTfilterquesttrials(stro,'PaperDefaults');

        if isempty(stro.sum.analog.sigid) || isempty(stro.trial)
            continue;
        end
        
        ntrials = size(stro.trial,1);
        lms = stro.sum.exptParams.RF_colors;
        lms = reshape(lms,[3,size(lms,1)/3]);
        lms(:,all(lms == 0)) = [];
        lms = mkbasis(lms);
        whichcol = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'color_dir'));
        stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_on'));
        stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_off'));
        targon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targ_on'));
        spatialPeriods = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'gabor_lambda'));
        correct = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'correct'));
        flashR = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_R'));
        flashG = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_G'));
        flashB = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_B'));
        bkgndrgb = [stro.sum.exptParams.bkgnd_r, stro.sum.exptParams.bkgnd_g, stro.sum.exptParams.bkgnd_b];
        M = reshape(stro.sum.exptParams.m_mtx, 3, 3);
        bkgndlms = M * bkgndrgb';
        x = 0:255; %the normal range of the gamma look up table
        xx = linspace(0, 255, 2^16); %the desired quantization of the gammaTable
        g1 = reshape(stro.sum.exptParams.gamma_table, 256, 3);
        gammaTable = [spline(x, g1(:,1), xx)', spline(x, g1(:,2), xx)', spline(x, g1(:,3), xx)'];
        
        sacstats = getSacData(stro);
        close;
        
        % Looping over trials
        % And getting data for precent correct analysis
        saccadeoccurred = nan*ones(ntrials,1);
        for j = 1:ntrials
            coloridx = find(abs(lms(:,whichcol(j))'*colordirs) > .99);  % Avoiding roundoff errors
            if (isempty(coloridx))
                disp('Got here 1');
                continue;
            end
            periodidx = find(spatialPeriods(j) == periods);
            if (isempty(periodidx))
                disp('Got here 2');
                %  keyboard
                continue;
            end
            st = sacstats.starttimes{j};
            L = logical(sacstats.endtimes{j} >= targon_t(j));
            st(L) = [];    % Omitting saccades that terminate after targets appear
            Lsac = (st > stimon_t(j)+offsets(1)) & (st < stimon_t(j)+offsets(2));
            if (any(Lsac))
                saccadeoccurred(j) = 1;
                if (correct(j))
                    x_sac(coloridx, periodidx) = x_sac(coloridx,periodidx)+1;
                end
                n_sac(coloridx,periodidx) = n_sac(coloridx,periodidx)+1;
            end
            if (~any(Lsac))
                saccadeoccurred(j) = 0;
                if (correct(j))
                    x_nosac(coloridx, periodidx) = x_nosac(coloridx,periodidx)+1;
                end
                n_nosac(coloridx,periodidx) = n_nosac(coloridx,periodidx)+1;
            end
        end
        
        % Stuff for getting data for threshold analysis
        for a = 1:size(colordirs,2)
            for b = 1:length(periods)
                coloridx = find(abs(lms'*colordirs(:,a)) > .99);  % Avoiding roundoff errors
                if (isempty(coloridx))
                    continue;
                end
                L = whichcol == coloridx;
                L = L & spatialPeriods == periods(b);
                L = L & stro.trial(:,strcmp(stro.sum.trialFields(1,:),'numframes')) > 0;
                % Above, need to eliminate 0 frame trials.  What are these?
                RGB =[flashR(L), flashG(L), flashB(L)];
                rgb = [gammaTable(RGB(:,1)+1, 1), gammaTable(RGB(:,2)+1, 2), gammaTable(RGB(:,3)+1, 3)];
                cc = ((M*rgb')'-repmat(bkgndlms',sum(L),1)) ./ repmat(bkgndlms',sum(L),1);
                ccnorms = sqrt(sum(cc.^2,2));
                if (any(ccnorms > 1.5))
                    keyboard
                end
                tmp = [saccadeoccurred(L), correct(L), ccnorms];
                data{b,a} = cat(1,data{b,a},tmp);
            end
        end
    end

   % Joint analysis
    sfs = stro.sum.exptParams.pixperdeg ./ periods;
    axes('Position',[5 2.5+3*(monkidx-1) AXESSIZE AXESSIZE]);
    hold on;
    relsup = flipud((x_nosac./n_nosac)-(x_sac./n_sac))*100;
    ses = sqrt((x_nosac./n_nosac).*(1-(x_nosac./n_nosac))./n_nosac +...
        (x_sac./n_sac).*(1-(x_sac./n_sac))./n_sac)*100;

    bar3(relsup);
    h = get(gca,'Children');
    set(h(1),'EdgeColor',[1 1 1]);
    set(gca,'XTick',[1 2 3]-.5,'XTickLabel',num2str(sfs',2));
    set(gca,'YTick',[1 2]-.5,'YTickLabel',fliplr({'Ach','L-M'}));
    set(gca,'Zlim',[.9*min(relsup(:)) 1.1*max(relsup(:)+ses(:))]);
    set(gca,'Ylim',[-1 3],'Xlim',[0 3]);
    set(gca,'View',[-28 24]);
    colormap(flipud(gray));
    zlabel('\Delta Percent Correct');
    for a = 1:size(colordirs,2)
        for b = 1:length(periods)
            whichcol = -a+3;
            h1 = plot3([b b],[whichcol whichcol],relsup(whichcol,b)+[0 ses(whichcol,b)],'k-','Linewidth',4,'Color',[.8 .8 .8]);
            h = equalproptest([x_sac(a,b) x_nosac(a,b)],[n_sac(a,b) n_nosac(a,b)],0.05);
            if (h)
                h2 = plot3(b,whichcol,relsup(whichcol,b),'k*');
                if (b == length(periods))
                    set(h2,'MarkerEdgeColor','white');
                end
            end
        end
    end
    

    % Thresholds
    % Doing a significance test using the asymptotic -2*log likelihood thing
    out = [];
    threshses = zeros(size(colordirs,1),length(periods));
    niter = 200;
    for a = 1:size(colordirs,2)
        for b = 1:length(periods)
            statforse = zeros(niter,1);
            for c = 1:niter
                tmp = data{b,a};
                Lsac = logical(tmp(:,1));
                if (c > 1)
                    Lsac = Lsac(randperm(length(Lsac)));
                end
                [fittedparams, success(1)] = weibullFitGH(tmp(:,3), tmp(:,2), 'sse', [mean(tmp(:,3)), 1]);
                [fittedparams_full, success(2), b_ses] = weibullFitGH([tmp(:,3), Lsac], tmp(:,2), 'mle', [fittedparams(1),fittedparams(2)]);
                [fittedparams_nested, success(2)] = weibullFitGH(tmp(:,3), tmp(:,2), 'mle', [fittedparams(1),fittedparams(2)]);
                pred_nested = 1-.5*exp(-(tmp(:,3)./fittedparams_nested(1)).^fittedparams_nested(2));
                pred_full = 1-.5*exp(-(tmp(:,3)./(fittedparams_full(1)+tmp(:,1)*fittedparams_full(4))).^fittedparams_full(2));
                dev_nested = sum(log(sum([tmp(:,2).*pred_nested ~tmp(:,2).*(1-pred_nested)],2)));
                dev_full = sum(log(sum([tmp(:,2).*pred_full ~tmp(:,2).*(1-pred_full)],2)));
                teststat = -2*(dev_nested-dev_full);
                statforse(c) = log10((fittedparams_full(1)+fittedparams_full(4))./fittedparams_full(1)); 
                if (c == 1)
                    teststat = -2*(dev_nested-dev_full);
                    p(a,b) = 1-chi2cdf(teststat,1);
                    out(b,a,1) = fittedparams_full(1)+fittedparams_full(4);  % Threshold with saccade would be high
                    out(b,a,2) = fittedparams_full(1);  % Threshold with saccade would be high
                    out(b,a,3) = b_ses(3); % Standard error of beta(1) (aka "delta" in weibullFitGH.m)
                end
            end
           ses(a,b) = std(statforse)
        end
    end
    out(:,:,1)-out(:,:,2) % estimate of beta(1)
    out(:,:,3) % se

    % Now in "log units"
    raweffectsize = flipud(log10(out(:,:,1)./out(:,:,2))')
    sesflipped = flipud(ses);
  %  effectsize = min(raweffectsize,0.2);  % upper bound for plotting
  effectsize = raweffectsize;
  axes('Position',[2 2.5+3*(monkidx-1) AXESSIZE AXESSIZE]);
    hold on;
    bar3(effectsize);
    h = get(gca,'Children');
    set(h(1),'EdgeColor',[1 1 1]);
    set(gca,'XTick',[1 2 3]-.5,'XTickLabel',num2str(sfs',2));
    set(gca,'YTick',[1 2]-.5,'YTickLabel',fliplr({'Ach','L-M'}));
    set(gca,'Zlim',[.9*min(effectsize(:)) 1.1*max(effectsize(:)+sesflipped(:))]);
    if (any(effectsize(:) ~= raweffectsize(:)))
    %    set(gca,'Ztick',[0 0.05 .1 .15 .2])
     %   set(gca,'Zticklabel',[0 0.05 .1 .15 round(max(raweffectsize(:))*10)/10]);
    end
    set(gca,'View',[-28 24]);
    colormap(flipud(gray));
    zlabel('\Delta Threshold (log units)');

    for a = 1:size(colordirs,2)
        for b = 1:length(periods)
            whichcol = -a+3;
            plot3([b b],[whichcol whichcol],effectsize(whichcol,b)+[0 sesflipped(whichcol,b)],'k-','Color',[.8 .8 .8],'LineWidth',4)
            if (p(a,b) < 0.05)
                h = plot3(b,whichcol,effectsize(whichcol,b),'k*');
                if (b == length(periods))
                    set(h,'MarkerEdgeColor','white');
                end
            end
        end
    end
    set(gca,'Ylim',[-1 3],'Xlim',[0 3]);
    
    % Monkey label
    axes('Position',[3.5 4.5+3*(monkidx-1) AXESSIZE 1],'Box','off','Visible','off');
    text(0.5,0,['Monkey ',filenames(1,1)],'HorizontalAlignment','center','FontSize',16)
end

%%
% Section 6.2
% Microsaccadic suppression summary.  Change in percent correct and 
% change in detection threshold.  As a line plot.

cd ('C:\Matlab Code\Analysis\Sandbox\Greg\DTemstuff');
filenamelist = {'N:\NexFiles\nexfilelists\Greg\DTEM\SednaQuestLrg.txt',...
    'N:\NexFiles\nexfilelists\Greg\DTEM\KaliQuestLrg.txt'};

figure('Units','inches','PaperPosition',[0 0 8.5 11],'Position',[0 0 8.5 11],'PaperOrientation','Portrait','DefaultAxesUnits','inches');
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultTextFontName','Helvetica','DefaultTextFontSize',12);
set(gcf,'DefaultAxesColorOrder',[0 0 0; 1 0 0]);
colordef white;
AXESSIZE = 1.5;

offsets = [.167 .5];   % relative to stim on, when to look for saccades
colordirs = mkbasis([1 1 1; 1 -1 0]');
periods = [206 58 16];
alpha = 0.05;
for monkidx = 1:length(filenamelist)
    x_sac = zeros(size(colordirs,2), length(periods));
    n_sac = zeros(size(colordirs,2), length(periods));
    x_nosac = zeros(size(colordirs,2), length(periods));
    n_nosac = zeros(size(colordirs,2), length(periods));
    data = cell(length(periods), size(colordirs,2));
    
    filenames = fnamesFromTxt(filenamelist{monkidx});
    for b = 1:size(filenames,1)
        stro = nex2stro(findfile(filenames(b,:)));
        stro = DTfilterquesttrials(stro,'PaperDefaults');

        if isempty(stro.sum.analog.sigid) || isempty(stro.trial)
            continue;
        end
        
        ntrials = size(stro.trial,1);
        lms = stro.sum.exptParams.RF_colors;
        lms = reshape(lms,[3,size(lms,1)/3]);
        lms(:,all(lms == 0)) = [];
        lms = mkbasis(lms);
        whichcol = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'color_dir'));
        stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_on'));
        stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_off'));
        targon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targ_on'));
        spatialPeriods = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'gabor_lambda'));
        correct = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'correct'));
        flashR = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_R'));
        flashG = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_G'));
        flashB = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_B'));
        bkgndrgb = [stro.sum.exptParams.bkgnd_r, stro.sum.exptParams.bkgnd_g, stro.sum.exptParams.bkgnd_b];
        M = reshape(stro.sum.exptParams.m_mtx, 3, 3);
        bkgndlms = M * bkgndrgb';
        x = 0:255; %the normal range of the gamma look up table
        xx = linspace(0, 255, 2^16); %the desired quantization of the gammaTable
        g1 = reshape(stro.sum.exptParams.gamma_table, 256, 3);
        gammaTable = [spline(x, g1(:,1), xx)', spline(x, g1(:,2), xx)', spline(x, g1(:,3), xx)'];
        
        sacstats = getSacData(stro);
        close;
        
        % Looping over trials
        % And getting data for precent correct analysis
        saccadeoccurred = nan*ones(ntrials,1);
        for j = 1:ntrials
            coloridx = find(abs(lms(:,whichcol(j))'*colordirs) > .99);  % Avoiding roundoff errors
            if (isempty(coloridx))
                disp('Got here 1');
                continue;
            end
            periodidx = find(spatialPeriods(j) == periods);
            if (isempty(periodidx))
                disp('Got here 2');
                %  keyboard
                continue;
            end
            st = sacstats.starttimes{j};
            L = logical(sacstats.endtimes{j} >= targon_t(j));
            st(L) = [];    % Omitting saccades that terminate after targets appear
            Lsac = (st > stimon_t(j)+offsets(1)) & (st < stimon_t(j)+offsets(2));
            if (any(Lsac))
                saccadeoccurred(j) = 1;
                if (correct(j))
                    x_sac(coloridx, periodidx) = x_sac(coloridx,periodidx)+1;
                end
                n_sac(coloridx,periodidx) = n_sac(coloridx,periodidx)+1;
            end
            if (~any(Lsac))
                saccadeoccurred(j) = 0;
                if (correct(j))
                    x_nosac(coloridx, periodidx) = x_nosac(coloridx,periodidx)+1;
                end
                n_nosac(coloridx,periodidx) = n_nosac(coloridx,periodidx)+1;
            end
        end
        
        % Stuff for getting data for threshold analysis
        for a = 1:size(colordirs,2)
            for b = 1:length(periods)
                coloridx = find(abs(lms'*colordirs(:,a)) > .99);  % Avoiding roundoff errors
                if (isempty(coloridx))
                    continue;
                end
                L = whichcol == coloridx;
                L = L & spatialPeriods == periods(b);
                L = L & stro.trial(:,strcmp(stro.sum.trialFields(1,:),'numframes')) > 0;
                % Above, need to eliminate 0 frame trials.  What are these?
                RGB =[flashR(L), flashG(L), flashB(L)];
                rgb = [gammaTable(RGB(:,1)+1, 1), gammaTable(RGB(:,2)+1, 2), gammaTable(RGB(:,3)+1, 3)];
                cc = ((M*rgb')'-repmat(bkgndlms',sum(L),1)) ./ repmat(bkgndlms',sum(L),1);
                ccnorms = sqrt(sum(cc.^2,2));
                if (any(ccnorms > 1.5))
                    keyboard
                end
                tmp = [saccadeoccurred(L), correct(L), ccnorms];
                data{b,a} = cat(1,data{b,a},tmp);
            end
        end
    end

    colors = [0 0 0; 1 0 0];
    sfs = stro.sum.exptParams.pixperdeg ./ periods;
    axes('Position',[5 2.5+3*(monkidx-1) AXESSIZE AXESSIZE]);
    hold on;
    relsup = ((x_nosac./n_nosac)-(x_sac./n_sac))*100;
    ses = sqrt((x_nosac./n_nosac).*(1-(x_nosac./n_nosac))./n_nosac +...
        (x_sac./n_sac).*(1-(x_sac./n_sac))./n_sac)*100;
    plot(sfs,relsup','linewidth',2);
    set(gca,'XTick',sfs,'XScale','log','XLim',[.9*sfs(1) 1.1*sfs(end)])
    ylabel('\Delta Percent Correct');
    xlabel('Spatial frequency (cpd)');
    for a = 1:size(colordirs,2)
        for b = 1:length(periods)
            plot([sfs(b), sfs(b)],[relsup(a,b)-ses(a,b) relsup(a,b)+ses(a,b)],'Color',colors(a,:),'LineWidth',2)
            h = equalproptest([x_sac(a,b) x_nosac(a,b)],[n_sac(a,b) n_nosac(a,b)],0.05);
            if (h)
                h = text(sfs(b),relsup(a,b)+ses(a,b),'*','HorizontalAlignment','center');
            end
        end
    end

    % Thresholds
    % Doing a significance test using the asymptotic -2*log likelihood thing
    out = [];
    threshses = zeros(size(colordirs,1),length(periods));
    niter = 20;
    for a = 1:size(colordirs,2)
        for b = 1:length(periods)
            statforse = zeros(niter,1);
            for c = 1:niter
                tmp = data{b,a};
                Lsac = logical(tmp(:,1));
                if (c > 1)
                    Lsac = Lsac(randperm(length(Lsac)));
                end
                [fittedparams, success(1)] = weibullFitGH(tmp(:,3), tmp(:,2), 'sse', [mean(tmp(:,3)), 1]);
                [fittedparams_full, success(2), b_ses] = weibullFitGH([tmp(:,3), Lsac], tmp(:,2), 'mle', [fittedparams(1),fittedparams(2)]);
                [fittedparams_nested, success(2)] = weibullFitGH(tmp(:,3), tmp(:,2), 'mle', [fittedparams(1),fittedparams(2)]);
                pred_nested = 1-.5*exp(-(tmp(:,3)./fittedparams_nested(1)).^fittedparams_nested(2));
                pred_full = 1-.5*exp(-(tmp(:,3)./(fittedparams_full(1)+tmp(:,1)*fittedparams_full(4))).^fittedparams_full(2));
                dev_nested = sum(log(sum([tmp(:,2).*pred_nested ~tmp(:,2).*(1-pred_nested)],2)));
                dev_full = sum(log(sum([tmp(:,2).*pred_full ~tmp(:,2).*(1-pred_full)],2)));
                teststat = -2*(dev_nested-dev_full);
                statforse(c) = log10((fittedparams_full(1)+fittedparams_full(4))./fittedparams_full(1)); 
                if (c == 1)
                    teststat = -2*(dev_nested-dev_full);
                    p(a,b) = 1-chi2cdf(teststat,1);
                    out(b,a,1) = fittedparams_full(1)+fittedparams_full(4);  % Threshold with saccade would be high
                    out(b,a,2) = fittedparams_full(1);  % Threshold with saccade would be high
                    out(b,a,3) = b_ses(3); % Standard error of beta(1) (aka "delta" in weibullFitGH.m)
                end
            end
           ses(a,b) = std(statforse)
        end
    end
    out(:,:,1)-out(:,:,2) % estimate of beta(1)
    out(:,:,3) % se
    
    % Now in "log units"
    raweffectsize = log10(out(:,:,1)./out(:,:,2))';
    effectsize = min(raweffectsize,2);  % upper bound for plotting
    axes('Position',[2 2.5+3*(monkidx-1) AXESSIZE AXESSIZE]);
    hold on;
    plot(sfs,effectsize','linewidth',2);
    set(gca,'XTick',sfs,'XScale','log','XLim',[.9*sfs(1) 1.1*sfs(end)])
    ylabel('\Delta Threshold (log units)');
    xlabel('Spatial frequency (cpd)');
    
    for a = 1:size(colordirs,2)
        for b = 1:length(periods)
            plot([sfs(b), sfs(b)],[effectsize(a,b)-ses(a,b) effectsize(a,b)+ses(a,b)],'Color',colors(a,:),'LineWidth',2)
            if (p(a,b) < 0.05)
                h = text(sfs(b),effectsize(a,b)+ses(a,b),'*','HorizontalAlignment','center');
            end
        end
    end
    
    % Monkey label
    axes('Position',[3.5 4.5+3*(monkidx-1) AXESSIZE 1],'Box','off','Visible','off');
    text(0.5,0,['Monkey ',filenames(1,1)],'HorizontalAlignment','center','FontSize',16)
end

%%
% Section 6.3
% Microsaccadic suppression summary.  Change in percent correct and 
% change in detection threshold.  As a 2-D bar plot with SEs.

cd ('C:\Matlab Code\Analysis\Sandbox\Greg\DTemstuff');
filenamelist = {'N:\NexFiles\nexfilelists\Greg\DTEM\SednaQuestLrg.txt',...
    'N:\NexFiles\nexfilelists\Greg\DTEM\KaliQuestLrg.txt'};

figure('Units','inches','PaperPosition',[0 0 8.5 11],'Position',[0 0 8.5 11],'PaperOrientation','Portrait','DefaultAxesUnits','inches');
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultTextFontName','Helvetica','DefaultTextFontSize',12);
set(gcf,'DefaultAxesColorOrder',[0 0 0; 1 0 0]);
colordef white;
AXESSIZE = 1.5;

offsets = [.167 .5];   % relative to stim on, when to look for saccades
colordirs = mkbasis([1 1 1; 1 -1 0]');
periods = [206 58 16];
alpha = 0.05;
niter = 200;
for monkidx = 1:length(filenamelist)
    x_sac = zeros(size(colordirs,2), length(periods));
    n_sac = zeros(size(colordirs,2), length(periods));
    x_nosac = zeros(size(colordirs,2), length(periods));
    n_nosac = zeros(size(colordirs,2), length(periods));
    data = cell(length(periods), size(colordirs,2));
    
    filenames = fnamesFromTxt(filenamelist{monkidx});
    for b = 1:size(filenames,1)
        stro = nex2stro(findfile(filenames(b,:)));
        stro = DTfilterquesttrials(stro,'PaperDefaults');

        if isempty(stro.sum.analog.sigid) || isempty(stro.trial)
            continue;
        end
        
        ntrials = size(stro.trial,1);
        lms = stro.sum.exptParams.RF_colors;
        lms = reshape(lms,[3,size(lms,1)/3]);
        lms(:,all(lms == 0)) = [];
        lms = mkbasis(lms);
        whichcol = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'color_dir'));
        stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_on'));
        stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_off'));
        targon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targ_on'));
        spatialPeriods = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'gabor_lambda'));
        correct = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'correct'));
        flashR = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_R'));
        flashG = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_G'));
        flashB = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_B'));
        bkgndrgb = [stro.sum.exptParams.bkgnd_r, stro.sum.exptParams.bkgnd_g, stro.sum.exptParams.bkgnd_b];
        M = reshape(stro.sum.exptParams.m_mtx, 3, 3);
        bkgndlms = M * bkgndrgb';
        x = 0:255; %the normal range of the gamma look up table
        xx = linspace(0, 255, 2^16); %the desired quantization of the gammaTable
        g1 = reshape(stro.sum.exptParams.gamma_table, 256, 3);
        gammaTable = [spline(x, g1(:,1), xx)', spline(x, g1(:,2), xx)', spline(x, g1(:,3), xx)'];
        
        sacstats = getSacData(stro);
        close;
        
        % Looping over trials
        % And getting data for precent correct analysis
        saccadeoccurred = nan*ones(ntrials,1);
        for j = 1:ntrials
            coloridx = find(abs(lms(:,whichcol(j))'*colordirs) > .99);  % Avoiding roundoff errors
            if (isempty(coloridx))
                disp('Got here 1');
                continue;
            end
            periodidx = find(spatialPeriods(j) == periods);
            if (isempty(periodidx))
                disp('Got here 2');
                %  keyboard
                continue;
            end
            st = sacstats.starttimes{j};
            L = logical(sacstats.endtimes{j} >= targon_t(j));
            st(L) = [];    % Omitting saccades that terminate after targets appear
            Lsac = (st > stimon_t(j)+offsets(1)) & (st < stimon_t(j)+offsets(2));
            if (any(Lsac))
                saccadeoccurred(j) = 1;
                if (correct(j))
                    x_sac(coloridx, periodidx) = x_sac(coloridx,periodidx)+1;
                end
                n_sac(coloridx,periodidx) = n_sac(coloridx,periodidx)+1;
            end
            if (~any(Lsac))
                saccadeoccurred(j) = 0;
                if (correct(j))
                    x_nosac(coloridx, periodidx) = x_nosac(coloridx,periodidx)+1;
                end
                n_nosac(coloridx,periodidx) = n_nosac(coloridx,periodidx)+1;
            end
        end
        
        % Stuff for getting data for threshold analysis
        for a = 1:size(colordirs,2)
            for b = 1:length(periods)
                coloridx = find(abs(lms'*colordirs(:,a)) > .99);  % Avoiding roundoff errors
                if (isempty(coloridx))
                    continue;
                end
                L = whichcol == coloridx;
                L = L & spatialPeriods == periods(b);
                L = L & stro.trial(:,strcmp(stro.sum.trialFields(1,:),'numframes')) > 0;
                % Above, need to eliminate 0 frame trials.  What are these?
                RGB =[flashR(L), flashG(L), flashB(L)];
                rgb = [gammaTable(RGB(:,1)+1, 1), gammaTable(RGB(:,2)+1, 2), gammaTable(RGB(:,3)+1, 3)];
                cc = ((M*rgb')'-repmat(bkgndlms',sum(L),1)) ./ repmat(bkgndlms',sum(L),1);
                ccnorms = sqrt(sum(cc.^2,2));
                if (any(ccnorms > 1.5))
                    keyboard
                end
                tmp = [saccadeoccurred(L), correct(L), ccnorms];
                data{b,a} = cat(1,data{b,a},tmp);
            end
        end
    end
    
    % Joint analysis
    sfs = stro.sum.exptParams.pixperdeg ./ periods;
    axes('Position',[2 2.5+3*(monkidx-1) AXESSIZE AXESSIZE]);

    hold on;
    relsup = ((x_nosac./n_nosac)-(x_sac./n_sac))*100;
    ses = sqrt((x_nosac./n_nosac).*(1-(x_nosac./n_nosac))./n_nosac +...
        (x_sac./n_sac).*(1-(x_sac./n_sac))./n_sac)*100;
    
    h = bar(relsup');
    set(h(1),'FaceColor',[0 0 0]);
    set(h(2),'FaceColor',[1 0 0]);
    set(gca,'XTick',[1 2 3],'XTicklabel',round(sfs*100)/100)
    ylabel('\Delta Percent Correct');
    xlabel('Spatial Frequency (cpd)');
    for a = 1:size(colordirs,2)
        for b = 1:length(periods)
            x = b+(a-1.5)/3.5;
            h1 = plot([x x],relsup(a,b)+[-ses(a,b) ses(a,b)],'k-','Linewidth',2,'Color',[.5 .5 .5]);
            h = equalproptest([x_sac(a,b) x_nosac(a,b)],[n_sac(a,b) n_nosac(a,b)],0.05);
            if (h)
                h2 = plot(x,relsup(a,b)*1.1,'k*');
            end
        end
    end
    

    % Thresholds
    % Doing a significance test using the asymptotic -2*log likelihood thing
    out = [];
    threshses = zeros(size(colordirs,1),length(periods));
    for a = 1:size(colordirs,2)
        for b = 1:length(periods)
            statforse = zeros(niter,1);
            for c = 1:niter
                tmp = data{b,a};
                Lsac = logical(tmp(:,1));
                if (c > 1)
                    Lsac = Lsac(randperm(length(Lsac)));
                end
                [fittedparams, success(1)] = weibullFitGH(tmp(:,3), tmp(:,2), 'sse', [mean(tmp(:,3)), 1]);
                [fittedparams_full, success(2), b_ses] = weibullFitGH([tmp(:,3), Lsac], tmp(:,2), 'mle', [fittedparams(1),fittedparams(2)]);
                [fittedparams_nested, success(2)] = weibullFitGH(tmp(:,3), tmp(:,2), 'mle', [fittedparams(1),fittedparams(2)]);
                pred_nested = 1-.5*exp(-(tmp(:,3)./fittedparams_nested(1)).^fittedparams_nested(2));
                pred_full = 1-.5*exp(-(tmp(:,3)./(fittedparams_full(1)+tmp(:,1)*fittedparams_full(4))).^fittedparams_full(2));
                dev_nested = sum(log(sum([tmp(:,2).*pred_nested ~tmp(:,2).*(1-pred_nested)],2)));
                dev_full = sum(log(sum([tmp(:,2).*pred_full ~tmp(:,2).*(1-pred_full)],2)));
                teststat = -2*(dev_nested-dev_full);
                statforse(c) = log10((fittedparams_full(1)+fittedparams_full(4))./fittedparams_full(1)); 
                if (c == 1)
                    teststat = -2*(dev_nested-dev_full);
                    p(a,b) = 1-chi2cdf(teststat,1);
                    out(b,a,1) = fittedparams_full(1)+fittedparams_full(4);  % Threshold with saccade would be high
                    out(b,a,2) = fittedparams_full(1);  % Threshold with saccade would be high
                    out(b,a,3) = b_ses(3); % Standard error of beta(1) (aka "delta" in weibullFitGH.m)
                end
            end
           ses(a,b) = std(statforse)
        end
    end
    out(:,:,1)-out(:,:,2) % estimate of beta(1)
    out(:,:,3) % se

    % Now in "log units"
    effectsize = log10(out(:,:,1)./out(:,:,2))';
    axes('Position',[5 2.5+3*(monkidx-1) AXESSIZE AXESSIZE]);
    hold on;
    h = bar(effectsize');
    set(h(1),'FaceColor',[0 0 0]);
    set(h(2),'FaceColor',[1 0 0]);
    set(gca,'XTick',[1 2 3],'XTicklabel',round(sfs*100)/100)
    ylabel('\Delta Threshold (log units)');

    for a = 1:size(colordirs,2)
        for b = 1:length(periods)
            x = b+(a-1.5)/3.5;
            h1 = plot([x x],effectsize(a,b)+[-ses(a,b) ses(a,b)],'k-','Linewidth',2,'Color',[.5 .5 .5]);
            if (p(a,b) < 0.05)
                h = plot(x,effectsize(a,b)*1.2,'k*');
            end
        end
    end
    xlabel('Spatial Frequency (cpd)');

    % Monkey label
    axes('Position',[3.5 4.5+3*(monkidx-1) AXESSIZE 1],'Box','off','Visible','off');
    text(0.5,0,['Monkey ',filenames(1,1)],'HorizontalAlignment','center','FontSize',16)
end

%%
% Section 7
% STAs and saccade-triggered PSTHs for a few example neurons
BINWIDTH = .05;
timebins = [-.25:BINWIDTH:.5];
COcells = {{'K031808005'};{'K051909001','K051909003'};{'K032608004'};{'K043008002'};
    {'K061809006'};{'K080808006'};{'K042109005','K042109009','K042109012'};{'K110408003','K110408006'};
    {'K042309004'};{'K111908002'};{'K040108002'};{'K052308001'};
    {'K081208002'};{'K050409004','K050409007'};{'K012309002'};{'K090909001','K090909004','K090909009'}};

COspikeidx = ones(1,length(COcells));
lumcells = {{'K102809001'};{'K102209008'};{'K101509001'};{'K082709001'};{'K080709004'};
    {'K071708001'};{'K040808001'};{'K032608005'};{'K101409001'};{'K082809001'};{'K110708002'};
    {'K082608005'};{'K052108002'};{'K022208008'};{'K050608002'};{'K090108008'}};
lumspikeidx = ones(1,length(lumcells));

AXESSIZE = 0.5;
AXESMARGIN = 0.25;
SCREENWIDTHPIX = 1024;
SCREENWIDTHCM = 35;
MONDISTCM = 100;

for COLOR = 0:1
    
    figure('Units','inches','PaperPosition',[0 0 8.5 11],'Position',[0 0 8.5 11],'PaperOrientation','Portrait','DefaultAxesUnits','inches');
    set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
    set(gcf,'DefaultTextFontName','Helvetica','DefaultTextFontSize',12);
    colordef white;
    
    
    if (COLOR == 1)
        cell_list = COcells;
        spikeIdx = COspikeidx;
    else
        cell_list = lumcells;
        spikeIdx = lumspikeidx;
    end

    for a = 1:length(cell_list);
        stro = {};
        for i = 1:size(cell_list{a},2)
            tmpstro = nex2stro(findfile(char(cell_list{a}(i))));
            if (isempty(stro))
                stro = tmpstro;
            else
                stro = strocat(stro, tmpstro);
            end
        end
        if (isempty(stro.sum.analog.sigid))
            continue;
        end
        maxT = 8;
        framerate = stro.sum.exptParams.framerate;
        nstixperside = stro.sum.exptParams.nstixperside;
        npixperstix = stro.sum.exptParams.npixperstix;
        ntrials = length(stro.sum.absTrialNum);
        stimonidx = find(strcmp(stro.sum.trialFields(1,:),'stim_on'));
        stimoffidx = find(strcmp(stro.sum.trialFields(1,:),'all_off'));
        nframesidx = find(strcmp(stro.sum.trialFields(1,:),'num_frames'));
        noisetypeidx = find(strcmp(stro.sum.trialFields(1,:),'noise_type'));
        sigmaidxs = strmatch('sigma',stro.sum.trialFields(1,:));
        
        L = stro.trial(:,noisetypeidx) == 1;
        stro.ras(~L,:) = [];
        stro.trial(~L,:) = [];
        
        out = getWhtnsStats(stro,maxT,'STCOVmex', {nstixperside^2, 3, maxT}, stro.sum.rasterCells{spikeIdx(a)});
        tmpstro = [];
        STAs = out{1};
        
        energy = sum(STAs.^2);
        whichframe = find(energy == max(energy));
        
        sacstats = getSacData(stro);
        close;
        
        PSTH = zeros(1,length(timebins));
        for j = 1:size(stro.trial,1)
            stimon_t = stro.trial(j,stimonidx);
            numframes = stro.trial(j,nframesidx);
            stimoff_t = stimon_t+numframes/stro.sum.exptParams.framerate;
            st = sacstats.starttimes{j};
            Lsac = (st > stimon_t) & (st < stimoff_t-.1); %.1 to avoid a saccade that leaves the window
            if any(Lsac)
                Lsac(sacstats.amplitudes{j}(Lsac) > 2) = [];
                spiketimes = [];
                for k = find(Lsac')
                    tmp = stro.ras{j,spikeIdx(a)}-sacstats.starttimes{j}(k);
                    spiketimes = [spiketimes; tmp((tmp > timebins(1)-BINWIDTH/2) & (tmp < timebins(end)+BINWIDTH/2))];
                end
                [n,x] = hist(spiketimes, timebins);
                PSTH = PSTH + n./(BINWIDTH*sum(Lsac));
            end
        end
        PSTH = PSTH./sum(L);
        
        % Preparing STA for plotting
        maxSTA = max(STAs(:,whichframe));
        minSTA = min(STAs(:,whichframe));
        potentialnormfactors = [(.5-eps)./(maxSTA-.5); (-.5+eps)./(minSTA-.5)];
        potentialnormfactors(potentialnormfactors < 0) = []; % if min > mu or max < mu
        normfactor = min(potentialnormfactors);
        muvect = reshape(repmat([.5 .5 .5],nstixperside^2,1),nstixperside^2*3,1);
        STA = normfactor*(STAs(:,whichframe)-muvect)+muvect;
        STA = reshape(STA,[nstixperside nstixperside 3]);
        
        % Creating axes for STA
        [col,row] = ind2sub([4 4],a);
        axes('Position',[(AXESSIZE+AXESMARGIN)*col 6-(AXESSIZE+AXESMARGIN)*row AXESSIZE AXESSIZE],'Box','off','Visible','off');
        image(STA);
        set(gca,'XTick',[],'YTick',[]);
        
        % Creating axes for PSTH
        axes('Position',[4+(AXESSIZE+AXESMARGIN)*col 6-(AXESSIZE+AXESMARGIN)*row AXESSIZE AXESSIZE],'Box','off','Visible','off');
        plot(timebins,PSTH,'k-');
        set(gca,'XLim',[min(timebins) max(timebins)],'XTick',[0],'XTickLabel',[]);
        set(gca,'YLim',[0 round(max(PSTH)*1.2)],'YTick',[0 round(max(PSTH)*1.2)],'FontSize',10);
        hold on;
        plot([0 0],get(gca,'Ylim'),'k:');
        drawnow;
    end
    % Drawing scale bars
    pixpercm = SCREENWIDTHPIX/SCREENWIDTHCM;
    theta = atan2(SCREENWIDTHCM/2, MONDISTCM)*180/pi;
    cmperdeg = SCREENWIDTHCM/(2*theta);
    pixperdeg = pixpercm*cmperdeg/2; % /2 because doublewide pixels
    pixperstim = nstixperside*npixperstix;
    stimperdeg = pixperdeg/pixperstim;
    [col,row] = ind2sub([4 4],a);
    axes('Position',[(AXESSIZE+AXESMARGIN)*col 3-(AXESSIZE*AXESMARGIN)*row AXESSIZE AXESSIZE]);
    plot([0 stimperdeg],[0 0],'k-','linewidth',2);
    set(gca,'XLim',[0 1],'Box','off','Visible','off');
    text(0,-1,'1 deg');
    
    axes('Position',[4+(AXESSIZE+AXESMARGIN)*col 3-(AXESSIZE*AXESMARGIN)*row AXESSIZE AXESSIZE]);
    plot([0 0.5],[0 0],'k-','linewidth',2);
    set(gca,'XLim',[min(timebins) max(timebins)],'Box','off','Visible','off');
    text(0,-1,'500 ms');

end

%%
% Section 8
% Average saccade-triggered PSTHs for a population of color-opponent and
% non-opponent neurons.
BINWIDTH = .01;
timebins = [-.25:BINWIDTH:.5];
BASELINETIME = [-.25 0];
DIPTIME = [.05 .1];
LBASELINETIME = timebins >= BASELINETIME(1) & timebins <= BASELINETIME(2);
LDIPTIME = timebins >= DIPTIME(1) & timebins <= DIPTIME(2);

filenamelists = {'N:\NexFiles\nexfilelists\Greg\LvsM.txt',...
    'N:\NexFiles\nexfilelists\Greg\Lum.txt'};
data = [];
whichlist = [];
for a = 1:length(filenamelists)
    [filenames, spikenums] = fnamesFromTxt2(filenamelists{a});
    tmpdata = nan*ones(length(filenames), length(timebins));
    for i = 1:size(filenames,1)
        stro = {};
        for f = 1:size(filenames{i},1)
            tmpstro = nex2stro(findfile(char(filenames{i}(f))));
            if (isempty(stro))
                stro = tmpstro;
            else
                stro = strocat(stro, tmpstro);
            end
        end
        
        noisetypeidx = find(strcmp(stro.sum.trialFields(1,:),'noise_type'));
        stimonidx = find(strcmp(stro.sum.trialFields(1,:),'stim_on'));
        nframesidx = find(strcmp(stro.sum.trialFields(1,:),'num_frames'));
        
        if (isempty(stro.sum.analog.sigid))
            continue;
        end
        sacstats = getSacData(stro);
        close;
        L = logical(stro.trial(:,noisetypeidx) == 1);
        
        PSTH = zeros(1,length(timebins));
        for j = find(L')
            stimon_t = stro.trial(j,stimonidx);
            numframes = stro.trial(j,nframesidx);
            stimoff_t = stimon_t+numframes/stro.sum.exptParams.framerate;
            st = sacstats.starttimes{j};
            Lsac = (st > stimon_t) & (st < stimoff_t-timebins(end));
            if any(Lsac)
                Lsac(sacstats.amplitudes{j}(Lsac) > 2) = [];
                spiketimes = [];
                for k = find(Lsac')
                    tmp = stro.ras{j,spikenums(i)}-sacstats.starttimes{j}(k);
                    spiketimes = [spiketimes; tmp((tmp > timebins(1)-BINWIDTH/2) & (tmp < timebins(end)+BINWIDTH/2))];
                end
                [n,x] = hist(spiketimes, timebins);
                PSTH = PSTH + n./(BINWIDTH*sum(Lsac));
            end
        end
        tmpdata(i,:) = PSTH./sum(L);
    end
    whichlist = [whichlist; repmat(a,length(filenames),1)];
    data = [data; tmpdata]; 
end

figure('Units','inches','PaperPosition',[0 0 8.5 11],'Position',[0 0 8.5 11],'PaperOrientation','Portrait','DefaultAxesUnits','inches');
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultTextFontName','Helvetica','DefaultTextFontSize',12);
colordef white;

% inserting a yellow line
idx = find(whichlist == 1,1,'last');
PSTHs = data(1:idx,:)./repmat(max(data(1:idx,:),[],2),1,size(data,2));
%PSTHs = data(1:idx,:)./repmat(mean(data(1:idx,:),2),1,size(data,2));

PSTHs = [PSTHs; 100*ones(2,size(data,2))];
PSTHs = [PSTHs; data(idx+1:end,:)./repmat(max(data(idx+1:end,:),[],2),1,size(data,2))];
%PSTHs = [PSTHs; data(idx+1:end,:)./repmat(mean(data(idx+1:end,:),2),1,size(data,2))];

axes('Position',[2 5 4 2]); hold on;
image(PSTHs*(99/max(PSTHs(PSTHs < 100))));
colormap([gray(100); 1 1 0 ])
axis fill;
set(gca,'Ytick',[],'XTick',[]);
% y = [timebins(1):.25:timebins(end)];
% x = (y-timebins(1))/BINWIDTH+1;
% set(gca,'XTick',x);
% set(gca,'XTickLabel',y);
text(-.3,.4,'R-G','rotation',90,'HorizontalAlignment','center','units','inches')
text(-.3,1.5,'Ach','rotation',90,'HorizontalAlignment','center','units','inches')
set(gca,'XLim',[.5 size(data,2)],'YLim',[1 size(data,1)])
axes('Position',[2 4.8 4 .2]); hold on;
plot(timebins(LBASELINETIME),zeros(1,sum(LBASELINETIME)),'k-','LineWidth',5);
plot(timebins(LDIPTIME),zeros(1,sum(LDIPTIME)),'k-','LineWidth',5);

set(gca,'XLim',[timebins(1) timebins(end)],'Visible','off')


normdata = data./repmat(max(data,[],2),1,size(data,2));
%normdata = data./repmat(mean(data,2),1,size(data,2));

axes('Position',[2 2 4 2]); hold on;
L = whichlist == 2; % Luminance cells
mn = mean(normdata(L,:));
se = std(normdata(L,:))/sqrt(sum(L));
h = patch([timebins fliplr(timebins)],[mn fliplr(mn)]+[-se se],[.5 .5 .5]);
plot(timebins,mn,'k-','Linewidth',2);
%set(h,'FaceAlpha',.5);

L = whichlist == 1; % Color opponent cells
mn = mean(normdata(L,:));
se = std(normdata(L,:))/sqrt(sum(L));
h = patch([timebins fliplr(timebins)],[mn fliplr(mn)]+[-se se],[1 .5 .5]);
plot(timebins,mn,'k-','Linewidth',2,'Color','red');
set(gca,'XLim',[min(timebins) max(timebins)]);
xlabel('Time (sec)'); ylabel('sp/sec');
legend({'Lum','L-M'});
set(gca,'YScale','linear')

% A little analysis
sacmodidx = sum(normdata(:,LDIPTIME),2);
[h,p] = ttest2(sacmodidx(whichlist == 1), sacmodidx(whichlist == 2))
sacmodidx = sum(normdata(:,timebins>=.12 & timebins<.3),2);
[h,p] = ttest2(sacmodidx(whichlist == 1), sacmodidx(whichlist == 2))

% Scatterplots for Reviewer 1
figure('Units','inches','PaperPosition',[0 0 8.5 11],'Position',[0 0 8.5 11],'PaperOrientation','Portrait','DefaultAxesUnits','inches');
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultTextFontName','Helvetica','DefaultTextFontSize',12);
colordef white;
axes('Position',[2 5 4 1]); hold on;
%plot(mean(data(whichlist == 1,timebins<0),2), mean(data(whichlist == 1,timebins>=.05 & timebins<.1),2),'k.')
%plot(mean(data(whichlist == 2,timebins<0),2), mean(data(whichlist == 2,timebins>=.05 & timebins<.1),2),'r.')
plot(mean(data(whichlist == 1,LBASELINETIME),2), mean(data(whichlist == 1,LDIPTIME),2)./mean(data(whichlist == 1,LBASELINETIME),2),'ro','MarkerFaceColor','red','MarkerSize',5);
plot(mean(data(whichlist == 2,LBASELINETIME),2), mean(data(whichlist == 2,LDIPTIME),2)./mean(data(whichlist == 2,LBASELINETIME),2),'ko','MarkerFaceColor','black','MarkerSize',5);
xlabel('Baseline (sp/sec)');
ylabel('Fractional suppression of firing rate');
set(gca,'XScale','log','YLim',[0 2.2]);

% Estimated percent change in response
mn1 = geomean(mean(data(:,LDIPTIME),2)./mean(data(:,LBASELINETIME),2))
med = median(mean(data(:,LBASELINETIME),2))
L = mean(data(:,LBASELINETIME),2) > med;
mn2 = geomean(mean(data(L,LDIPTIME),2)./mean(data(L,LBASELINETIME),2))
mn3 = geomean(mean(data(~L,LDIPTIME),2)./mean(data(~L,LBASELINETIME),2))
plot([.1 60],[mn1 mn1],'k:','LineWidth',1);
plot([.1 60],[1 1],'k-','LineWidth',1);


% Comparing high firing cells and low firing cells
[h,p] = ttest2(mean(data(L,LDIPTIME),2)./mean(data(L,LBASELINETIME),2), mean(data(~L,LDIPTIME),2)./mean(data(~L,LBASELINETIME),2))
% One sample t-test
[h,p] = ttest(mean(data(:,LDIPTIME),2)./mean(data(:,LBASELINETIME),2)-1)
[b,bint] = regress(mean(data(:,LDIPTIME),2)./mean(data(:,LBASELINETIME),2), [mean(data(:,LBASELINETIME),2) ones(size(data,1),1)],.05)
[rho,p] = corr(mean(data(:,LDIPTIME),2)./mean(data(:,LBASELINETIME),2), mean(data(:,LBASELINETIME),2))
% Estimated percent change in response separately for Lum and L-M cells
L = mean(data(:,LDIPTIME),2)./mean(data(:, LBASELINETIME),2) > 1.5;  % Option to take out "outliers"
for i = 1:2 % 1 = L-M, 2 = Lum
    mn = geomean(mean(data(~L & whichlist == i,LDIPTIME),2)./mean(data(~L & whichlist == i, LBASELINETIME),2));
    1-mn
    [h,p] = ttest(mean(data(~L & whichlist == i,LDIPTIME),2)./mean(data(~L & whichlist == i, LBASELINETIME),2)-1)
end
% ttest2: is suppression the same for Lum and L-M cells?
[h,p] = ttest2(mean(data(~L & whichlist == 1,LDIPTIME),2)./mean(data(~L & whichlist == 1,LBASELINETIME),2),...
    mean(data(~L & whichlist == 2,LDIPTIME),2)./mean(data(~L & whichlist == 2,LBASELINETIME),2))

%ylabel('Dip (sp/sec)');
%set(gca,'XScale','log','YScale','log')

%[b,bint] = regress(log(mean(data(:,timebins>=.05 & timebins<.1),2)), [log(mean(data(:,timebins<0),2)), ones(size(data,1),1)])
%[b,bint] = regress(mean(data(:,timebins>=.05 & timebins<.1),2), [mean(data(:,timebins<0),2), ones(size(data,1),1)]);

%plot([.1 100],[.1 100],'k-');
%lsendpoints = b(1)*[.1 100]+b(2);
%plot([(.1-b(2))/b(1) 100],[.1 lsendpoints(2)],'k:');
%plot([(.1-b(2))/b(1) 100],[.1 lsendpoints(2)],'k:');

%axes('Position',[5 2 2 2]); hold on;
%plot(mean(data(whichlist == 1,timebins<0),2), mean(data(whichlist == 1,timebins>=.12 & timebins<.3),2),'k.')
%plot(mean(data(whichlist == 2,timebins<0),2), mean(data(whichlist == 2,timebins>=.12 & timebins<.3),2),'r.')
%xlabel('Baseline (sp/sec)');
%ylabel('Rebound (sp/sec)');
%set(gca,'XScale','log','YScale','log')
%plot([.1 100],[.1 100],'k-');



%%
% Section 8.1
% Saccade-triggered PSTHs in DTspot for L-M stimulation (only cells that
% were reasonably responsive to this).

BINWIDTH = .01;
timebins = [-.25:BINWIDTH:.5];

filenamelist = 'N:\NexFiles\nexfilelists\Greg\DTEM\DTspotLvsM.txt';
[filenames, spikenums] = fnamesFromTxt2(filenamelist);
data = nan*ones(length(filenames), length(timebins));
for i = 1:size(filenames,1)
    stro = {};
    for f = 1:size(filenames{i},1)
        tmpstro = nex2stro(findfile(char(filenames{i}(f))));
        if (isempty(stro))
            stro = tmpstro;
        else
            stro = strocat(stro, tmpstro);
        end
    end
    
    if (isempty(stro.sum.analog.sigid))
        continue;
    end
    sacstats = getSacData(stro);
    close;
    
    whichColor = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'color_dir'));
    colors = mkbasis(reshape(stro.sum.exptParams.RF_colors,[3 3]));
    Lcol = whichColor == find(abs([1/sqrt(2) -1/sqrt(2) 0]*colors)>.99);
    
    flashside = sign(stro.sum.exptParams.rf_x) == sign(stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_x')));
    % flashside == 1 means "in RF"
    stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_on'));
    stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_off'));
    targon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targ_on'));
    PSTH = zeros(1,length(timebins));
    L = Lcol;
    for j = find(L')
        st = sacstats.starttimes{j};
        Lsac = (st > stimon_t(j)-.5) & (st < targon_t(j));
        if any(Lsac)
            spiketimes = [];
            for k = find(Lsac')
                tmp = stro.ras{j,1}-sacstats.starttimes{j}(k);
                spiketimes = [spiketimes; tmp((tmp > timebins(1)-BINWIDTH/2) & (tmp < timebins(end)+BINWIDTH/2))];
            end
            [n,x] = hist(spiketimes, timebins);
            PSTH = PSTH + n./(BINWIDTH*sum(Lsac));
        end
    end
    data(i,:) = PSTH./sum(L);
end
PSTHs = data./repmat(max(data,[],2),1,size(data,2));
%PSTHs = data./repmat(mean(data,2),1,size(data,2));
se = std(PSTHs)./sqrt(size(PSTHs,1))
figure('Units','inches','PaperPosition',[0 0 8.5 11],'Position',[0 0 8.5 11],'PaperOrientation','Portrait','DefaultAxesUnits','inches');
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultTextFontName','Helvetica','DefaultTextFontSize',12);
colordef white;
axes('Position',[2 2 4 2]); hold on;
h = patch([timebins fliplr(timebins)],[mean(PSTHs)+se fliplr(mean(PSTHs)-se)]','r');
set(h,'FaceColor',[1 .5 .5]);
plot(timebins,mean(PSTHs),'r.-','Linewidth',2);
set(gca,'Ylim',[0.15 .55]);
set(gca,'Xlim',[timebins(1) timebins(end)])
set(gca,'XLim',[min(timebins) max(timebins)]);
xlabel('Time (sec)'); ylabel('sp/sec');
%set(gca,'YScale','log');

%%
% Section 8.2 CHARLIE'S VERSION OF GREG'S CODE FROM SECTION 8.1
% Saccade-triggered PSTHs in DTspot for L-M stimulation (only cells that
% were reasonably responsive to this).

NSHUFFS = 200;
BINWIDTH = .020;
preTime = 0.250;
postTime = 0.500;
timebins = [-preTime:BINWIDTH:postTime];

%filenamelist = 'N:\NexFiles\nexfilelists\Greg\DTEM\DTspotLvsM.txt'; %the old list
filenamelist = 'N:\NexFiles\nexfilelists\Greg\DTEM\DTspotLvsM_charlie.txt'; %the new list

[filenames, spikenums] = fnamesFromTxt2(filenamelist);
data = nan(length(filenames), length(timebins));
shuffData = nan(length(filenames), length(timebins));
shuffNormToMaxData = nan(length(filenames), length(timebins));
totalNumExpts = size(filenames,1)
for i = 1:size(filenames,1)
    fprintf('computing file <%d>\n', i);
    stro = dtobj(filenames{i}{1});
    stro = stripOutGratingTrials(stro);
    if (isempty(stro.sum.analog.sigid))
        fprintf('Expt <%s> has no analog data', filenames{i}{1});
        continue;
    end
    
    
    whichColor = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'color_dir'));
    colors = mkbasis(reshape(stro.sum.exptParams.RF_colors,[3 3]));
    Lcol = whichColor == find(abs([1/sqrt(2) -1/sqrt(2) 0]*colors)>.99);
    stro.trial(~Lcol,:) = [];
    stro.ras(~Lcol,:) = [];
    
    
    %get the saccade data.
    sacstats = getSacData(stro);
    close;
    
    %setup the suffeled data correspondences
    shuffInd = nan(NSHUFFS+1,size(stro.trial,1));
    shuffInd(1,:) = 1:size(stro.trial,1); % First row is unshuffled
    for a = 2:NSHUFFS+1;
        shuffInd(a,:) = randperm(size(stro.trial,1));
    end
    
    frameon_t = stro.trial(:,stro.idx.frameOn);
    targon_t = stro.trial(:, stro.idx.targOn);
    choice_t = stro.trial(:,stro.idx.choiceTime);
    trl_psth = nan(size(stro.trial,1), length(timebins));
    expt_psth = nan(NSHUFFS+1, length(timebins));
    for shuff = 1:NSHUFFS+1;
        for j = 1:size(stro.trial,1);
            anlyStart = frameon_t(j);
            anlyEnd = frameon_t(j)+0.200+0.666;
            
            sacStarts = sacstats.starttimes{shuffInd(shuff,j)};
            sacStarts = sacStarts - stro.ras{shuffInd(shuff,j), stro.idx.anlgStart}; %sacStarts from time zero
            sacStarts = sacStarts + stro.ras{j, stro.idx.anlgStart}; %sacStarts from trial(trl) anlgStart
            sacStarts(sacStarts < anlyStart) = [];
            sacStarts(sacStarts > anlyEnd) = [];
            if any(sacStarts)
                sac_psth = [];
                for k = 1:length(sacStarts)
                    spiketimes = stro.ras{j,1}-sacStarts(k);
                    spiketimes(spiketimes < timebins(1)-BINWIDTH/2) = [];
                    spiketimes(spiketimes > timebins(end)+BINWIDTH/2) = [];
                    sac_psth(end+1,:) = hist(spiketimes, timebins);
                    sac_psth(end,timebins+sacStarts(k) > choice_t(j)) = nan;
                end
                trl_psth(j,:) = nanmean(sac_psth,1);
            end
        end
        expt_psth(shuff,:) = nanmean(trl_psth,1)./BINWIDTH;
    end
    data(i,:) = expt_psth(1,:);
    tmp = expt_psth(2:end,:);
    shuffData(i,:) = nanmean(tmp,1);
    shuffNormToMaxData(i,:) = nanmean(tmp./repmat(max(tmp,[],2),1,size(tmp,2)),1);
end
normPsth = data./repmat(max(data,[],2),1,size(data,2));
figure; hold on;
plot(timebins,mean(shuffNormToMaxData))
plot(timebins,mean(normPsth))

% Mean PSTH across monkeys.
NORMMETH = 'baseline';
switch lower(NORMMETH)
    case 'max'
        normPsth = data./repmat(max(data,[],2),1,size(data,2));
        shuffNormPsth = shuffNormToMaxData; %normalized to max during for loop
    case 'baseline'
        preBins = timebins<0;
        normPsth = data./repmat(mean(data(:,preBins),2),1,size(data,2));
        shuffNormPsth = shuffData./repmat(mean(shuffData(:,preBins),2),1,size(shuffData,2));
end


%normalized and averaged across monkeys
figure
hold on,
sem = std(normPsth,[],1)./sqrt(size(normPsth,1));
shuffSem = std(shuffNormPsth,[],1)./sqrt(size(shuffNormPsth,1));
patch([timebins, fliplr(timebins)], [mean(shuffNormPsth,1)+shuffSem, fliplr(mean(shuffNormPsth,1)-shuffSem)], 'r', 'facecolor', [1, 0.5, 0.5]);
plot(timebins, mean(shuffNormPsth,1), 'r-', 'Linewidth', 2);
patch([timebins fliplr(timebins)],[mean(normPsth,1)+sem fliplr(mean(normPsth,1)-sem)],'k', 'FaceColor',[.6 .6 .6]);
plot(timebins, mean(normPsth,1), 'k-', 'Linewidth', 2);
xlim([min(timebins) max(timebins)])
xlabel('Time (sec)');
ylabel('Normalized Rate');
legend('Shuff Data', '', 'Non-Shuff Data')
title('Average Across Monkeys')
hold off

% Now the average PSTH for each monkey seperately.
monkNames = char([filenames{:}]);
l_sedna = monkNames(:,1) == 'S';
l_kali = monkNames(:,1) == 'K';
figure
subplot(1,2,1)%sedna
hold on,
s_psth = normPsth(l_sedna,:);
s_avg = mean(s_psth,1);
s_sem = std(s_psth,[],1)./sqrt(sum(l_sedna));
s_shuffPsth = shuffNormPsth(l_sedna,:);
s_shuffAvg = mean(s_shuffPsth,1);
s_shuffSem = std(s_shuffPsth,[],1)./sqrt(size(s_shuffPsth,1));
patch([timebins, fliplr(timebins)], [s_shuffAvg+s_shuffSem, fliplr(s_shuffAvg-s_shuffSem)], 'r', 'facecolor', [1, 0.5, 0.5])
plot(timebins, s_shuffAvg, 'r', 'linewidth', 2)
patch([timebins, fliplr(timebins)], [s_avg+s_sem, fliplr(s_avg-s_sem)], 'k', 'facecolor', [0.6, 0.6, 0.6])
plot(timebins, s_avg, 'k', 'linewidth', 2)
legend('Shuff Data', '', 'Non-Shuff Data')
xlim([min(timebins), max(timebins)])
title('Sedna')
xlabel('time')
ylabel('Normalized Rate')
hold off
subplot(1,2,2)%kali
hold on,
k_psth = normPsth(l_kali,:);
k_avg = mean(k_psth,1);
k_sem = std(k_psth,[],1)./sqrt(sum(l_kali));
k_shuffPsth = shuffNormPsth(l_kali,:);
k_shuffAvg = mean(k_shuffPsth,1);
k_shuffSem = std(k_shuffPsth,[],1)./sqrt(size(k_shuffPsth,1));
patch([timebins, fliplr(timebins)], [k_shuffAvg+k_shuffSem, fliplr(k_shuffAvg-k_shuffSem)], 'r', 'facecolor', [1, 0.5, 0.5])
plot(timebins, k_shuffAvg, 'r', 'linewidth', 2)
patch([timebins, fliplr(timebins)], [k_avg+k_sem, fliplr(k_avg-k_sem)], 'k', 'facecolor', [0.6, 0.6, 0.6])
plot(timebins, k_avg, 'k', 'linewidth', 2)
legend('Shuff Data', '', 'Non-Shuff Data')
xlim([min(timebins), max(timebins)])
title('Kali')
xlabel('time')
ylabel('Normalized Rate')
hold off


%% 
%Section 8.2
% Charlie's attempt at demonstrating reductions in firing rate for sac
% trials over no-sac trials during DT


filenamelist = '/Users/charliehass/LabStuff/NexFiles/nexfilelists/Greg/DTEM/DTspotLvsM_charlie.txt'; %the new list
[filenames, spikenums] = fnamesFromTxt2(filenamelist);
sacZscores = [];
nosacZscores = [];
junkInd = [];
for a = 1:length(filenames);
    %open a file
    filenames{a}{1}
    DT = dtobj(filenames{a}{1});
    DT = stripOutGratingTrials(DT);
    [monk, cell, expt] = DTunpack(DT,1);
    close all
    l_inRF = DT.trial(:, DT.idx.flashX) == DT.sum.exptParams.rf_x;
    contrast = DT.trial(:, DT.idx.cntrstLev);
    colorDir = DT.trial(:, DT.idx.colorDir);
    gaborOn_t = DT.trial(:, DT.idx.flashOn);
    gaborOff_t = DT.trial(:,DT.idx.flashOff);
    spikes = DT.ras(:, DT.idx.spikes);
    tStart = mat2cell(gaborOn_t, ones(length(gaborOn_t),1));
    tEnd = mat2cell(gaborOff_t, ones(length(gaborOff_t),1));
    nSpikes = cellfun(@(x,y,z)(sum((x>y)&(x<=z))), spikes, tStart, tEnd);
    
    %get the analog data. if it's not there than skip this file.
    if ~isempty(DT.sum.analog.sigid)
        sacstats = getSacData(DT);
        nSacDuringGabor = cellfun(@(x,y,z)(sum((x>y)&(x<=z))), sacstats.starttimes', tStart, tEnd);
        sacTrials = nSacDuringGabor>0;
        close(gcf);
    else
        fprintf('   ****** File <%d> has no analog data *****', a)
        continue
    end
    
    %determine which color dir corresponds to L-M.
    testedColors = reshape(DT.sum.exptParams.RF_colors, 3, 3)';
    lvmIdx1 = ismember(sign(testedColors), [1 -1 0], 'rows');
    lvmIdx2 = ismember(sign(testedColors), [-1 1 0], 'rows');
    lvmIdx = find(lvmIdx1 | lvmIdx2);
    
    for clr = 1:2
        if isnan(cell.alpha(clr,1))
            disp('nan neurofun')
            continue
        end
        
        
        for cntrst = 2:max(contrast);
            if (clr==1)&&(cntrst==1) %i.e., just do this once
                sacTList = (~l_inRF | colorDir==1) & sacTrials;
                noSacTList = (~l_inRF | colorDir==1) & ~sacTrials;
            else
                sacTList = l_inRF & (colorDir==clr) & (contrast==cntrst) & sacTrials;
                noSacTList = l_inRF & (colorDir==clr) & (contrast==cntrst) & ~sacTrials;
            end
            sacCounts = nSpikes(sacTList);
            noSacCounts = nSpikes(noSacTList);
            
            %create a pooled estimate of sigma and mu.
            nSac = numel(sacCounts);
            nNoSac = numel(noSacCounts);
            if ~nSac
                sigma = std(noSacCounts);
                mu = mean(noSacCounts);
            elseif ~nNoSac
                sigma = std(sacCounts);
                mu = mean(sacCounts);
            else
                sigma = sqrt(((nSac-1)*var(sacCounts) + (nNoSac-1)*var(noSacCounts))/(nSac+nNoSac-2));
                mu = (nSac.*mean(sacCounts)+nNoSac.*mean(noSacCounts))./ (nSac+nNoSac);
            end
            tmpSacZ = (sacCounts-mu) ./ (sigma+eps);
            tmpNoSacZ = (noSacCounts-mu) ./ (sigma+eps);
            if any([tmpSacZ(:);tmpNoSacZ(:)]>50)
                junkInd(end+1,:)=[a, lvmIdx, cntrst]
                continue
            end
            sacZscores = [sacZscores(:);tmpSacZ(:)];
            nosacZscores = [nosacZscores(:);tmpNoSacZ(:)];
        end
    end
end

%a simple t test on the z scores
[h,p] = ttest2(sacZscores, nosacZscores)

%a permutation test
NPERMS = 10000;
permScores = nan(1,NPERMS);
combDat = [sacZscores(:); nosacZscores(:)];
for a = 1:NPERMS
    permInd = randperm(length(combDat));
    tmpSac = combDat(permInd(1:length(sacZscores)));
    tmpNoSac = combDat(permInd(length(sacZscores)+1:end));
    sigma = sqrt(((length(tmpSac)-1)*var(tmpSac) + (length(tmpNoSac)-1)*var(tmpNoSac))/(length(tmpSac)+length(tmpNoSac)-2));
    permScores(a) = (mean(tmpNoSac)-mean(tmpSac))./sigma;
end
exptPoolSigma = sqrt(((length(sacZscores)-1)*var(sacZscores) + (length(nosacZscores)-1)*var(nosacZscores))/(length(sacZscores)+length(nosacZscores)-2));
exptScore = (mean(nosacZscores)-mean(sacZscores))./exptPoolSigma
critVals = prctile(permScores, [0.025 0.975])
p = sum((permScores<-abs(exptScore)) | (permScores>abs(exptScore)))./length(permScores)

%%
% Section 9
% Plot showing how percent correct varies with microsaccade amplitude.
% OUTDATED - Delete this once analyses are moved over

filenamelist = {'N:\NexFiles\nexfilelists\Greg\DTEM\KaliQuestLrg.txt',...
    'N:\NexFiles\nexfilelists\Greg\DTEM\SednaQuestLrg.txt'};

offsets = [.167 .5];   % relative to stim on, when to look for saccades
%offsets = [.666 1.5];   % relative to stim on, when to look for saccades

BINPRCTLS = [0:8:100];
data = zeros(3,length(BINPRCTLS)-1,2);
nosacdata = zeros(2,2);
for a = 1:length(filenamelist)
    tmpdata = [];
    [filenames, spikenums] = fnamesFromTxt2(filenamelist{a});  
    for b = 1:size(filenames,1)
        stro = nex2stro(findfile(filenames{b}));
        stro = DTfilterquesttrials(stro,'PaperDefaults');
        
        if isempty(stro.sum.analog.sigid) || isempty(stro.trial)
            continue;
        end
        
        ntrials = size(stro.trial,1);
        lms = stro.sum.exptParams.RF_colors;
        lms = reshape(lms,[3,size(lms,1)/3]);
        lms(:,all(lms == 0)) = [];
        lms = mkbasis(lms);
        whichcol = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'color_dir'));
        stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_on'));
        stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_off'));
        targon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targ_on'));
        correct = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'correct'));
        
        sacstats = getSacData(stro);
        close;
        
        % Looping over trials
        % And getting data for precent correct analysis
        for j = 1:ntrials
            coloridx = find(abs(lms(:,whichcol(j))'*[.5774 .5774 .5774]') > .99);  % Avoiding roundoff errors
            if (coloridx ~= 1)
                continue;
            end
            st = sacstats.starttimes{j};
            L = logical(sacstats.endtimes{j} >= targon_t(j)-.2);
            st(L) = [];    % Omitting saccades that terminate after targets appear
            Lsac = (st > stimon_t(j)+offsets(1)) & (st < stimon_t(j)+offsets(2));
            if (any(Lsac))
                k = find(Lsac);
                tmpdata = [tmpdata; correct(j) max(sacstats.amplitudes{j}(k))];
            else
                tmpdata = [tmpdata; correct(j) nan];
            end
        end
    end
    

    binedges = prctile(tmpdata(:,2),BINPRCTLS);
    binedges = linspace(.9*nanmin(tmpdata(:,2)),1.1*nanmax(tmpdata(:,2)),length(BINPRCTLS));
    for j = 2:length(binedges)
        L = tmpdata(:,2)>binedges(j-1) & tmpdata(:,2)<binedges(j);
        binmiddle = (binedges(j)+binedges(j-1))/2;
        data(:,j-1,a) = [binmiddle sum(tmpdata(L,1)) sum(L)];
    end
 
    nosacdata(1,a) = sum(tmpdata(Lnosac,1));
    nosacdata(2,a) = sum(Lnosac);
end

figure('Units','inches','PaperPosition',[0 0 8.5 11],'Position',[0 0 8.5 11],'PaperOrientation','Portrait','DefaultAxesUnits','inches');
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultTextFontName','Helvetica','DefaultTextFontSize',12);
colordef white;
clear h;
h(1) = axes('Position',[2 5 4 2]); hold on;
h(2) = axes('Position',[2 2 4 2]); hold on;

for a = 1:length(filenamelist)
    axes(h(a));
    % getting rid of bins with too few data points
    L = data(3,:,a) < 10
    processeddata = data(:,~L,a);
    
    pcor = processeddata(2,:)./processeddata(3,:);
    patch([processeddata(1,:) fliplr(processeddata(1,:))],...
        [pcor + sqrt((pcor.*(1-pcor))./processeddata(3,:))...
        fliplr(pcor - sqrt((pcor.*(1-pcor))./processeddata(3,:)))],[.5 .5 .5]);
        plot(processeddata(1,:),pcor,'k.-','LineWidth',2)

    pcornosac = sum(nosacdata(1,a))/sum(nosacdata(2,a));
    plot([0 .85],[pcornosac pcornosac],'k:','LineWidth',1);
    set(gca,'XLim',[0 .85]);
    ylabel('P(correct)');
    xlabel('Amplitude (deg)');
    if (a == 1)
        title('Monkey K','FontSize',16);
    else
        title('Monkey S','FontSize',16);
    end
    % Linear regression on binned data
end

%%
% Section 10
% Percent correct as a function of microsaccade amplitude.  Also
% percent correct as a function of time in the trial separately for large
% and small amplitude saccades.  Achromatic only.

figure('Units','inches','PaperPosition',[0 0 8.5 11],'Position',[0 0 8.5 11],'PaperOrientation','Portrait','DefaultAxesUnits','inches');
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultTextFontName','Helvetica','DefaultTextFontSize',12);
colordef white;

AXESSIZE = 1.7;
filenamelist = {'N:\NexFiles\nexfilelists\Greg\DTEM\KaliQuestLrg.txt',...
    'N:\NexFiles\nexfilelists\Greg\DTEM\SednaQuestLrg.txt'};
timebinwidth = .05; %sec
timebins = [-.2:timebinwidth:.767];

offsets = [.167 .5];   % relative to stim on, when to look for saccades
NSACBINEDGES = 11;
sacdata = zeros(3,NSACBINEDGES,2);
nosacdata = zeros(2,2);
AMPLITUDETHRESHOLD = 0.5;  % Reviewer 1's suggested threshold
MINNSACS = 10;
for a = 1:length(filenamelist)
    x = zeros(length(timebins),2);
    n = zeros(length(timebins),2);
    tmpdata = [];
    filenames = fnamesFromTxt(filenamelist{a});
    for b = 1:size(filenames,1)
        stro = nex2stro(findfile(filenames(b,:)));
        stro = DTfilterquesttrials(stro,'PaperDefaults');

        if isempty(stro.sum.analog.sigid) || isempty(stro.trial)
            continue;
        end
        
        ntrials = size(stro.trial,1);
        lms = stro.sum.exptParams.RF_colors;
        lms = reshape(lms,[3,size(lms,1)/3]);
        lms(:,all(lms == 0)) = [];
        lms = mkbasis(lms);
        whichcol = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'color_dir'));
        stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_on'));
        stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_off'));
        targon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targ_on'));
        correct = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'correct'));
        
        sacstats = getSacData(stro);
        close;
        
        % Looping over trials
        for j = 1:ntrials
        %    coloridx = find(abs(lms(:,whichcol(j))'*[.5774 .5774 .5774]') > .99);  % Avoiding roundoff errors
        %    if (coloridx ~= 1)
        %        continue;
        %    end
            dirs = sacstats.directions{j};
            amps = sacstats.amplitudes{j};
            st = sacstats.starttimes{j};
            L = logical(sacstats.starttimes{j} >= targon_t(j));
            st(L) = [];    % Omitting saccades that begin after targets appear
            amps(L) = [];  % Omitting saccades that begin after targets appear
            dirs(L) = [];  % Omitting saccades that begin after targets appear
            Lsac = (st > stimon_t(j)+offsets(1)) & (st < stimon_t(j)+offsets(2));
            if (any(Lsac))
                tmpdata = [tmpdata; correct(j) max(sacstats.amplitudes{j}(Lsac))];
            else
                tmpdata = [tmpdata; correct(j) nan];
            end
            for k = 1:length(timebins)
                Lsac = (st > stimon_t(j)+timebins(k)-timebinwidth/2) & (st < stimon_t(j)+timebins(k)+timebinwidth/2);
                if (any(Lsac))
                    ampidx = max(amps(Lsac) > AMPLITUDETHRESHOLD) + 1;
                    if (correct(j))
                        x(k,ampidx) = x(k,ampidx)+1;
                    end
                    n(k,ampidx) = n(k,ampidx)+1;
                end
            end
        end
    end
    
    Lnosac = isnan(tmpdata(:,2));
    nosacdata(1,a) = sum(tmpdata(Lnosac,1));
    nosacdata(2,a) = sum(Lnosac);
    pcornosac = sum(nosacdata(1,a))/sum(nosacdata(2,a));
    
    % P(correct) v time
    axes('Position',[1+3*(a-1) 3 AXESSIZE AXESSIZE]);
    hold on;
    colors = [1 0 0; 0 0 1];
    for ampidx = 1:2
        p = x(:,ampidx)./n(:,ampidx);
        se = sqrt((p.*(1-p))./n(:,ampidx));
        h = patch([timebins fliplr(timebins)],[p+se; flipud(p-se)]','b');
        set(h,'FaceColor',colors(ampidx,:));
        h = plot(timebins, x(:,ampidx)./n(:,ampidx),'k-','Linewidth',2);
        set(h,'Color',colors(ampidx,:));
    end
    ylabel('P(cor | saccade at t)','FontSize',12);
    xlabel('Time (s)','FontSize',12);
    set(gca,'Xlim',[timebins(1) timebins(end)+eps])
    if (a == 1)
        set(gca,'YLim',[.3 1]);
    else
       set(gca,'YLim',[.5 1]);
    end
    plot([0 .85],[pcornosac pcornosac],'k:','LineWidth',1);


    % P(correct) v amplitude
    axes('Position',[1+3*(a-1) 5.5 AXESSIZE AXESSIZE]);
    hold on;
    binedges = linspace(.9*nanmin(tmpdata(:,2)),1.1*nanmax(tmpdata(:,2)),NSACBINEDGES);
    for j = 2:length(binedges)
        L = tmpdata(:,2)>binedges(j-1) & tmpdata(:,2)<binedges(j);
        binmiddle = (binedges(j)+binedges(j-1))/2;
        sacdata(:,j-1,a) = [binmiddle sum(tmpdata(L,1)) sum(L)];
    end
    L = sacdata(3,:,a) < MINNSACS;
    processeddata = sacdata(:,~L,a);
    pcor = processeddata(2,:)./processeddata(3,:);
    patch([processeddata(1,:) fliplr(processeddata(1,:))],...
        [pcor + sqrt((pcor.*(1-pcor))./processeddata(3,:))...
        fliplr(pcor - sqrt((pcor.*(1-pcor))./processeddata(3,:)))],[.5 .5 .5]);
        plot(processeddata(1,:),pcor,'k.-','LineWidth',2)

    plot([0 .85],[pcornosac pcornosac],'k:','LineWidth',1);
    set(gca,'XLim',[0 .85]);
    ylabel('P(cor | amplitude)','FontSize',12);
    xlabel('Amplitude (deg)','FontSize',12);
    
    [b,bse] = lscov([processeddata(1,:)' ones(length(pcor),1)],pcor',diag(1./se(~L).^2));
    disp(['slope = ',num2str(b(1))]);
    p = tcdf(b(1)/bse(1),sum(~L)-2)
    
    % Monkey label
    axes('Position',[1+3*(a-1) 7.5 AXESSIZE 1],'Box','off','Visible','off');
    text(0.5,0,['Monkey ',filenames(1,1)],'HorizontalAlignment','center','FontSize',16);
    
    % Doing some stats --------------------------
    % T-test on amplitudes preceding correct vs. incorrect
    disp(['Doing stats on monkey ',num2str(a)]);
    Lnosac = isnan(tmpdata(:,2));
    mean(tmpdata(~Lnosac & tmpdata(:,1) == 0,2))
    mean(tmpdata(~Lnosac & tmpdata(:,1) == 1,2))
    [h,p] = ttest2(tmpdata(~Lnosac& tmpdata(:,1) == 0,2),tmpdata(~Lnosac&tmpdata(:,1) == 1,2))
    
    % GLM
    %[b,dev,stats] = glmfit(tmpdata(:,2),tmpdata(:,1),'binomial');
    %stats.p
    
    % Equal proportion test of big saccade vs small saccade trials
    L_bigsac = logical(tmpdata(:,2) >= 0.5);
    sum(L_bigsac & ~Lnosac)/length(L_bigsac& ~Lnosac)  % What fraction of saccade trials have big microsaccades?
    sum(tmpdata(L_bigsac& ~Lnosac,1))./sum(L_bigsac& ~Lnosac)
    sum(tmpdata(~L_bigsac& ~Lnosac,1))./sum(~L_bigsac& ~Lnosac)
    [h,p] = equalproptest([sum(tmpdata(L_bigsac& ~Lnosac,1)) sum(tmpdata(~L_bigsac& ~Lnosac,1))],[sum(L_bigsac& ~Lnosac) sum(~L_bigsac& ~Lnosac)],0.05)
    disp('End of stats');
    % ----------------------------

end
%%
% Section 11
% Charlie's code to look at differences b/w RG and LUM cells' post saccade
% spike densities. Playing around with different ways to normalize the
% PSDs.
BINWIDTH = .01;
timebins = [-.25:BINWIDTH:.5];
BASELINETIME = [-.25 0];
DIPTIME = [.05 .1];
LBASELINETIME = timebins >= BASELINETIME(1) & timebins <= BASELINETIME(2);
LDIPTIME = timebins >= DIPTIME(1) & timebins <= DIPTIME(2);

filenamelists = {'N:\NexFiles\nexfilelists\Greg\LvsM.txt',...
    'N:\NexFiles\nexfilelists\Greg\Lum.txt'};
data = [];
whichlist = [];
for a = 1:length(filenamelists)
    [filenames, spikenums] = fnamesFromTxt2(filenamelists{a});
    tmpdata = nan*ones(length(filenames), length(timebins));
    for i = 1:size(filenames,1)
        stro = {};
        for f = 1:size(filenames{i},1)
            tmpstro = nex2stro(findfile(char(filenames{i}(f))));
            if (isempty(stro))
                stro = tmpstro;
            else
                stro = strocat(stro, tmpstro);
            end
        end
        
        noisetypeidx = find(strcmp(stro.sum.trialFields(1,:),'noise_type'));
        stimonidx = find(strcmp(stro.sum.trialFields(1,:),'stim_on'));
        nframesidx = find(strcmp(stro.sum.trialFields(1,:),'num_frames'));
        
        if (isempty(stro.sum.analog.sigid))
            continue;
        end
        sacstats = getSacData(stro);
        close;
        L = logical(stro.trial(:,noisetypeidx) == 1);
        
        PSTH = zeros(1,length(timebins));
        for j = find(L')
            stimon_t = stro.trial(j,stimonidx);
            numframes = stro.trial(j,nframesidx);
            stimoff_t = stimon_t+numframes/stro.sum.exptParams.framerate;
            st = sacstats.starttimes{j};
            Lsac = (st > stimon_t) & (st < stimoff_t-timebins(end));
            if any(Lsac)
                Lsac(sacstats.amplitudes{j}(Lsac) > 2) = [];
                spiketimes = [];
                for k = find(Lsac')
                    tmp = stro.ras{j,spikenums(i)}-sacstats.starttimes{j}(k);
                    spiketimes = [spiketimes; tmp((tmp > timebins(1)-BINWIDTH/2) & (tmp < timebins(end)+BINWIDTH/2))];
                end
                [n,x] = hist(spiketimes, timebins);
                PSTH = PSTH + n./(BINWIDTH*sum(Lsac));
            end
        end
        tmpdata(i,:) = PSTH./sum(L);
    end
    whichlist = [whichlist; repmat(a,length(filenames),1)];
    data = [data; tmpdata]; 
end

%parse the data between the RG and LUM cells.
lvmpsth = data((whichlist == 1),:);
lumpsth = data((whichlist == 2), :);

%determine the percent change in firing rate due to a saccade.
preBins = (timebins >= -0.250) & (timebins < 0);
postBins = (timebins > 0) & (timebins <= 0.400);
dipBins = (timebins > 0.05) & (timebins <= 0.100);
peakBins = (timebins > 0.120) & (timebins <= .300);
lvmDeltaRate = (mean(lvmpsth(:, preBins),2) - mean(lvmpsth(:, dipBins),2))./mean(lvmpsth(:, preBins),2);
lumDeltaRate = (mean(lumpsth(:, preBins),2) - mean(lumpsth(:, dipBins),2))./mean(lumpsth(:, preBins),2);

lvmDipRate = mean(lvmpsth(:, dipBins),2)./mean(lvmpsth(:, preBins),2);
lumDipRate = mean(lumpsth(:, dipBins),2)./mean(lumpsth(:, preBins),2);

lvmPeakRate = mean(lvmpsth(:, peakBins),2)./mean(lvmpsth(:, preBins),2);
lumPeakRate = mean(lumpsth(:, peakBins),2)./mean(lumpsth(:, preBins),2);

figure
subplot(1,2,1)
hold on,
hist(lvmDeltaRate)
plot(mean(lvmDeltaRate), 4, 'mv', 'markerfacecolor', 'm')
title('lvmDeltaRate')
xlabel('%change')
subplot(1,2,2)
hold on
hist(lumDeltaRate)
plot(mean(lumDeltaRate), 2, 'mv', 'markerfacecolor', 'm')
title('lumDeltaRate')
xlabel('%change')
[hPrcntChange,pPrcntChange] = ttest2(lvmDeltaRate(:), lumDeltaRate(:))



%redo the analysis on normalized firing rates, but here I'm just calculating
%the pre/post diff in normalized rate...
lvmNorm = lvmpsth./repmat(max(lvmpsth,[],2),1,size(lvmpsth,2));
lumNorm = lumpsth./repmat(max(lumpsth,[],2),1,size(lumpsth,2));
preBins = (timebins >= -0.250) & (timebins < 0);
postBins = (timebins > 0) & (timebins <= 0.400);
lvmDeltaRate = mean(lvmNorm(:, preBins),2) - mean(lvmNorm(:, postBins),2);
lumDeltaRate = mean(lumNorm(:, preBins),2) - mean(lumNorm(:, postBins),2);


lvmPeakNormToMaxRate = mean(lvmNorm(:, peakBins),2);
lumPeakNormToMaxRate = mean(lumNorm(:, peakBins),2);

lvmDipNormToMaxRate = mean(lvmNorm(:, dipBins),2);
lumDipNormToMaxRate = mean(lumNorm(:, dipBins),2);


figure
subplot(1,2,1)
hold on,
hist(lvmDeltaRate)
plot(mean(lvmDeltaRate), 4, 'mv', 'markerfacecolor', 'm')
title('lvmDeltaRate')
xlabel('delta norm rate')
subplot(1,2,2)
hold on
hist(lumDeltaRate)
plot(mean(lumDeltaRate), 2, 'mv', 'markerfacecolor', 'm')
title('lumDeltaRate')
xlabel('delta norm rate')
[hNorm,pNorm] = ttest2(lvmDeltaRate(:), lumDeltaRate(:))


% recreate greg's fig
figure, hold on, %sould be exactly the same as greg's fig
plot(timebins, mean(lumpsth./repmat(max(lumpsth, [], 2), 1, size(lumpsth,2))), 'k')
plot(timebins, mean(lvmpsth./repmat(max(lvmpsth, [], 2), 1, size(lvmpsth,2))), 'r')
plot(timebins(preBins), repmat(mean(mean(lvmNorm(:, preBins),2)), 1, length(timebins(preBins))), 'r--')
plot(timebins(preBins), repmat(mean(mean(lumNorm(:, preBins),2)), 1, length(timebins(preBins))), 'k--')
plot(timebins(postBins), repmat(mean(mean(lvmNorm(:, postBins),2)), 1, length(timebins(postBins))), 'r--')
plot(timebins(postBins), repmat(mean(mean(lumNorm(:, postBins),2)), 1, length(timebins(postBins))), 'k--')
title('normalized firing rate')
xlabel('time')
ylabel('normalized rate')
hold off


% mean rate vs. time
figure, hold on, 
plot(timebins, mean(lumpsth,1), 'k')
plot(timebins, mean(lvmpsth,1)+1.5, 'r')
title('Firing Rate vs. Time')
xlabel('time')
ylabel('firing rate')


%percent change in firing rate 1-(rate./mean(baseline))
lvmBaseline = mean(lvmpsth(:, preBins),2);
lumBaseline = mean(lumpsth(:, preBins),2);
figure, hold on,
plot(timebins, mean((lvmpsth./repmat(lvmBaseline, 1, size(lvmpsth,2))),1), 'r')
plot(timebins, mean((lumpsth./repmat(lumBaseline, 1, size(lumpsth,2))),1), 'k')

%%
% Section 12
% An animated GIF showing fixational eye movements.

% Whitenoise file
filename = 'K031808005.nex';
stro = nex2stro(['/Users/greghorwitz/NexFiles/Greg/Kali/2008/',filename]);
stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_on'));
stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'all_off'));

% DTspot file
filename = 'S051810003.nex';
stro = nex2stro(['/Users/greghorwitz/NexFiles/Greg/Sedna/',filename]);
stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'frame_on'));
stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'choice_time'));


H = stro.ras(:,strcmp(stro.sum.rasterCells(1,:),'AD11'));
V = stro.ras(:,strcmp(stro.sum.rasterCells(1,:),'AD12'));
T = [stro.ras{:,strcmp(stro.sum.rasterCells(1,:),'anlgStartTime')}]'



whichtrial = 7;  % 7, 12 are good
h = H{whichtrial};
v = V{whichtrial};
t = T(whichtrial)+.002*[0:length(h)-1];
L = t > stimon_t(whichtrial) & t < stimoff_t(whichtrial)-1;
h = h(L)-mean(h(L)); v = v(L)-mean(v(L));
kernel = [.2 .5 1 .5 .2];
h = conv(h,kernel/sum(kernel));
v = conv(v,kernel/sum(kernel));
h = h(1:length(kernel):end);  % Subsampling
v = v(1:length(kernel):end);  % Subsampling

figure;
ax = axes; set(ax,'Units','pixels','position',[0 0 96 96]);
im = uint8(zeros(96,96,1,length(h)));
offset = 48;
gain = 1500;
[fpx,fpy] = meshgrid(offset+[-5:5]);
gauskern = normpdf([-4:.4:4],0,1)'*normpdf([-4:.4:4],0,1);
for i = 1:length(h)
    tmp1 = uint8(zeros(96,96));
    tmp2 = zeros(96,96);
    tmp1(fpx,fpy) = 1;
    tmp2(round(gain*v(i)+offset),round(gain*h(i)+offset)) = 1;
    tmp2 = round(conv2(tmp2,100*gauskern,'same'));
    tmp2(tmp2 <2) = 0;
    
    im(:,:,1,i) = tmp1+uint8(tmp2);
end
cmap = jet(30);
cmap(1,:) = [.5 .5 .5];
cmap(2,:) = [0 0 0];
close;

imwrite(im,cmap,'EP.gif','gif','DelayTime',0,'LoopCount',inf)

