% CompareSTAandGrating
%
% A function for comparing tuning properties assessed with the grating
% paradigm and with the white noise paradigm
% The filelist that goes with this function is "Gratings&WN.txt".
% The reason that this code is in the form of a function instead of a
% script is so I can declare subfunctions that can be used as callbacks
% when you click the data points.
%
% GDLH  8/20/09
function CompareSTAandGrating()
    if (exist('tmpGH.mat','file'))
        data = struct2cell(load('tmpGH'));
        data = data{1};
    else
        data = DoIt;
    end

    % Setting up filters for various parameters up front
    % rather than independently for each analysis.
    Lcolor = ones(length(data),1);
    LSTAenergy = ones(length(data),1);
    Loriented = ones(length(data),1);
    Lrespmag = ones(length(data),1);
    for i = 1:length(data)
        if (isnan(data{i}.gratingstruct.color.prefcolor))
            Lcolor(i) = 0;
        end
        energy = data{i}.gaborstruct.energy;
        if (max(energy) < 1.5*energy(1))
            LSTAenergy(i) = 0;
        end
        if (data{i}.gaborstruct.prefSF <= 0.6 | data{i}.gratingstruct.sf.prefSF <= 0.3)
           Loriented(i) = 0; 
        end
        if (max(data{i}.gratingstruct.color.colresp) < 5);
           Lrespmag(i) = 0; 
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Making the figures
    % Orientation tuning
    figure; axes; hold on;
    title('ORIENTATION');
    plot([0 pi],[0 pi],'k:');
    STA_lambda = cellfun(@(x) x.gaborstruct.lambda,data);
    STA_theta = mod(cellfun(@(x) x.gaborstruct.theta,data),pi);
    Grating_theta = mod(cellfun(@(x) x.gratingstruct.orient.preforient,data),pi);
    for i =1:size(STA_theta,2)
        if (~LSTAenergy(i) | ~Loriented(i) | ~Lrespmag(i))
            continue
        end
        h = plot(STA_theta(i),Grating_theta(i),'k.','MarkerSize',20);
        set(h,'ButtonDownFcn',{@showSTAGrat,data{i},'orient'});
    end
    xlabel('STA orientation'); ylabel('Grating pref orient'); axis square;
    
    % Cone weights
    figure; axes; hold on;
    title('CONE WEIGHTS');
    set(gca,'XLim',[-1 1],'Ylim',[-1 1]);
    plot([1 0 -1 0 1],[0 1 0 -1 0],'k:');
    axis square;

    for i = 1:length(data)
        if (~Lcolor(i) | ~LSTAenergy(i) | ~Lrespmag(i))
            continue
        end
        data{i}.gratingstruct.color.prefcolor = data{i}.gratingstruct.color.prefcolor./sum(abs(data{i}.gratingstruct.color.prefcolor));
        plot([data{i}.gaborstruct.lms(1),data{i}.gratingstruct.color.prefcolor(1)],...
            [data{i}.gaborstruct.lms(2),data{i}.gratingstruct.color.prefcolor(2)],'k-');
        h(1) = plot(data{i}.gaborstruct.lms(1),data{i}.gaborstruct.lms(2),'ko','MarkerSize',6);  % white noise
        h(2) = plot(data{i}.gratingstruct.color.prefcolor(1),data{i}.gratingstruct.color.prefcolor(2),'mo','MarkerSize',6); % gratings
        for j = 1:2
            if (j == 1)
                s = data{i}.gaborstruct.lms(3);
            else
                s = data{i}.gratingstruct.color.prefcolor(3);
            end
            edgecol = get(h(j),'Color');
            if (s < 0)
                set(h(j),'MarkerFacecolor',edgecol);
            else
                set(h(j),'MarkerEdgeColor',edgecol/2+[.5 .5 .5],'MarkerFaceColor',edgecol/2+[.5 .5 .5]);
            end
            set(h(j),'ButtonDownFcn',{@showSTAGrat,data{i},'color'});
        end
    end
    legend([h(1);h(2)],{'STA','Gratings'});
    
    % Spatial frequency tuning
    figure; axes; hold on;
    title('SPATIAL FREQUENCY');
    plot([0 4],[0 4],'k:');
    for i = 1:length(data)
        if (~LSTAenergy(i) | ~Lrespmag(i))
            continue
        end
        h = plot(data{i}.gaborstruct.prefSF,data{i}.gratingstruct.sf.prefSF,'k.','MarkerSize',20);
        set(h,'ButtonDownFcn',{@showSTAGrat,data{i},'sf'});
    end
    xlabel('STA spatial freq (dva)'); ylabel('Grating spatial freq (dva)');
    axis square;
    set(gca,'XLim',[0 4],'YLim',[0 4]);

    
    %%%%%%%%%%%%%%%%%%%%%%%
    % Nested subfunctions
    function showSTAGrat(e,h,data,whichanal)
        h = figure;
        set(gcf,'NumberTitle','off')
        fname = data.gaborstruct.fname;
        idx = find(fname == abs('\'),1,'last');
        fname(1:idx) = [];
        set(gcf,'name',fname)
        set(h,'Units','normalized','position',[.2 .4 .2 .2],'Toolbar','none');

        im = data.gaborstruct.im./(2*max(abs(data.gaborstruct.im(:))));
        im = reshape(im,[sqrt(size(im,1)), sqrt(size(im,1)), 3]);
        im = im +.5*ones(size(im));
        subplot(1,2,1);
        imagesc(im);
        xlabel(['MR = ',num2str(data.gratingstruct.modulationratio)]);
        axis square;
        if (data.gaborstruct.noisetype == 2)
            title('Cone Noise Image');
        end
        axis square;
        set(gca,'XTick',[],'YTick',[])
        if (strcmp(whichanal,'orient'))
            subplot(1,2,2);
            axis square;
            set(gca,'XTick',[],'YTick',[]);
            resp = data.gratingstruct.orient.resp;
            stim = data.gratingstruct.orient.stim;
            sems = resp(:,2)./sqrt(resp(:,3)); 
            polar([stim; stim(1)],[resp(:,1); resp(1,1)]+[sems; sems(1)],'b:');
            hold on;
            h = polar([stim; stim(1)],[resp(:,1); resp(1,1)],'b-');
            set(h,'Linewidth',2);
            polar([stim; stim(1)],[resp(:,1); resp(1,1)]-[sems; sems(1)],'b:');
            title('orient');
        elseif (strcmp(whichanal,'sf'))
            subplot(1,2,2);
            set(gca,'XTick',[],'YTick',[]);
            resp = data.gratingstruct.sf.resp;
            stim = data.gratingstruct.sf.stim;
            sems = resp(:,2)./sqrt(resp(:,3)); 
            h = errorbar(stim, resp(:,1), sems, 'b-');
            set(h,'Linewidth',2);
            axis square;
            set(gca,'XScale','log');
            title('sf');
        elseif (strcmp(whichanal,'color'))
            subplot(1,2,2);
            axis square;
            set(gca,'XTick',[],'YTick',[]);
            disp('Color key')
            disp([[1:9]', data.gratingstruct.color.colors])
            bar(data.gratingstruct.color.colresp(:,1));
            hold on;
            plot([0 10],repmat(data.gratingstruct.baselines(1),1,2),'k:');
        end
    end

    % This is the subfunction that actually does the bulk of the work
    % (creating the 'data' structure.
    function data = DoIt()
        
        maxT = 8;
        data = {};
        [fnames, spikeIdx] = fnamesFromTxt2();
        MINRESP = 10;  % At least 10 sp/sec to preferred grating.
        
        for a = 1:size(fnames,1)
            WN = {}; GT = {};
            for i = 1:size(fnames{a},2)
                tmpstro = nex2stro(findfile(char(fnames{a}(i))));
                if tmpstro.sum.paradigmID == 150
                    GT = tmpstro;
                elseif (isempty(WN))
                    WN = tmpstro;
                else
                    WN = strocat(WN, tmpstro);
                end
            end
            
            % Useful for debugging
            disp(['saving... ',num2str(a)]);
            disp(GT.sum.fileName);
            save tmpGH data
            
            nstixperside = WN.sum.exptParams.nstixperside;
            noisetypeidx = find(strcmp(WN.sum.trialFields(1,:),'noise_type'));
            mondist = 100; % cm
            screenwidth = 36; %cm
            screenwidthinpix = 1024; % singlewide pixels
            pixpercm = screenwidthinpix/screenwidth;
            cmperdeg = screenwidth/(2*atan2(screenwidth/2, mondist)*180/pi);
            pixperdeg = pixpercm*cmperdeg;
            stixperdeg = pixperdeg/(2*WN.sum.exptParams.npixperstix);
            % Reconstructing the M matrix and gamma table
            fundamentals = WN.sum.exptParams.fundamentals;
            fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
            mon_spd = WN.sum.exptParams.mon_spd;
            mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
            mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
            M = fundamentals'*mon_spd;
            
            Lgun = WN.trial(:,noisetypeidx) == 1;
            Lcone = WN.trial(:,noisetypeidx) == 2;
            
            if (sum(Lgun) > sum(Lcone))
                noisetype = 1;
                WN.ras(~Lgun,:) = [];
                WN.trial(~Lgun,:) = [];
                
                out = getWhtnsStats(WN,maxT,'STCOVmex', {nstixperside^2, 3, maxT}, WN.sum.rasterCells{spikeIdx(a)});
                STAs = out{1};
                
                energy = sum(STAs.^2);
                whichframe = energy == max(energy);
                
                % fitting the gabor
                STAim = reshape(STAs(:,whichframe),nstixperside*nstixperside,3);
                [u,s,v] = svd(STAim');
                im = reshape(v(:,1),[nstixperside nstixperside]);
                gaborstruct = gaborfit(im);
                gaborstruct.sf = stixperdeg./gaborstruct.lambda;
                gaborstruct.rgb = u(:,1);
                gaborstruct.lms = inv(M')*u(:,1);
                gaborstruct.lms = gaborstruct.lms./sum(abs(gaborstruct.lms));
            else
                noisetype = 2;
                WN.ras(~Lcone,:) = [];
                WN.trial(~Lcone,:) = [];
                
                conesigmas = [WN.trial(1,strcmp(WN.sum.trialFields(1,:),'sigma1'));...
                    WN.trial(1,strcmp(WN.sum.trialFields(1,:),'sigma2'));...
                    WN.trial(1,strcmp(WN.sum.trialFields(1,:),'sigma3'))];
                out = getWhtnsStats(WN,maxT,'STCOVmex', {nstixperside^2, 3, maxT}, WN.sum.rasterCells{spikeIdx(a)});
                STAs = out{1};
                
                energy = sum(STAs.^2);
                whichframe = find(energy == max(energy));
                
                % fitting the gabor
                STAim = reshape(STAs(:,whichframe),nstixperside*nstixperside,3);
                [u,s,v] = svd(STAim');
                if (sum(v(:,1)) < 0)
                    u = -u;
                end
                im = reshape(v(:,1),[nstixperside nstixperside]);
                gaborstruct = gaborfit(im);
                gaborstruct.sf = stixperdeg./gaborstruct.lambda;
                u(:,1) = u(:,1)./conesigmas;
                gaborstruct.lms = u(:,1)./sum(abs(u(:,1)));
                gaborstruct.rgb = inv(M')*u(:,1);
                gaborstruct.rgb = gaborstruct.lms./sum(abs(gaborstruct.lms));
            end
            gaborstruct.energy = energy;
            gaborstruct.im = STAim;
            gaborstruct.noisetype = noisetype;
            mondist = 100; % cm
            screenwidth = 36; %cm
            screenwidthinpix = 1024; % singlewide pixels
            pixpercm = screenwidthinpix/screenwidth;
            cmperdeg = screenwidth/(2*atan2(screenwidth/2, mondist)*180/pi);
            pixperdeg = pixpercm*cmperdeg;
            gaborstruct.stixperdeg = pixperdeg/(2*WN.sum.exptParams.npixperstix);
            gaborstruct.prefSF = gaborstruct.stixperdeg./gaborstruct.lambda;
            gaborstruct.fname = WN.sum.fileName;
            
            gratingstruct = getGratingTuning(GT, spikeIdx(a));
            if (gratingstruct.color.prefcolor * gaborstruct.lms < 0)
                gratingstruct.color.prefcolor = -gratingstruct.color.prefcolor;
            end
            data{a}.gaborstruct = gaborstruct;
            data{a}.gratingstruct = gratingstruct;     
        end
    end
end