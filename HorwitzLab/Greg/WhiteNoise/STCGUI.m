% STCGUI.m
%
% Code for looking at various projections of the spike-triggered
% stimuli.  This will create a little interactive GUI for exploring
% the PCs of spike-triggered stimuli.
% Use the left mouse button to select and deselect pixels.
% Use the right mouse button to plot STA, PC1, and PC2 in 
% (x,y) (x,t) (y,t) and (t) slices.
%
% GDLH 3/5/09

function STCGUI(STA, STC, nstixperside, maxT, filename)
    
    global axesh;
    global data;
    global uicontrols;
    
    if (nargin == 0)
       [STA, STC, nstixperside, maxT, filename] = loadData; 
    end
    
    % Initializing some variables
    STAs = STA./(2*max(abs(STA(:))))+.5;
    STAs = reshape(STAs,[nstixperside nstixperside 3 maxT]);
    pixelmask = ones(nstixperside,nstixperside);
    SetUpFigs;
  
    % If data aren't passed in as arguments to STCgui, load them in here.
    function [STA, STC, nstixperside, maxT, filename] = loadData()

        [fnames, pathname] = uigetfile('*.nex', 'Select NeuroExplorer file(s)','MultiSelect','on');
        if (~iscell(fnames))
            stro = nex2stro([pathname,'/',fnames]);
        else
            stro = {};
            for i = 1:length(fnames)
                tmpstro = nex2stro([pathname,'/',char(fnames{i})]);
                if (isempty(stro))
                    stro = tmpstro;
                else
                    stro = strocat(stro, tmpstro);
                end
            end
        end

        maxT = 9;
        nstixperside = stro.sum.exptParams.nstixperside;
        noisetypeidx = find(strcmp(stro.sum.trialFields(1,:),'noise_type'));
        
        % Getting spike number and noise type
        spikename = getSpikenum(stro);
        Lgunnoise = stro.trial(:,noisetypeidx) == 1;
        Lconenoise = stro.trial(:,noisetypeidx) == 2;
        choice = 0;
        if (sum(Lgunnoise) > 0 && sum(Lconenoise) > 0)
            liststring = {['Gun (n = ',num2str(sum(Lgunnoise)),')'] ,['Cone (n = ',num2str(sum(Lconenoise)),')']};
            choice = listdlg('PromptString','Which noise type?','SelectionMode','Single','ListSize',[150 50],...
                'ListString',liststring);
        end
        if (choice == 1 || sum(Lconenoise) == 0)
            stro.ras(Lconenoise,:) = [];
            stro.trial(Lconenoise,:) = [];
        elseif (choice == 2 || sum(Lgunnoise) == 0)
            stro.ras(Lgunnoise,:) = [];
            stro.trial(Lgunnoise,:) = [];
        end
        filename = stro.sum.fileName;
        out = getWhtnsStats(stro,maxT,'STCOVfull', {nstixperside^2, 3, 1, maxT}, spikename);
        STA = out{1};
        STC = out{2};
    end

    function SetUpFigs()
        figure(1); clf;
        set(gcf,'NumberTitle','off','Name', filename);
        set(gcf,'Position', [352 793 599 165]);
        for i = 1:maxT
            subplot(1,maxT,i);
            h = image(squeeze(STAs(:,:,:,i)));
            set(h,'ButtonDownFcn',{@UpdatePlots,i});
            set(gca,'XTick',[],'YTick',[]);
            axis square;
        end
        figure(2); clf;
        set(gcf,'Position', [352 229 599 513],'DefaultAxesColorOrder',[1 0 0; 0 1 0; 0 0 1]);
        uicontrols.projorthosta = uicontrol('style','radio','Position',[515,445,10,10],'Value',0,'CallBack','');
        uicontrols.projtext = uicontrol('style','text','Position',[450,460,150,15],'String','Project orthogonal to STA');
        uicontrols.loworder = uicontrol('style','radio','Position',[515,400,10,10],'Value',0,'CallBack','');
        uicontrols.lowordertext = uicontrol('style','text','Position',[475,415,100,15],'String','Low order PCs');

        expansionfactor = 8;
        [axesh.STA.xy, axesh.STA.xt, axesh.STA.ty, axesh.STA.t] = SetUp4Axes([50 370], 10);
        [axesh.D.xy, axesh.D.xt, axesh.D.ty, axesh.D.t] = SetUp4Axes([50 150], 10);
        [axesh.V1.xy, axesh.V1.xt, axesh.V1.ty, axesh.V1.t] = SetUp4Axes([300 370], 10);
        [axesh.V2.xy, axesh.V2.xt, axesh.V2.ty, axesh.V2.t] = SetUp4Axes([300 150], 10);
        
        function [xy, xt, ty, t] = SetUp4Axes(xypos, margin)
            xy = axes('Units','pixels');
            set(xy,'Position',[xypos(1), xypos(2), expansionfactor*nstixperside, expansionfactor*nstixperside]);
            xt = axes('Units','pixels');
            set(xt,'Position',[xypos(1), xypos(2)-expansionfactor*maxT-margin, expansionfactor*nstixperside, expansionfactor*maxT]);
            ty = axes('Units','pixels');
            set(ty,'Position',[xypos(1)+expansionfactor*nstixperside+margin, xypos(2), expansionfactor*maxT, expansionfactor*nstixperside]);
            t = axes('Units','pixels');
            set(t,'Position',[xypos(1)+expansionfactor*nstixperside+margin, xypos(2)-expansionfactor*maxT-margin, expansionfactor*maxT, expansionfactor*maxT]);
        end
    end

    function UpdatePlots(h, ev, whichframe)
        whichpt = get(gca,'CurrentPoint');
        whichpt = round(whichpt(1,[1 2]));
        whichpt = min([whichpt; nstixperside nstixperside]);
        whichbutton = find(strcmp(get(gcf,'SelectionType'),{'normal','alt'}));
        if (whichbutton == 1)
            StartPixelSelect(whichpt(1), whichpt(2), whichframe);
        elseif (whichbutton == 2)
            ComputeEigenvectors(whichpt(1), whichpt(2), whichframe);
            UpdateSTAPlots(whichpt(1), whichpt(2), whichframe);
            UpdateEigenvaluePlots;
            UpdateEigenvectorPlots;
        end
    end

    function StartPixelSelect(x, y, whichframe)
        lastselectedpoint = [nan nan];
        drawmode = ~pixelmask(y,x); % 1 = add, 0 = erase 
        PixelSelect(nan, nan);
        set(1,'WindowButtonMotionFcn',@PixelSelect);
        set(1,'WindowButtonUpFcn',@StopPixelSelect);
        
        function PixelSelect(h,ev)
            whichpt = get(gca,'CurrentPoint');
            whichpt = round(whichpt(1,[1 2]));
            whichpt = min([whichpt; nstixperside nstixperside]);
            if any(whichpt < 1)
                return
            end
            if (~all(whichpt == lastselectedpoint))
                pixelmask(whichpt(2),whichpt(1)) = drawmode;
                lastselectedpoint = whichpt;
                UpdateMaskedSTA(whichframe);
            end
        end
        function StopPixelSelect(h,ev)
            set(1,'WindowButtonMotionFcn','')
            set(1,'WindowButtonUpFcn','')
        end
        
        function UpdateMaskedSTA(whichframe)
            h = image(squeeze(STAs(:,:,:,whichframe)).*repmat(pixelmask,[1 1 3]));
            set(h,'ButtonDownFcn',{@UpdatePlots,whichframe});
            set(gca,'XTick',[],'YTick',[]);
            axis square;
        end
    end
    
    function UpdateSTAPlots(x, y, whichframe)
        image(squeeze(STAs(:,:,:,whichframe)),'Parent',axesh.STA.xy);
        set(axesh.STA.xy, 'NextPlot','add');
        plot(x,y,'y.','Parent',axesh.STA.xy);
        set(axesh.STA.xy,'NextPlot','replace','XTick',[],'YTick',[],'XLim',[.5 nstixperside+.5],'YLim',[.5 nstixperside+.5]);
        
        image(permute(squeeze(STAs(y,:,:,:)),[3 1 2]),'Parent',axesh.STA.xt);
        set(axesh.STA.xt, 'NextPlot','add');
        plot(x,whichframe,'y.','Parent',axesh.STA.xt);
        set(axesh.STA.xt,'NextPlot','replace','XTick',[],'YTick',[],'XLim',[.5 nstixperside+.5],'YLim',[.5 maxT+.5]);
        
        image(permute(squeeze(STAs(:,x,:,:)),[1 3 2]),'Parent',axesh.STA.ty);
        set(axesh.STA.ty, 'NextPlot','add');
        plot(whichframe,y,'y.','Parent',axesh.STA.ty);
        set(axesh.STA.ty,'NextPlot','replace','XTick',[],'YTick',[],'XLim',[.5 maxT+.5],'YLim',[.5 nstixperside+.5]);
        
        plot(squeeze(STAs(y,x,:,:))','LineWidth',2,'Parent',axesh.STA.t);
        set(axesh.STA.t,'TickLength',[0 0],'YTick',[],'XTickLabel',[],'XLim',[.5 maxT+.5]);
    end

    function ComputeEigenvectors(x,y,whichframe)
        template = reshape(1:size(STC,1),nstixperside, nstixperside, 3, maxT);
        template = template .* repmat(pixelmask, [1 1 3 maxT]);
        
        [data.vxy, data.dxy] = EigStripCov(template(:,:,:,whichframe));
        [data.vxt, data.dxt] = EigStripCov(template(y,:,:,:));
        [data.vty, data.dty] = EigStripCov(template(:,x,:,:));
        [data.vt, data.dt] = EigStripCov(template(y,x,:,:));
        
        m = max([abs(data.vxy(:)); abs(data.vxt(:)); abs(data.vty(:))]);
        normfactor = (0.5-eps)./m;
        data.vxy = normfactor.*data.vxy;
        data.vxt = normfactor.*data.vxt;
        data.vty = normfactor.*data.vty;
        data.vt = normfactor.*data.vt;
        
        function [v,d] = EigStripCov(template)
            PROJORTHOSTA = get(uicontrols.projorthosta,'Value');
            LOWORDER = get(uicontrols.loworder,'Value');
            keepidxs = template(:);
            chuckidxs = keepidxs == 0;
            keepidxs(chuckidxs) = [];
            if (PROJORTHOSTA)
                subSTA = STA(:);
                subSTA = subSTA(keepidxs);
                P = eye(length(keepidxs),length(keepidxs))-subSTA/(subSTA'*subSTA)*subSTA';
                [tmpv,tmpd] = eig(P*STC(keepidxs,keepidxs)*P');
            else
                [tmpv,tmpd] = eig(STC(keepidxs,keepidxs));
            end
            v = zeros(size(chuckidxs,1),size(tmpv,1));
            v(~chuckidxs,:) = tmpv;
            if (LOWORDER)
                [d, idxs] = sort(diag(tmpd),1,'descend');  
            else
                [d, idxs] = sort(diag(tmpd),1,'ascend');            
            end
            v = v(:,idxs);
        end
    end

    function UpdateEigenvaluePlots
        plot(flipud(data.dxy),'k.','Parent',axesh.D.xy)
        set(axesh.D.xy,'XLim',[1 10],'XTick',[],'YTick',[]);
        
        plot(flipud(data.dxt),'k.','Parent',axesh.D.xt)
        set(axesh.D.xt,'XLim',[1 10],'XTick',[],'YTick',[]);
        
        plot(flipud(data.dty),'k.','Parent',axesh.D.ty)
        set(axesh.D.ty,'XLim',[1 10],'XTick',[],'YTick',[]);
        
        plot(flipud(data.dt),'k.','Parent',axesh.D.t)
        set(axesh.D.t,'XLim',[1 10],'XTick',[],'YTick',[]);
    end

    function UpdateEigenvectorPlots()
        whichvect = 1;
        
        for whichvect = [1,2]
            if whichvect == 1
                h = axesh.V1;
            else
                h = axesh.V2;
            end
            UpdateOneSetEigenvectorPlots;
        end
        
        function UpdateOneSetEigenvectorPlots
            try
                if (~isempty(data.vxy))
                    image(reshape(data.vxy(:,end-whichvect+1)+.5, [nstixperside, nstixperside 3]),'Parent',h.xy);
                    set(h.xy,'XTick',[],'YTick',[],'XLim',[.5 nstixperside+.5],'YLim',[.5 nstixperside+.5]);
                end
                if (~isempty(data.vxt))
                    image(permute(reshape(data.vxt(:,end-whichvect+1)+.5, [nstixperside, 3, maxT]),[3 1 2]),'Parent',h.xt);
                    set(h.xt,'XTick',[],'YTick',[],'XLim',[.5 nstixperside+.5],'YLim',[.5 maxT+.5]);
                end
                if (~isempty(data.vty))
                    image(permute(reshape(data.vty(:,end-whichvect+1)+.5, [nstixperside, 3, maxT]),[1 3 2]),'Parent',h.ty);
                    set(h.ty,'XTick',[],'YTick',[],'XLim',[.5 maxT+.5],'YLim',[.5 nstixperside+.5]);
                end
                if (~isempty(data.vt))
                    tmp = reshape(data.vt(:,end-whichvect+1)+.5,[3 maxT])';
                else
                    tmp = .5*ones(maxT, 3);
                end
                plot(tmp,'LineWidth',2,'Parent',h.t);
                set(h.t,'TickLength',[0 0],'YTick',[],'XTickLabel',[],'XLim',[.5 maxT+.5], 'YLim',[.9*min(tmp(:)) 1.1*max(tmp(:))]);
            catch
                keyboard
            end
        end
    end
end
