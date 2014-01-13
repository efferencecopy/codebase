%% LINEAR MECHANISMS MODEL
%cycle through the valid neurons and compare the actual and predicted TRs
%of the intermediate color direction based on a linear model
l_validConds = ~(out.errors(:,cardTieInd) | out.errors(:,intTieInd) | out.errors(:, lt10spikesInd) | out.errors(:, lt24DTtrialsInd) | out.errors(:,lt5GTtrialsInd));
validCells = find(l_validConds);
for a = 1:length(validCells)
    cellInd = validCells(a);
    
    cardNT = out.dat(cellInd).c.alpha(out.maskCardTR(cellInd,:)');
    cardColor = out.dat(cellInd).expt.standColors(out.maskCardTR(cellInd,:)', :);
    cardUnit(a,:) = cardColor ./ norm(cardColor);
    cardVect = cardUnit(a,:)*cardNT;
    
    measIntNT(a) = out.dat(cellInd).c.alpha(out.maskIntTR(cellInd,:)');
    intColor = out.dat(cellInd).expt.standColors(out.maskIntTR(cellInd,:)', :);
    intUnit(a,:) = intColor ./ norm(intColor);
    intVect = intUnit(a,:)*measIntNT(a);
     
    % Under the linear model, a light is detected when the dot product of any light 
    % onto the mechanism, cardUnit(a,:), equals cardNT
    
    % It works for the cardinal light at threshold...
    [cardVect*cardUnit(a,:)'  cardNT]
    
    % ...but it doesn't work for the intermediate light at threshold.
    [abs(intVect*cardUnit(a,:)')  cardNT]
    
    % What do we have to multiply intVect by so that its
    % dot product onto the mechanism, cardUnit(a,:), is cardNT?
    % (x*light2*mech) = light1*mech
    % x = (light1*mech)/(light2*mech)
    intScaleFactor = (cardVect*cardUnit(a,:)') / abs(intVect*cardUnit(a,:)');
    
    % intScaleFactor is how much do we have to multiply intVect by so that
    % the dotproduct between intVect and cardUnit is same as dot prod
    % between cardVect and cardUnit.
    % The two expressions below should evaluate to the same number.
    % They are (1) the projection of the cardinal light onto the
    % cardinal mechanism and (2) the projection of the scaled intermediate
    % light onto the cardinal mechanism.
    [cardVect*cardUnit(a,:)' abs(intScaleFactor*intVect*cardUnit(a,:)')]
    
    % We can use this to derive a prediction about the threshold in the
    % intermediate color direction.
  	predIntNT(a) = intScaleFactor*norm(intVect)
    
    % To make sure this analysis is invariant to how we represent the data,
    % we are going to linearly transformation of the space by multiplying 
    % everything by some random matrix, A.
    A = normrnd(0,1,3,3);
    mech = inv(A')*(cardUnit(a,:))';
    light1 = A*cardVect';
    light2 = A*intVect';

    % Now just repeating the above code with the new mechanism and new
    % lights.
    newScaleFactor(a) = (mech'*light1) / abs(mech' * light2);
    newpredIntNT(a) = intScaleFactor*norm(light2);
    newmeasIntNT(a) = norm(light2);
end


% Hopefully, the ratio of predicted to actual neurometric thresholds
% is the same with and without the random transformation of the space.  
% Otherwise, the analysis is lame.
newpredIntNT./newmeasIntNT
predIntNT./measIntNT
% On the other hand, we expect the absolute values of the thresholds to
% depend on the linear transformation (like converting from
% inches to centimeters - the absolute values of the numbers change but the
% ratios don't).

%plot the data
figure
subplot(1,2,1), hold on
plot(measIntNT,predIntNT , 'b.')
xlabel('measured Int NT')
ylabel('predicted Int NT')
maxVal = max([measIntNT(:) ; predIntNT(:)]) .* 1.1;
minVal = min([measIntNT(measIntNT>0)' ; predIntNT(predIntNT>0)']) .* 0.90;
plot([minVal maxVal], [minVal, maxVal], 'k-')
nNans = sum(isnan(measIntNT(:)+predIntNT(:))); %either the pred or meas TR is a nan
title(sprintf('n = %d, nNan = %d', sum(l_validConds), nNans));
axis tight
hold off

subplot(1,2,2), hold on
plot(predIntNT./measIntNT, newpredIntNT./newmeasIntNT, 'b.')
xlabel('Pred/Meas NT in cone space');
ylabel('Pred/Meas NT in ???? space');
x = get(gca,'XLim');
plot(x.*[.9 1.1], x.*[.9 1.1],'k-');
axis tight


%%
% Fitting a half-squaring (with a threshold) to simulated data
b = [3 1 2];
contrasts = linspace(0,4,6);
lambda = b(1)+b(3)*(max(contrasts-b(2),0)).^2;

ntrials = 7;
nspikes = nan*ones(ntrials,length(contrasts));
for i = 1:ntrials
    nspikes(i,:) = poissrnd(lambda,1,length(lambda));
end

% OK, I've got my fake data. Now I just need to fit it.
x = repmat(contrasts,ntrials,1); % making a design matrix
%beta = glmfit(x(:),nspikes(:),'poisson','link','identity');
beta = regress(nspikes(:),[ones(length(x(:)),1) x(:).^2]);
%

options = optimset;
options = optimset(options,'Diagnostics','off','Display','off');
options = optimset(options,'LargeScale','off');
vlb = [0 0 0];
vub = [50 50 50];

params0 = [beta(1) .1 beta(2)];
c = repmat(contrasts,ntrials,1);
CAHfit(params0,c(:),nspikes(:))

f1 = fmincon('CAHfit',params0,[],[],[],[],vlb,vub,[],options,c(:),nspikes(:));

figure; axes; hold on;
plot(contrasts,nspikes,'k.');
x = linspace(contrasts(1),contrasts(end),30);
% plot(x,beta(1)+beta(2)*x.^2,'b-');  % Initial fit
plot(x,b(1)+b(3)*(max(x-b(2),0).^2),'k:'); % truth
plot(x,f1(1)+f1(3)*(max(x-f1(2),0).^2),'m-'); % fmincon fit

% Here's the -llik  of the fit
CAHfit(f1,c(:),nspikes(:))

% ---------------------------------------
% Now trying a *pair* of color directions
% ---------------------------------------
b_int = [b(1) b(2) b(3)]
contrast_scaling = 1;  % bounded between .1 and 10
contrasts_int = [0:1:4]/contrast_scaling;
lambda_int = b_int(1)+b_int(3)*(max(contrasts_int.*contrast_scaling-b_int(2),0).^2);

ntrials_int = 7;
nspikes_int = nan*ones(ntrials_int,length(contrasts_int));
for i = 1:ntrials_int
    nspikes_int(i,:) = poissrnd(lambda_int,1,length(lambda_int));
end

% --------------------
% Fitting them individually
x = repmat(contrasts,ntrials,1); % making a design matrix
beta1 = regress(nspikes(:),[ones(length(x(:)),1) x(:).^2]);
params0 = [beta1(1) 0 beta1(2)];
c = repmat(contrasts,ntrials,1);
f1 = fmincon('CAHfit',params0,[],[],[],[],vlb,vub,[],options,c(:),nspikes(:));
dev1 = CAHfit(f1,c(:),nspikes(:));

x = repmat(contrasts_int,ntrials_int,1); % making a design matrix
beta2 = regress(nspikes_int(:),[ones(length(x(:)),1) x(:).^2]);
params0 = [beta2(1) 0 beta2(2)];
c = repmat(contrasts_int,ntrials_int,1);
f2 = fmincon('CAHfit',params0,[],[],[],[],vlb,vub,[],options,c(:),nspikes_int(:));
dev1(2)  = CAHfit(f2,c(:),nspikes_int(:));
dev1 = sum(dev1);

% plotting
figure; axes; hold on;
plot(contrasts,nspikes','k.');
plot(contrasts_int,nspikes_int','m.')
x = linspace(0,max([contrasts, contrasts_int]),100);
plot(x,f1(1)+f1(3)*(max(x-f1(2),0).^2),'k-');
plot(x,f2(1)+f2(3)*(max(x-f2(2),0).^2),'m-');
set(gca,'Ylim',[0 max([nspikes_int(:);nspikes(:)])]);

% --------------------
% Trying to fit both yoked
vlb = [0 0 0 .1];
vub = [50 50 50 10];

params1 = [mean([f1(1); f2(1)]) f1(2) f1(3) min([sqrt(f2(3)./f1(3)), f1(2)./f2(2)])];
params1 = [mean([f1(1); f2(1)]) f1(2) f1(3) sqrt(f2(3)./f1(3))];
params1 = [mean([f1(1); f2(1)]) f1(2) f1(3) contrast_scaling]; % cheating
c1 = repmat(contrasts,ntrials,1);
c2 = repmat(contrasts_int,ntrials_int,1);
X = [c1(:);c2(:)];
X(:,2) = [zeros(1,numel(c1)) ones(1,numel(c2))];
y = [nspikes(:); nspikes_int(:)];

f = fmincon('CAHfit',params1,[],[],[],[],vlb,vub,[],options,X,y);
dev2 = CAHfit(f,X,y);

% plotting
figure; axes; hold on;
plot(contrasts,nspikes','k.');
plot(contrasts_int,nspikes_int','m.')
x = linspace(0,max([contrasts, contrasts_int]),30);
plot(x,f(1)+f(3)*(max(x-f(2),0).^2),'k-'); % fmincon fit
plot(x,f(1)+f(3)*(max(x*f(4)-f(2),0).^2),'m-'); % fmincon fit
set(gca,'Ylim',[0 max([nspikes(:);nspikes_int(:)])]);

devstat = 2*(dev2-dev1);  % "full model" goes second. dev1 and dev2 are -1*llik
p = 1- chi2cdf(devstat,2);
title(['p = ',num2str(p)]);

%[b(1) b(2) b(3) b_int(3)]
%params1


%%
% L, M, and S cone weights for an example neuron stimulated with Gaussian
% RGB noise. This is for a figure to give to the DTspot manuscript reviewers
whichframes = [5 6 7];

WN=nex2stro;
framerate = WN.sum.exptParams.framerate;
nstixperside = WN.sum.exptParams.nstixperside;
ntrials = length(WN.sum.absTrialNum);
stimonidx = find(strcmp(WN.sum.trialFields(1,:),'stim_on'));
stimoffidx = find(strcmp(WN.sum.trialFields(1,:),'all_off'));
nframesidx = find(strcmp(WN.sum.trialFields(1,:),'num_frames'));
noisetypeidx = find(strcmp(WN.sum.trialFields(1,:),'noise_type'));
sigmaidxs = strmatch('sigma',WN.sum.trialFields(1,:));

hepidx = find(strcmp(WN.sum.rasterCells(1,:),'AD11'));
vepidx = find(strcmp(WN.sum.rasterCells(1,:),'AD12'));
anlgStartTimeidx = find(strcmp(WN.sum.rasterCells(1,:),'anlgStartTime'));
eyestart_t = [WN.ras{:,anlgStartTimeidx}]';
eyesampperiod = 1/WN.sum.analog.storeRates{1};
gammaTable = WN.sum.exptParams.gamma_table;
gammaTable = reshape(gammaTable, length(gammaTable)/3, 3);
gammaTable1 = interp1(linspace(0,255,256),gammaTable,linspace(0,255,65536), 'spline');
invgamma = InvertGamma(gammaTable, 0);

% Reconstructing the M matrix
fundamentals = WN.sum.exptParams.fundamentals;
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = WN.sum.exptParams.mon_spd;
mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals'*mon_spd;

% Getting the background rgb/lms
ridx = find(strcmp(WN.sum.trialFields(1,:),'bkgnd_r'));
gidx = find(strcmp(WN.sum.trialFields(1,:),'bkgnd_g'));
bidx = find(strcmp(WN.sum.trialFields(1,:),'bkgnd_b'));
bkgndRGB = [mode(WN.trial(:,ridx)), mode(WN.trial(:,gidx)), mode(WN.trial(:,bidx))];
bkgndrgb = [gammaTable(bkgndRGB(1)+1,1); gammaTable(bkgndRGB(2)+1,2); gammaTable(bkgndRGB(3)+1,3)];
bkgndlms = M*bkgndrgb;
spikename = getSpikenum(WN);
spikeidx = find(strcmp(WN.sum.rasterCells(1,:),spikename));
maxT = 9;
Lgunnoise = WN.trial(:,noisetypeidx) == 1;
Lconenoise = WN.trial(:,noisetypeidx) == 2;
WN.ras(Lconenoise,:) = [];
WN.trial(Lconenoise,:) = [];
out = getWhtnsStats(WN,maxT,'STCOVmex', {nstixperside^2, 3, maxT}, spikename);


nspikes = out{3};
STAs = [];
for i = 1:length(whichframes)
    lms = inv(M')*reshape(out{1}(:,whichframes(i)),100,3)';
    lms = lms';
  %  lms = reshape(out{1}(:,whichframes(i)),100,3); % debugging
    STAs(:,i) = lms(:);
end
STAs = mean(STAs,2);
maxes = max(STAs(:));
mins = min(STAs(:));

potentialnormfactors = [(1-[.5; .5; .5]-eps)./(maxes-[.5; .5; .5]); (-[.5; .5; .5]+eps)./(mins-[.5; .5; .5])];
% 'eps' in above line is a kludge that is required for avoiding
% out of bounds errors.
potentialnormfactors(potentialnormfactors < 0) = []; % if min > mu or max < mu
normfactor = min(potentialnormfactors);

%muvect = reshape(repmat(bkgndrgb',nstixperside^2,1),nstixperside^2*3,1);
muvect = reshape(repmat([.5 .5 .5],nstixperside^2,1),nstixperside^2*3,1);

% Plotting
figure;
for i = 1:size(STAs,2)
    STA = normfactor*(STAs(:,i)-muvect)+muvect;
    STA = reshape(STA,[nstixperside nstixperside 3]);
    for j = 1:3
        subplot(3,size(STAs,2),j);
        image(255*STA(:,:,j)); colormap(gray(255));
        set(gca,'XTick',[],'YTick',[]); axis square;
    end
end