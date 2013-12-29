%%
% Section 1
% Plotting the spectrum of a just barely detectable L-M light (and the
% spectrum of the background)

stro = nex2stro(findfile('K040609001'));
bkgndrgb = [stro.sum.exptParams.bkgnd_r stro.sum.exptParams.bkgnd_g stro.sum.exptParams.bkgnd_b];
guns = reshape(stro.sum.exptParams.mon_spect,81,3);
x = linspace(380,780,81);
bkgndspect = guns*bkgndrgb'
M = reshape(stro.sum.exptParams.m_mtx,3,3);
bkgndlms = M*bkgndrgb';
threshlms = bkgndlms.*(1+[.01 -.01 0]');
(threshlms-bkgndlms)./bkgndlms  % just a sanity check
threshrgb = inv(M)*threshlms;
threshspect = guns*threshrgb;

figure; axes; hold on;
plot(x,bkgndspect);
plot(x,threshspect,'k-');
xlabel('nm');
ylabel('w/sr/m^2');

[thresholds, colorDirs, sfs, QuestTrajectories] = DTquestUnpackGH(stro, 'mode');
thresholdlms = bkgndlms


%%
% Section 2
% Continuing from previous section.  Stealing code from Psychtoolbox's
% "IsomerizationInEyeDemo.m"
whatCalc = 'LivingHumanFovea';
photoreceptors = DefaultPhotoreceptors(whatCalc);
photoreceptors = FillInPhotoreceptors(photoreceptors);
% GH converting watts/sr-m^2-wlinterval to  watts/sr-m^2
radianceWatts = SplineSpd([380 5 81],bkgndspect,[380 1 401]);
load T_xyz1931	
T_xyz = SplineCmf(S_xyz1931,683*T_xyz1931,[380 1 401]);
theXYZ = T_xyz*radianceWatts; theLuminance = theXYZ(2);
[nil,pupilAreaMM] = PupilDiameterFromLum(theLuminance,photoreceptors.pupilDiameter.source);
% GH This seems awfully big: 7.5 mm

radianceWatts = SplineSpd([380 5 81],bkgndspect,[380 1 401]);
irradianceWatts = RadianceToRetIrradiance(radianceWatts,[380 1 401], ...
			pupilAreaMM,photoreceptors.eyeLengthMM.value);
irradianceQuanta = EnergyToQuanta([380 1 401],irradianceWatts);
figure; set(gcf,'Position',[100 400 700 300]);
subplot(1,2,1); hold on
set(plot(SToWls([380 1 401]),irradianceQuanta,'r'),'LineWidth',1);
set(title('Light Spectrum'),'FontSize',14);
set(xlabel('Wavelength (nm)'),'FontSize',12);
set(ylabel('Quanta/sec-um^2-wlinterval'),'FontSize',12);
[isoPerConeSecBkgnd,absPerConeSec,photoreceptors] = ...
	RetIrradianceToIsoRecSec(irradianceWatts,[380 1 401],photoreceptors);

radianceWatts = SplineSpd([380 5 81],threshspect,[380 1 401]);
irradianceWatts = RadianceToRetIrradiance(radianceWatts,[380 1 401], ...
			pupilAreaMM,photoreceptors.eyeLengthMM.value);
irradianceQuanta = EnergyToQuanta([380 1 401],irradianceWatts);
subplot(1,2,1); hold on
set(plot(SToWls([380 1 401]),irradianceQuanta,'k'),'LineWidth',1);
[isoPerConeSecThresh,absPerConeSec,photoreceptors] = ...
	RetIrradianceToIsoRecSec(irradianceWatts,[380 1 401],photoreceptors);
%GH above line: what's the difference between the first two output
%arguments?  It depends on the difference between
%photoreceptors.effectiveAbsorbtance and photoreceptors.isomerizationAbsorbtance
subplot(1,2,2);
bar([isoPerConeSecBkgnd, isoPerConeSecThresh])
ylabel('Photoisomerizations/cone/sec');

%%
% Section 3
% For Juan Anguerya's general exam.  Looking at how receptoral and
% post-receptoral noise affect detection contours.
darknoise = 1;

sigmu = 1;
x = [-10:.1:10];
noi = normpdf(x,0,1);
sig = normpdf(x,sigmu,1);
plot(x,noi,'b-')
hold on
plot(x,sig,'g-');
% percent correct in seen/not seen with d' = 1 and optimal criterion
normcdf(sigmu/2,0,1)
% Looks like 69% correct with d'=1 and 84% with d'=2
[X,Y] = meshgrid(x,x);
% l_out = normcdf(X,1,receptornoise);
% m_out = normcdf(Y,1,receptornoise);
% lm_both = l_out .* m_out;
% % P(A or B) = P(A)+P(B)-P(A and B)
% imagesc(l_out+m_out-(l_out .* m_out))
% axis xy

%contour(X,Y,data,3)

%%
% OK, just brute forcing a simulation

x = [-10:.1:10];
l_noise = 1;
m_noise = 1;
lumnoise = 3;
rgnoise = 0;
lumgain = 1;
rggain = 4;
[X,Y] = meshgrid(x,x);
niter = 1000;
data = zeros(size(X));
noise = zeros(1,niter);
s_lum = sqrt(lumgain^2*l_noise^2+m_noise^2+lumnoise^2);
s_rg = sqrt(rggain^2*2*l_noise^2+m_noise^2+rgnoise^2);
% Above equations checked out empirically 9/23/10

tmp = [];
for i = 1:niter
    l = X+normrnd(0,l_noise,size(X));
    m = Y+normrnd(0,m_noise,size(Y));

    lum = lumgain*(l+m)+normrnd(0,lumnoise,size(X));
    rg = rggain*(l-m)+normrnd(0,rgnoise,size(X));
    det = abs(lum) > s_lum | abs(rg) > s_rg;
    %det = (abs(lum)./s_lum).^2+(abs(rg)./s_rg).^2 > 4;  % arbitrary
    data = data+det;
    
    tmp = [tmp; lum(50,50), rg(50,50)];
end

% predicted correlation
denom = sqrt(((l_noise^2+m_noise^2)+(rgnoise^2/rggain^2))*((l_noise^2+m_noise^2)+(lumnoise^2/lumgain^2)));
%(l_noise^2-m_noise^2)/(l_noise^2+m_noise^2)
(l_noise^2-m_noise^2)/denom

plot(tmp(:,1),tmp(:,2),'k.')
corr(tmp)

figure;
imagesc(data>.82*niter)
axis xy;
colormap(gray);
axis square

%%
% Section 4
% Trying to recapitulate the model that proposed by Petri and Fred in their
% manuscript "Nonlinear integration of single photon responses in the inner
% retina sets absolute visual threshold". Evidently, a threshold
% nonlinearity should be sufficient to cause two points of inflection in
% the mean output curve (as input is increased).

x = [0:.01:20];
y = [];
THRESH = 5;
ints = [0:1:4*x(end)];
L = ints < THRESH;
for i = 1:length(x)
    tmppdf = poisspdf(ints,x(i));
    zeroresp = sum(tmppdf(L));
    tmppdf(L) = 0;
    tmppdf(1) = zeroresp
    y(i) = mean(ints*tmppdf');
end
plot(x,y,'k.');
set(gca,'Yscale','log','Xscale','log')