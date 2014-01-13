function [weights_gun, weights_cone, btStrp_gun, btStrp_cone, bs_diffs] = ConeComparisons;
%This script should batch process experiments with interleaved gun and cone
%noise and plot cone weights from each experiment type. It will also test
%(on a cell by cell basis) for significant differences between the two cone
%weight estimates. Generation of cone weights was taken from Greg's
%ConeAnalysis code.

%*****   CONSTANTS   ******%
NSTRAPS = 31;

cd('C:\NexFiles\Kali') %go to the directory with all the nex files
fnames = fnamesFromTxt('ConeVGun.txt'); %load a matrix of names to the expts

%for development purposes
%fnames = ['K032608005.nex'];
numExpts = size(fnames, 1);

%loop over the experiments and compile the cone weights
[weights_gun, weights_cone] = deal(nan(3, numExpts)); %just initializing
for a = 1:numExpts;
    stro = nex2stro(fnames(a, :));
    disp(stro.sum.fileName)
    
    %determine a few indicies and reconstruct approptiate calibration data
    noisetypeidx = find(strcmp(stro.sum.trialFields(1,:),'noise_type'));
    sigmaidxs = strmatch('sigma',stro.sum.trialFields(1,:));
    [M, bkgndlms] = getCalsData();
    maxT = 8;
    nstixperside = stro.sum.exptParams.nstixperside;
    
    %get the regular cone weights
    [STAs_cone, STAs_gun] = getSTAs(stro);
    [STAgunmat, STAconemat] = getExptMats();
    [weights_gun(:,a), weights_cone(:,a)] = getWeights(stro);
    
    %now find the bootstrap weights
    for n = 1:NSTRAPS;
        btStrpStro  = resampleStro();
        [STAs_cone, STAs_gun] = getSTAs(btStrpStro);
        [STAgunmat, STAconemat] = getExptMats();
        [btStrp_gun(:, n, a), btStrp_cone(:, n, a)] = getWeights(btStrpStro);
    end
    
    %now make all pairwise comparisons
    inds = [fullfact([NSTRAPS, NSTRAPS])]';
    bs_diffs(:,:,a) = btStrp_gun(:,inds(1,:), a) - btStrp_cone(:, inds(2,:), a);
end
        


    function [gunWeight, coneWeight] = getWeights(inStro);
        %needed for scaling??
        conesigmas = inStro.trial(inStro.trial(:,noisetypeidx) == 2, sigmaidxs)/1000;  % cone excitation units
        lmscontrasts = mode(conesigmas)./bkgndlms';
        
        % Finding whichpixels are different from bkgnd
        sdthresh = 3;
        whichpixgun = logical(sum(STAgunmat.^2)>sdthresh*std(sum(STAgunmat.^2)));
        whichpixcone = logical(sum(STAconemat.^2)>sdthresh*std(sum(STAconemat.^2)));
        whichpix = find(whichpixcone | whichpixgun);

        % First, the gun noise STA
        [u,s,v] = svd(STAgunmat(:,whichpix));
        if (sum(v(:,1)) < 0)
            u = -u;
        end
        oldv = v(:,1);
        coneweights = inv(M')*u(:,1);
        %coneweights = coneweights.*bkgndlms;   % Cone specific adaptation (on the weights?)
        %coneweights = coneweights.*(M*[1 1 1]');
        gunWeight = coneweights./sum(abs(coneweights));

        % Cone weights from cone noise expressed in contrast units
        [u,s,v] = svd(STAconemat(:,whichpix)./repmat(lmscontrasts',[1 length(whichpix)]));%
        [u,s,v] = svd(STAconemat(:,whichpix)./repmat(mode(conesigmas)',[1 length(whichpix)]));
        if (v(:,1)'*oldv< 0 | sum(oldv) < 0)
            u = -u;
        end
        coneWeight = u(:,1)./sum(abs(u(:,1)));
    end

    function [gunMat, coneMat] = getExptMats();
        %deal with the cone noise trials
        tmpSTA = reshape(STAs_cone, [nstixperside^2 3 maxT]);
        tmpSTA = permute(tmpSTA, [2 1 3]);
        coneMat = reshape(tmpSTA,[3 nstixperside^2*maxT]);


        %deal with the gun noise trials
        tmpSTA = reshape(STAs_gun, [nstixperside^2 3 maxT]);
        tmpSTA = permute(tmpSTA, [2 1 3]);
        gunMat = reshape(tmpSTA,[3 nstixperside^2*maxT]);
    end

    function [M, bkgndlms] = getCalsData();
        %reconstructing the gammatable
        gammaTable = stro.sum.exptParams.gamma_table;
        gammaTable = reshape(gammaTable, length(gammaTable)/3, 3);
        gammaTable1 = interp1(linspace(0,255,256),gammaTable,linspace(0,255,65536), 'spline');
        
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
    end

    function [STAs_cone, STAs_gun] = getSTAs(inStro);
       % Getting STAs
        for noisetype = 1:2
            tmpstro = inStro;
            L = inStro.trial(:,noisetypeidx) == noisetype;
            tmpstro.ras(~L,:) = [];
            tmpstro.trial(~L,:) = [];
            out = getWhtnsStats(tmpstro,maxT,'STCOVmex', {nstixperside^2, 3, maxT});
            if (noisetype == 1)
                STAs_gun = out{1};
                STCs_gun = out{2};
                nspikes_gun = out{3};
            elseif (noisetype == 2)
                STAs_cone = out{1};
                STCs_cone = out{2};
                nspikes_cone = out{3};        
            end
        end
        tmpstro = []; 
    end
        
    function btStrpStro  = resampleStro();
        btStrpStro = stro;
        btStrpStro.trial = [];
        btStrpStro.ras = {};
        
        %resample with replacement
        Lgun = find(stro.trial(:,noisetypeidx) == 1);
        gunInds = Lgun(unidrnd(length(Lgun), length(Lgun), 1));
        Lcone = find(stro.trial(:,noisetypeidx) == 2);
        coneInds = Lcone(unidrnd(length(Lcone), length(Lcone), 1));
        
        %remake the .trial and .ras arrays in the btStrpStro structure
        btStrpStro.trial = [btStrpStro.trial; stro.trial(gunInds,:); stro.trial(coneInds,:)];
        btStrpStro.ras = [btStrpStro.ras; stro.ras(gunInds,:); stro.ras(coneInds,:)];
    end
        
        
end