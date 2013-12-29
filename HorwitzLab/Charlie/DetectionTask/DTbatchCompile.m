function out = DTbatchCompile


%retrieve the file names.
if ispc
    cd ('C:\NO BACKUP\NexFiles\Charlie')
elseif ismac
    cd('/Users/charliehass/LabStuff/NexFiles/Charlie')
end
[fname, fpath] = uigetfile('*.txt', 'Pick a text file');
[out.fnames out.spikenums] = fnamesFromTxt2([fpath, fname]);
cd (fpath); %this is where the batch data will be saved

%open up the 10deg fundamentals to aid in identifying expts that use the
%10deg versions of SMJ
load ('T_cones_smj10')
fund10deg = T_cones_smj10;


% maintain a matrix of error types. Change the zero entries to ones for
% every infraction that occurs.
out.errorTypes = {'nanNeuroThresh1', 'nanNeuroThresh2', 'nanPsyThresh1',...
    'nanPsyThresh2', '<16DTtrials', '<20DTtrials', 'infCSI', '<5GTtrials',...
    '<5sp/sec', 'cardMismatch', 'intMismatch', 'prefTie', 'cardTie',...
    'intTie', 'orientMismatch', 'sfMismatch', 'posMismatch', 'contMismatch',...
    'GT10degFund', 'DT10degFund'};
neuroThresh1Ind = strcmp('nanNeuroThresh1', out.errorTypes);
neuroThresh2Ind = strcmp('nanNeuroThresh2', out.errorTypes);
psyThresh1Ind = strcmp('nanPsyThresh1', out.errorTypes);
psyThresh2Ind = strcmp('nanPsyThresh2', out.errorTypes);
lt16DTtrialsInd = strcmp('<16DTtrials', out.errorTypes);
lt20DTtrialsInd = strcmp('<20DTtrials', out.errorTypes);
infCSIInd = strcmp('infCSI', out.errorTypes);
lt5GTtrialsInd = strcmp('<5GTtrials', out.errorTypes);
lt5spikesInd = strcmp('<5sp/sec', out.errorTypes);
cardMismatchInd = strcmp('cardMismatch', out.errorTypes);
intMismatchInd = strcmp('intMismatch', out.errorTypes);
prefTieInd = strcmp('prefTie', out.errorTypes);
cardTieInd = strcmp('cardTie', out.errorTypes);
intTieInd = strcmp('intTie', out.errorTypes);
orientMismatchInd = strcmp('orientMismatch', out.errorTypes);
sfMismatchInd = strcmp('sfMismatch', out.errorTypes);
rfPosMismatchInd = strcmp('posMismatch', out.errorTypes);
contMismatchInd = strcmp('contMismatch', out.errorTypes);
GT10degFundInd = strcmp('GT10degFund', out.errorTypes);
DT10degFundInd = strcmp('DT10degFund', out.errorTypes);
out.errors(1:size(out.fnames, 1), 1:size(out.errorTypes,2)) = 0;


%common analysis params for all expts:
defaultParams.cellNum = 1;          % forces the DT and GT analysis to the first unit
defaultParams.start = 'Gabor On';        %from gabor onset
defaultParams.end = 'Gabor Off';          %from gabor offset
defaultParams.lowCutoff = 3;        %for choice probability
defaultParams.meth = 'rate';        % gets updated on the basis of mod ratio...

%iterate over the data files.
for a = 1:size(out.fnames,1);
    fprintf('File <%d> out of <%d>: %s\n',a, length(out.fnames), out.fnames{a}{1});
    
    GT = gtobj(out.fnames{a}{2});
    [out.dat(a).grating] = getGratingTuning(GT, defaultParams.cellNum);

    DT = dtobj(out.fnames{a}{1});
    if out.dat(a).grating.modulationratio > 1;
        defaultParams.meth = 'F1'; %i.e, f1 resp amplitude
    elseif out.dat(a).grating.modulationratio < 1 || isnan(out.dat(a).grating.modulationratio)
        defaultParams.meth = 'rate'; %i.e., firing rate
    end
    [out.dat(a).m, out.dat(a).c, out.dat(a).expt] = DTmocsUnpack(DT, defaultParams);
    out.dat(a).analysisParams = defaultParams;
    close all %don't accumulate the figures that DTunpack prints out
    
    
    
    %
    % Preferred colors: intermediate and cardinal
    %
    % standardize the color directions so that the L cone sign is
    % always positive. The S-iso case is always all cones positive (for
    % the grating paradigm) so we don't need to worry about the
    % ambiguous case here (like in DTmocsUnpack)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    colsFromGT = out.dat(a).grating.color.colors;
    norms = sqrt(sum(colsFromGT.^2,2));
    colsFromGT = colsFromGT./repmat(norms, 1, 3);
    colsFromGT(colsFromGT(:,1)<0, :) = colsFromGT(colsFromGT(:,1)<0,:) .* -1; %standardize the color directions
    GTcolorSigns = sign(colsFromGT);
    l_Siso = ismember(GTcolorSigns, [0 0 1], 'rows');
    l_LvM = ismember(GTcolorSigns, [1 -1 0], 'rows');
    l_SwM = ismember(GTcolorSigns, [1 -1 -1], 'rows');
    l_SwL = ismember(GTcolorSigns, [1 -1 1], 'rows');
    l_LM = ismember(GTcolorSigns, [1 1 0], 'rows');
    
    % pref Isolum Color
    l_isolum = l_Siso | l_LvM | l_SwM | l_SwL;
    rate = max(out.dat(a).grating.color.colresp(l_isolum, 1)); %max resp to isolum color
    idx = (out.dat(a).grating.color.colresp(:,1) == rate) & l_isolum;
    if sum(idx)>1
        out.errors(a, prefTieInd) = 1;
        out.dat(a).prefIsolum = [nan nan nan];
    else
        out.dat(a).prefIsolum = colsFromGT(idx,:);
    end
    
    %pref cardinal
    l_cardinal = l_Siso | l_LvM;
    rate = max(out.dat(a).grating.color.colresp(l_cardinal,1));
    idx = (out.dat(a).grating.color.colresp(:,1) == rate) & l_cardinal;
    if sum(idx)>1
        out.errors(a, cardTieInd) = 1;
        out.dat(a).prefCard = [nan nan nan];
    else
        out.dat(a).prefCard = colsFromGT(idx,:);
    end
    
    %pref intermediate
    l_intermediate = l_SwM | l_SwL;
    rate = max(out.dat(a).grating.color.colresp(l_intermediate, 1));
    idx = (out.dat(a).grating.color.colresp(:,1) == rate) & l_intermediate;
    if sum(idx)>1
        out.errors(a, intTieInd) = 1;
        out.dat(a).prefInt = [nan nan nan];
    else
        out.dat(a).prefInt = colsFromGT(idx,:);
    end
    
    %
    %     CSI
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    isoLumRate = max(out.dat(a).grating.color.colresp(l_isolum, 1)); %max resp to isolum color
    lumRate = out.dat(a).grating.color.colresp(l_LM,1);
    LvMRate = out.dat(a).grating.color.colresp(l_LvM,1);
    out.dat(a).csi = isoLumRate./lumRate;
    out.dat(a).jscsi = LvMRate./lumRate;
    
    
    %
    %   RF POS
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    out.dat(a).expt.rfpos = [DT.sum.exptParams.rf_x, DT.sum.exptParams.rf_y];
    
    %
    %  FUNDAMENTAL IN USE
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    GT_funds = reshape(GT.sum.exptParams.fundamentals, [], 3);
    if all(all(GT_funds == fund10deg')); %all(down cols), then all(across rows)
        out.errors(a, GT10degFundInd) = 1;
    end
    DT_Mmtx = reshape(DT.sum.exptParams.m_mtx, [],3);
    DT_monspect = reshape(DT.sum.exptParams.mon_spect, [], 3);
    Mmtx_assuming10degFunds = fund10deg * DT_monspect;
    if all(all(DT_Mmtx == Mmtx_assuming10degFunds))
        out.errors(a, DT10degFundInd) = 1;
    end
    
    
    %
    %   Threshold ratios
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %cardinal color
    idx = ismember(sign(out.dat(a).expt.standColors), sign(out.dat(a).prefCard), 'rows');
    switch sum(idx)
        case 0
            out.dat(a).cardTR = nan;
            if ~out.errors(a, cardTieInd)
                out.errors(a, cardMismatchInd) = 1;
            end
        case 1
            out.dat(a).cardTR = out.dat(a).c.alpha(idx)./out.dat(a).m.alpha(idx);
        case 2
            error('too many indicies')
    end
    
    %intermediate color
    idx = ismember(sign(out.dat(a).expt.standColors), sign(out.dat(a).prefInt), 'rows');
    switch sum(idx)
        case 0
            out.dat(a).intTR = nan;
            if ~out.errors(a,intTieInd)
                out.errors(a, intMismatchInd) = 1;
            end
        case 1
            out.dat(a).intTR = out.dat(a).c.alpha(idx)./out.dat(a).m.alpha(idx);
        case 2
            error('too many indicies')
    end
    
    %pref color
    idx = ismember(sign(out.dat(a).expt.standColors), sign(out.dat(a).prefIsolum), 'rows');
    switch sum(idx)
        case 0
            out.dat(a).prefTR = nan;
            out.dat(a).auxTR = nan;
            if ~out.errors(a,prefTieInd)
                if ~(out.errors(a, intMismatchInd) || out.errors(a, cardMismatch));
                    error('something wrong with expt <%d>', a);
                end
            end
        case 1
            out.dat(a).prefTR = out.dat(a).c.alpha(idx)./out.dat(a).m.alpha(idx);
            out.dat(a).auxTR = out.dat(a).c.alpha(~idx)./out.dat(a).m.alpha(~idx);
        case 2
            error('too many indicies')
    end
    
    %
    %   Check that the contrast ranges were set reasonably well for each color
    %   direction. All we're trying to do is measure signals at detection
    %   threshold, so make sure that the contrast ranges span detection
    %   threshold.
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cont = vertcat(out.dat(a).m.performance{:});
    contAbove = any(cont>=0.9, 2);
    contBelow = any(cont<=0.75,2);
    goodColor = contAbove & contBelow;
    if sum(goodColor) ~= 2
        out.errors(a, contMismatchInd) = 1;
    end
    
    
    % a list of additional checks for book keeping purposes
    if isinf(out.dat(a).csi); out.errors(a, infCSIInd) = 1; end %infinate CSI
    if isnan(out.dat(a).m.alpha(1)); out.errors(a, psyThresh1Ind) = 1; end %flat neuroFun
    if isnan(out.dat(a).m.alpha(2)); out.errors(a, psyThresh2Ind) = 1; end
    if isnan(out.dat(a).c.alpha(1)); out.errors(a, neuroThresh1Ind) = 1; end   %flat psyFun
    if isnan(out.dat(a).c.alpha(2)); out.errors(a, neuroThresh2Ind) = 1; end
    if min(horzcat(out.dat(a).c.nTrialsIn{:}))<8; out.errors(a, lt16DTtrialsInd) = 1; end %incomplete DT file
    if min(horzcat(out.dat(a).c.nTrialsIn{:}))<10; out.errors(a, lt20DTtrialsInd) = 1; end %incomplete DT file
    if min(out.dat(a).grating.color.colresp(:,3))<5; out.errors(a, lt5GTtrialsInd) = 1; end %fewer than 5 protocol4 trials per color
    if max(out.dat(a).grating.color.colresp(:,1))<5; out.errors(a, lt5spikesInd) = 1; end %driven response druing protocol4 is <5sp/sec
    
    diffTollerence = pi/6;
    orientDiff = abs(out.dat(a).grating.orient.preforient-unique(DT.trial(:, DT.idx.gaborTheta)));
    if orientDiff>pi; orientDiff = (2*pi) - orientDiff; end %allows for two orients near 0
    orientDiff = min(orientDiff, pi-orientDiff); %allows for non-direction selectivity
    if orientDiff > diffTollerence; out.errors(a, orientMismatchInd) = 1; end %the orientation used during GT is not matched to DT
    
    DTsf = 1/unique(DT.trial(:, DT.idx.gaborLambda)).*DT.sum.exptParams.pixperdeg;
    sfRatio = DTsf ./ out.dat(a).grating.sf.prefSF;
    if ((sfRatio < 0.50) || (sfRatio > 2)); out.errors(a, sfMismatchInd) = 1; end %DT has to be within 1/2 and 2x of GT
    
    xDiff = abs(GT.sum.exptParams.rf_x - DT.sum.exptParams.rf_x);
    yDiff = abs(GT.sum.exptParams.rf_y - DT.sum.exptParams.rf_y);
    if any([xDiff, yDiff] > 0); out.errors(a, rfPosMismatchInd) = 1; end %GT and DT use different RF locations.
end

%
%save a date stamped .mat file before going on to the analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
str = ['CardVsInt_' date];
eval(['save ' str ' out']);
presDir = pwd;
fprintf('**** Batch data saved to directory: %s ****\n', presDir);

end
