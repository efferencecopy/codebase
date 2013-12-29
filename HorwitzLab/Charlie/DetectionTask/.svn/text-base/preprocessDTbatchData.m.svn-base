% using this script to unpack a batch data set. I call it during every cell
% of analysis in the script "DTbatchAnalysis". I'm doing this so that each
% new analysis cell can define variables that won't get confused with other
% cells' variables...



%load the batch files.
load(cardVsIntBatchPath)
if exist('blpBatchPath', 'var')
    load(blpBatchPath)
end

%iterate over expts and compile the necessary data
for a = 1:length(out.dat);
    rawTRs(a,:) = out.dat(a).c.alpha(:)' ./ out.dat(a).m.alpha(:)'; %transpose so that colors go across rows.
    rawNTs(a,:) = out.dat(a).c.alpha(:)';
    rawPTs(a,:) = out.dat(a).m.alpha(:)';
    rawCardTRs(a,1) = out.dat(a).cardTR;
    rawIntTRs(a,1) = out.dat(a).intTR;
    rawPrefTRs(a,1)= out.dat(a).prefTR;
    rawAuxTRs(a,1) = out.dat(a).auxTR;
    rawSF(a,1) = out.dat(a).expt.sfs;
    rawModRatios(a,1) = out.dat(a).grating.modulationratio;
    rawEccentricity(a,1) = out.dat(a).grating.eccentricity;
    rawCSIs(a,1) = out.dat(a).csi;
    rawJSCSIs(a,1) = out.dat(a).jscsi;
    rawOrients(a,1) = out.dat(a).grating.orient.preforient;
    prefCards(a,:) = out.dat(a).prefCard;
    prefInts(a,:) = out.dat(a).prefInt;
    prefIsolum(a,:) = out.dat(a).prefIsolum;
    monkInitials(a,1) = out.fnames{a}{1}(1);
end

%retrieve the indices to elements in the 'errors' matrix.
neuroThresh1Ind = strmatch('nanNeuroThresh1', out.errorTypes);
neuroThresh2Ind = strmatch('nanNeuroThresh2', out.errorTypes);
psyThresh1Ind = strmatch('nanPsyThresh1', out.errorTypes);
psyThresh2Ind = strmatch('nanPsyThresh2', out.errorTypes);
lt16DTtrialsInd = strmatch('<16DTtrials', out.errorTypes);
lt20DTtrialsInd = strmatch('<20DTtrials', out.errorTypes);
infCSIInd = strmatch('infCSI', out.errorTypes);
lt5GTtrialsInd = strmatch('<5GTtrials', out.errorTypes);
lt5spikesInd = strmatch('<5sp/sec', out.errorTypes);
prefTieInd = strmatch('prefTie', out.errorTypes);
cardTieInd = strmatch('cardTie', out.errorTypes);
intTieInd = strmatch('intTie', out.errorTypes);
cardMismatchInd = strmatch('cardMismatch', out.errorTypes);
intMismatchInd = strmatch('intMismatch', out.errorTypes);
orientMismatchInd = strmatch('orientMismatch', out.errorTypes);
sfMismatchInd = strmatch('sfMismatch', out.errorTypes);
posMismatchInd = strmatch('posMismatch', out.errorTypes);
contMismatchInd = strcmp('contMismatch', out.errorTypes);
tenDegFundInd = strmatch('DT10degFund', out.errorTypes);


% calculate a few other things
commonExclusions = sum(out.errors(:,[lt16DTtrialsInd, orientMismatchInd, posMismatchInd]),2)>0;
l_sedna = monkInitials == 'S';
l_kali = monkInitials == 'K';

