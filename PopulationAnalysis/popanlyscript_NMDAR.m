%% IMPORT THE NECESSARY DATA INTO A datDB STRUCTURE

fin

% buid a structure of params from physiology notes. The globals are used
% during execution of the physiology_notes script
global GL_ADD_TO_MDB GL_SUPPRESS_ANALYSIS
GL_ADD_TO_MDB = true; %#ok<NASGU>
GL_SUPPRESS_ANALYSIS = true; %#ok<NASGU>
physiology_notes % run to store prams in MDB.

% cd to the directory of workbooks and import the spreadsheet of cells
workbookpath = [GL_DOCUPATH, filesep, 'Other_workbooks', filesep, 'Cell_Library'];
workbooksheet = 'NMDAR';
[~, ~, raw] = xlsread(workbookpath, workbooksheet);
header = raw(1,:);
mouseNames = raw(2:end,1);
siteNumber = raw(2:end,2);
goodNeurons(:,1) = cat(1,raw{2:end, 3});
goodNeurons(:,2) = cat(1,raw{2:end, 4});
goodNeurons = logical(goodNeurons);
neuronType = raw(2:end, 6:7);
HVA = raw(2:end, 5);


%
% build a population data structure by performing the appropriate analysis.
% I'll use the params.name field to build each new entry into the DB.
% Finish by saving the structure (possibly with a date string)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initalize things that I care about
[dat.ivcurve.mV, dat.ivcurve.mv_corrected, dat.ivcurve.pA] = deal(repmat({[] []}, numel(mouseNames), 1));

% iterate over the mice in the cell library (some mice get analyzed
% multiple times if there were multiple experiments per mouse)
for ex = 1:numel(mouseNames)
    
    ex_mouseName = mouseNames{ex};
    ex_siteNum = siteNumber{ex};
    mdb = initMouseDB(false, true);
    [~, idx] = mdb.search(ex_mouseName);
    
    params = mdb.mice{idx}.popAnly{ex_siteNum};
    params.fxns = {@anlyMod_optoIV, @anlyMod_NMDAR};
    params = invitroAnalysisOverview(params);
    close all; drawnow
    
    %
    % add the data to an array for the AMPA/NMDA
    % Excitation/Inhibition ratio stuff. Adding 'NMDAR' to the 'group'
    % array allows the script to pull out NMDAR IV curve data
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for g = 1:numel(group)
        
        if isfield(params.ivdat, 'NMDAR')
            
            for ch = 1:2
                
                if ~isempty(params.ivdat.NMDAR.ivcurve.mV{ch});
                    dat.ivcurve.mV{ex, ch} = params.ivdat.NMDAR.ivcurve.mV{ch};
                    dat.ivcurve.mv_corrected{ex, ch} = params.ivdat.NMDAR.ivcurve.mV_corrected{ch};
                    dat.ivcurve.pA{ex, ch} = params.ivdat.NMDAR.ivcurve.pA{ch};
                end
                
            end
        end
    end % iterate over groups
    
end % iterate over recording sites

% package all the useful things into a single structure
dat.hva = HVA;
dat.goodNeurons = goodNeurons;
dat.neuronType = neuronType;
dat.mice = mouseNames;
dat.siteNum = cat(1, siteNumber{:});


% create grouping lists for HVAs
l_valid = dat.goodNeurons(:);
hvas = repmat(dat.hva, 2,1);
hvaList.('pm') = cellfun(@(x) ~isempty(x), regexpi(hvas, 'pm'));
hvaList.('lm') = cellfun(@(x) ~isempty(x), regexpi(hvas, 'lm'));
hvaList.('al') = cellfun(@(x) ~isempty(x), regexpi(hvas, 'al'));
hvaList.('und') = cellfun(@(x) ~isempty(x), regexpi(hvas, 'und'));

%
% create grouping lists for neuron type
neuronType = dat.neuronType(:);
neuronType = cellfun(@num2str, neuronType, 'uniformoutput', false);
l_IN = cellfun(@(x) ~isempty(x), regexpi(neuronType, 'in'));
l_SOM = cellfun(@(x) ~isempty(x), regexpi(neuronType, 'som'));
typeList.IN = l_IN | l_SOM;
typeList.SOM = l_SOM;
typeList.PY = cellfun(@(x) ~isempty(x), regexpi(neuronType, 'py'));
typeList.und = ~typeList.IN & ~typeList.PY;


% save the data and the grouping lists
originalDir = pwd;
cd(GL_POPDATPATH);
save('popAnly_NMDAR.mat', 'dat', 'hvaList', 'typeList')
cd(originalDir);

% be nice and return these variables to their default values
GL_ADD_TO_MDB = false;
GL_SUPPRESS_ANALYSIS = false;
