%% IMPORT THE NECESSARY DATA INTO A datDB STRUCTURE

fin

% buid a structure of params from physiology notes
global GL_ADD_TO_MDB GL_SUPPRESS_ANALYSIS
GL_ADD_TO_MDB = true;
GL_SUPPRESS_ANALYSIS = true;
physiology_notes % run to store prams in MDB.

% cd to the directory of workbooks and import the spreadsheet of cells
workbookpath = [GL_DOCUPATH, filesep, 'Other_workbooks', filesep, 'Cell_Library'];
workbooksheet = 'E_I and A_N ratios';
[~, ~, raw] = xlsread(workbookpath, workbooksheet);
header = raw(1,:);
mouseNames = raw(2:end,1);
siteNumber = raw(2:end,2);
goodNeurons(:,1) = cat(1,raw{2:end, 3});
goodNeurons(:,2) = cat(1,raw{2:end, 4});
goodNeurons = logical(goodNeurons);
HVA = raw(2:end, 5);


% build a datDB_EIAN by performing the appropriate analysis. I'll use the
% params.name field to build each new entry into the DB. Finish by saving
% the DB (possibly with a date string)
% Things to save
%  * raw conductances
%  * raw currents?
%  * peak conductances
%  * list of errors
%      * Vclamp error too big

% iterate over the mice in the cell library (some mice get analyzed
% multiple times if there were multiple experiments per mouse)
datDB_EIAN = [];
[ampa.peak, nmda.peak, excit.peak, inhib.peak] = deal(nan(numel(mouseNames), 2));
for ex = 1:numel(mouseNames)
    
    ex_mouseName = mouseNames{ex};
    ex_siteNum = siteNumber{ex};
    mdb = initMouseDB(false, true);
    [~, idx] = mdb.search(ex_mouseName);
    
    params = mdb.mice{idx}.popAnly{ex_siteNum};
    params.fxns = {@anlyMod_optoIV, @anlyMod_EIbalance};
    params = invitroAnalysisOverview(params);
    close all
    
    % add the data to an array
    idx = goodNeurons(ex,:);
    if isfield(params.isolatedData, 'ampa')
        tmp = params.isolatedData.ampa.peak_nS;
        tmp(cellfun(@isempty, tmp)) = {nan};
        tmp = cat(2,tmp{:});
        ampa.peak(ex, idx) = tmp(idx);
    end
    
    if isfield(params.isolatedData, 'nmda')
        tmp = params.isolatedData.nmda.peak_nS;
        tmp(cellfun(@isempty, tmp)) = {nan};
        tmp = cat(2,tmp{:});
        nmda.peak(ex, idx) = tmp(idx);
    end
    
    if isfield(params.isolatedData, 'excit')
        tmp = params.isolatedData.excit.peak_nS;
        tmp(cellfun(@isempty, tmp)) = {nan};
        tmp = cat(2,tmp{:});
        excit.peak(ex, idx) = tmp(idx);
    end
    
    if isfield(params.isolatedData, 'inhib')
        tmp = params.isolatedData.inhib.peak_nS;
        tmp(cellfun(@isempty, tmp)) = {nan};
        tmp = cat(2,tmp{:});
        inhib.peak(ex, idx) = tmp(idx);
    end
        
end
