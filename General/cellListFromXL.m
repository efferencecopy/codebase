function [mdbidx, cellnum] = cellListFromXL(workbookname, workbooksheet, attribute)

global GL_DOCUPATH




% condition the inputs
if ~exist('workbooksheet', 'var') || isempty(workbooksheet)
    workbooksheet = 1; % take the first one.
end



% cd to the directory of workbooks and import the spreadsheet
cd([GL_DOCUPATH, filesep, 'Other_workbooks'])
[~, ~, raw] = xlsread(workbookname, workbooksheet);


% remove the first line and designate it as the header
header = raw(1,:);
l_nan = cellfun(@(x) all(isnan(x)), header, 'uniformoutput', true);
header(l_nan) = [];
raw(1,:) = [];

% iterate over the optional inputs (varargin) and find the neurons that
% satisfy these conditions.
if ~exist('attribute', 'var') || isempty(attribute)
    l_match = true(size(raw,1), 1);
else
    attributeNames = fieldnames(attribute);
    match = false(size(raw,1), numel(attributeNames));
    for a = 1:numel(attributeNames);
        
        % figure out which column to grab
        colIdx = regexpi(header, attributeNames{a});
        colIdx = cellfun(@(x) ~isempty(x), colIdx);
        
        % turn everything to text so that regexp works.
        rowText = cellfun(@num2str, raw(:,colIdx), 'uniformoutput', false);
        
        % find the rows that match the attribute
        rowIdx = regexpi(rowText, num2str(attribute.(attributeNames{a})));
        rowIdx = cellfun(@(x) ~isempty(x), rowIdx);
        
        % compile all the logical vectors
        match(:,a) = rowIdx;
    end
    
    l_match = sum(match, 2) == numel(attributeNames);
    
end


% now figure out the indicies of each neuron in the mouseDB
mdb = initMouseDB;
mdbidx = nan(sum(l_match), 1);
cellnum = nan(sum(l_match), 1);
hits = find(l_match);
for a = 1:sum(l_match)
    
    % mdb idx
    [~, tmp] = mdb.search(raw{hits(a),1});
    assert(sum(tmp)==1, 'ERROR: found zero or multiple instances of mouse name')
    mdbidx(a) = find(tmp);
    
    
    % cell num
    cellnum(a) = raw{hits(a), 2};
    
end
    
    


