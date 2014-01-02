function mdb = initMouseDB(option)

global GL_DOCUPATH

% tell the user what's happening
fprintf('Initializing the mouseView data base\n')

% deal with the optional input
if ~exist('option', 'var')
    option = '';
end
OVERWRITE = strcmpi(option, 'new');

% grab all the excel files
presDir = pwd; % cd back to this at the end of the routing...
cd(GL_DOCUPATH)
d = dir;

% determine if the database is already available
WBnames = {d(:).name};
DBavailable = any(strcmpi('mouseDB.mat', WBnames));
if DBavailable && ~OVERWRITE
    load('mouseDB.mat');
else
    mdb = []; % 'mouse data base' structure
end

% figure out which files are absent from the existing MDB
if numel(mdb) == 0
    l_new = 1:numel(WBnames);
else
    mdb_names = {mdb.name};
    l_new = cellfun(@(x,y) regexpi(x,y), WBnames, repmat({mdb_names}, 1, numel(WBnames)), 'uniformoutput', false);
    l_new = cellfun(@(x) ~any(cell2mat(x)), l_new); %notice the ~ ('not')!!
    l_new = find(l_new);
end


% figure out which files in the DOCUBASE are newer than the ones in
% matlab's Mouse Database


for a = l_new
    
    % exclude hidden files (things that begin with '.'), the excel template
    % (begins with '--') or the mouse data base ('mouseDB.mat')
    exclude = any(cell2mat(regexp(d(a).name, {'^[\.]', '^[-*]', 'mouseDB.mat'})));
    if exclude; continue; end
    
    % let the user know what you're doing
    fprintf('  Updating %s\n', d(a).name);
    
    % where in the mouse db should the new stuff go?
    idx = numel(mdb)+1;
    
    % pull out the name
    name = regexp(d(a).name, '\.\w+', 'split');
    mdb(idx).name = name{1};
    
    % pull out the histology information
    mdb(idx).histo = getHistologyInfo(d(a).name);
    
    
    
end

% save the database!
save([GL_DOCUPATH, 'mouseDB.mat'], 'mdb')
fprintf('Initialization complete\n')
cd(presDir)
end


function histo = getHistologyInfo(fName)
    
    % import the data from excel
    [~, ~, raw] = xlsread(fName, 'Histology');
    
    % define the keywords, and where the information is relative to the
    % keywords.
    % [<keyword>,<field name>,<data type>,<column for actual data>]
    keys = {'Notes on histology', 'notes', 'string', 2;...
            'Slice Thickness', 'thickness', 'num',  2;...
            'Fix method', 'method', 'string', 2;...
            'Number of plates', 'numPlates', 'num', 2;...
            'Slices per plate', 'slicesPerPlate', 'num', 2;...
            'Location of physical slides', 'plateLocation', 'string', 2};
    
    for a = 1:size(keys,1)
        % find the row where the keyword appears, and then grab the actual
        % data from the appropriate column
        row = regexp({raw{:,1}}', keys{a,1}, 'ignorecase'); % returns a cell array
        row = ~cellfun(@isempty, row); % now a logical vector (notice the 'not')
        col = keys{a,4};
        
        % pull out the raw info, and assign it to the apprpriate field in
        % the structure
        tmp = raw{row, col};
        switch keys{a,3}
            case 'string'
                if isnumeric(tmp)
                    tmp = num2str(tmp);
                end
            case 'num'
                if ischar(tmp)
                    tmp = str2num(tmp);
                end
        end
        histo.(keys{a,2}) = tmp;
        
    end
    
end
