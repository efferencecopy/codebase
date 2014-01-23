function mdb = initMouseDB(option)
 
 % %%%%%%%%%%%%%%%%%%
 %
 % to do
 %
 % Add functionality to import the other work sheets.
 %
 % Start making GUI to visualize data and to text searches.
 %
 % %%%%%%%%%%%%%%%%%%%%%
 
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
    
    % eliminate the hidden files (and other files that we should ignore)
    % from the directory tree.
    exclude = {'^[\.]', '^[-*]', 'mouseDB.mat'};
    l_exclude = cellfun(@(x,y) regexpi(x,y), {d.name}, repmat({exclude}, 1, numel(d)), 'uniformoutput', false);
    l_exclude = cellfun(@(x) any(cell2mat(x)), l_exclude);
    d(l_exclude) = [];
    WBnames = {d(:).name}; % up-date this variable after deleting some entries from the dir tree
    
    
    
    %
    % MAKE SURE THERE ARE NO DUPLICATES IN THE
    % MOUSE DATA BASE
    %
    %%%%%%%%%%%%%%%%%%%%
    uniqueNames = unique(WBnames');
    assert(numel(uniqueNames) <= numel(WBnames), 'Duplicate entries exist in the MDB');
    
    
    
    %
    % ADD FILES THAT ARE ABSENT FROM THE DATABASE
    %
    %%%%%%%%%%%%%%%%%%%%
    if ~DBavailable || OVERWRITE
        l_absent = true(numel(WBnames),1);
    else
        % pull out the names of things in the DB
        mdb_names = {mdb.mice(:).name};
        
        % find the Workbooks that are already present in the database
        l_present = cellfun(@(x,y) regexpi(x,y), WBnames, repmat({mdb_names}, 1, numel(WBnames)), 'uniformoutput', false);
        l_absent = cellfun(@(x) ~any(cell2mat(x)), l_present); %notice the ~ ('not')!!
    end


    % add the missing files 
    fprintf('  Adding %d files:\n', sum(l_absent));
    for a = find(l_absent(:))'
        
        % display the name of the file about to get added
        fprintf('    %s\n', d(a).name)

        % where in the mouse db should the new stuff go?
        if (OVERWRITE || ~DBavailable) && (a==1)
            idx = 1;
        else
            idx = numel(mdb.mice)+1;
        end
        
        
        % build mdb entry
        mdb.mice(idx) = build_DB_entry(d(a)); % make a column of structures.

    end
    
    
    
    %
    % UPDATE FILES IN THE MDB THAT HAVE NEWER VERSIONS ON THE HARDDRIVE.
    %
    %%%%%%%%%%%%%%%%%%%%
    fprintf('  Updating modified files:\n') 
    WBnames = cellfun(@(x) x(1:regexp(x, '.\.', 'once')), WBnames, 'uniformoutput', false); % removing the file extension
    for a = 1:numel(mdb.mice)
        
        % find the MDB entry in the file directory
        idx = find(strcmpi(mdb.mice(a).name, WBnames));
        assert(numel(idx)==1, 'Number of matches ~= 1');
        
        % update the mdb if the version in the directory is newer        
        if d(idx).datenum > mdb.mice(a).modDate
            fprintf('    %s\n', d(idx).name)
            mdb.mice(a) = build_DB_entry(d(idx));
        end        
    end
    

    
    %
    % SAVE THE DATABASE
    %
    %%%%%%%%%%%%%%%%%%%%
    save([GL_DOCUPATH, 'mouseDB.mat'], 'mdb')
    fprintf('Initialization complete\n')
    cd(presDir)
    
    
end

function out = build_DB_entry(fileInfo)
        
        % pull out the name and the date the file was modified.
        name = regexp(fileInfo.name, '\.\w+', 'split');
        out.name = name{1};
        out.modDate = fileInfo.datenum;

        % pull out the histology information
        out.histo = getHistologyInfo(fileInfo.name);

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
                    if isnan(tmp); tmp = ''; end
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
