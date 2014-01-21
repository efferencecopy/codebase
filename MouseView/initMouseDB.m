function mdb = initMouseDB(option)
 
 % %%%%%%%%%%%%%%%%%%
 %
 % to do
 %
 % I need to make sure there is a one to one match between files in the
 % DocuBase dir and entries in the MBD. For each file there should be only
 % one entry, and vice versa. Do this error checking b/4 doing anything
 % else.
 %
 % Then add workbooks that do not exist in the directory
 %
 % Then update MDB entries that are out of date. At this point in the code,
 % the number of entries in the MDB and the dir should be identical (after
 % ignoring the excluded files).
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
        mdb = {}; % 'mouse data base' structure
    end
    
    % eliminate the hidden files (and other files that we should ignore)
    % from the directory tree.
    exclude = {'^[\.]', '^[-*]', 'mouseDB.mat'};
    l_exclude = cellfun(@(x,y) regexpi(x,y), {d.name}, repmat({exclude}, 1, numel(d)), 'uniformoutput', false);
    l_exclude = cellfun(@(x) any(cell2mat(x)), l_exclude);
    d(l_exclude) = [];
    WBnames = {d(:).name}; % up-date this variable after deleting some entries from the dir tree
    
    % figure out which files are absent from the existing MDB
    if ~DBavailable || OVERWRITE
        l_absent = true(numel(WBnames),1);
    else
        % pull out the names of things in the DB
        disp('Need to deal with this section')
        keyboard
        mdb_names = {mdb{:}.name};
        
        % find the Workbooks that are already present in the database
        l_present = cellfun(@(x,y) regexpi(x,y), WBnames, repmat({mdb_names}, 1, numel(WBnames)), 'uniformoutput', false);
        l_absent = cellfun(@(x) ~any(cell2mat(x)), l_present); %notice the ~ ('not')!!
    end


    % add the missing files
    if any(l_absent); 
        fprintf('  Adding new files:\n');
    end
    for a = find(l_absent)'
        
        % display the name of the file about to get added
        fprintf('    %s\n', d(a).name)

        % where in the mouse db should the new stuff go?
        idx = numel(mdb)+1;
        
        
        % build mdb entry
        mdb{1,idx} = build_DB_entry(d(a)); % make a column of structures.

    end
    
    
    
    
%     % figure out which files in the DOCUBASE are newer than the ones in
%     % matlab's Mouse Database
%     if numel(mdb) == 0
%         l_new = true(numel(WBnames),1);
%     else
%         % for each entry in the DB, determine which file it corresponds to
%         % in the DocuBase directory. Then compare the modDates for the
%         % workbook and the MDB entry.
%         WBnames = cellfun(@(x) x(1:regexp(x, '.\.', 'once')), WBnames, 'uniformoutput', false); % removing the file extension
%         l_new = false(numel(WBnames), 1);
%         for a = 1:numel(mdb)
%             
%             idx = find(strcmpi(mdb_names{a}, WBnames));
%             if numel(idx) == 0; warning('Entry <%s> exists in the Mouse Database but there is no workbook in the DocuBase', mdb_names{a}); continue; end
%             if numel(idx) > 1; error('Entry <%s> from the MDB has multiple workbooks in the DocuBase', mdb_names{a}); end
%             
%             % compare the modification date
%             datestr(mdb(a).modDate)
%             datestr(d(idx).datenum)
%             l_new(a) = mdb(a).modDate ~= d(idx).datenum;
%             
%         end
%         
%     end

    % save the database!
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
