function mdb = initMouseDB(overwrite, suppressOutput)

%
% INITALIZE THE MOUSE DATABASE
%
%      mdb = initMouseDB(overwrite, suppressOutput)
%
%  new MDB      ->   mdb = initMouseDB('new');
%  update only  ->   mdb = initMouseDB('update');
%  non-verbose  ->   mdb = initMouseDB(~, 'notext');
%
%
% The variables [overwrite] and [suppressOutput] can also be logicals.
%
% C.Hass 2014
 
    global GL_DOCUPATH

    % deal with the optional input #1
    if ~exist('overwrite', 'var')
        overwrite = false;
    else
        switch overwrite
            case {'new', true}
                overwrite = true;
            case {'update', false}
                overwrite = false;
        end
    end
    
    % deal with optional input #2
    if ~exist('suppressOutput', 'var')
        suppressOutput = false;
    else
        switch suppressOutput
            case {true, 'notext'}
                suppressOutput = true;
            case {false, 'verbose'}
                suppressOutput = false;
        end
    end
    
    % tell the user what's happening
    if ~suppressOutput
        fprintf('Initializing the mouseView data base\n')
    end

    % grab all the excel files
    presDir = pwd; % cd back to this at the end of the routine...
    cd([GL_DOCUPATH, 'Mouse_workbooks'])
    d = dir;
    
    % determine if the database is already available
    WBnames = {d(:).name};
    DBavailable = any(strcmpi('mouseDB.mat', WBnames));
    if DBavailable && ~overwrite
        load('mouseDB.mat');
    else
        mdb = mdbobj; % 'mouse data base' object
    end
    
    % eliminate the hidden files (and other files that we should ignore)
    % from the directory tree.
    exclude = {'^[\.]', '^[-*]', '^[~]', '^[$]', 'mouseDB.mat', 'SyncToy'};
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
    if ~DBavailable || overwrite
        l_absent = true(numel(WBnames),1);
    else
        % pull out the names of things in the DB
        mdb_names = mdb.mouseNames;
        
        % find the Workbooks that are already present in the database
        l_present = cellfun(@(x,y) regexpi(x,y), WBnames, repmat({mdb_names}, 1, numel(WBnames)), 'uniformoutput', false);
        l_absent = cellfun(@(x) ~any(cell2mat(x)), l_present); %notice the ~ ('not')!!
    end


    % add the missing files
    if ~suppressOutput
        fprintf('  Adding %d files:\n', sum(l_absent));
    end
    for a = find(l_absent(:))'
        
        % display the name of the file about to get added
        if ~ suppressOutput
            fprintf('    %s\n', d(a).name)
        end

        % where in the mouse db should the new stuff go?
        if (overwrite || ~DBavailable) && (a==1)
            idx = 1;
        else
            idx = numel(mdb.mice)+1;
        end
        
        
        % build mdb entry. I'm using set field b/c the mdbobj construct
        % defines things as vectors at first, and matlab complains about
        % conversion from struct to double if I don't use setfield.
        mdb = setfield(mdb, 'mice',{idx},  {build_DB_entry(d(a))});
        mdb = setfield(mdb, 'searchText', {idx}, {build_DB_text(mdb.mice{idx})});

    end
    
    
    %
    % UPDATE FILES IN THE MDB THAT HAVE NEWER VERSIONS ON THE HARDDRIVE.
    %
    %%%%%%%%%%%%%%%%%%%%
    if ~suppressOutput
        fprintf('  Updating modified files:\n')
    end
    
    WBnames = cellfun(@(x) x(1:regexpi(x, '.\.', 'once')), WBnames, 'uniformoutput', false); % removing the file extension
    for a = 1:numel(mdb.mice)
        
        % find the MDB entry in the file directory
        idx = find(strcmpi(mdb.mice{a}.name, WBnames));
        if numel(idx) == 0;
            error('Could not locate file <%s>', mdb.mice{a}.name)
        elseif numel(idx) > 1;
            error('Too many copies of <%s> in the database', mdb.mice{a}.name)
        end
        
        % update the mdb if the version in the directory is newer
        if d(idx).datenum > mdb.mice{a}.modDate
            fprintf('    %s\n', d(idx).name)
            mdb = setfield(mdb, 'mice',{a},  {build_DB_entry(d(idx))});
            mdb = setfield(mdb, 'searchText', {a}, {build_DB_text(mdb.mice{a})});
        end
    end
    

    
    %
    % SAVE THE DATABASE
    %
    %%%%%%%%%%%%%%%%%%%%
    save([GL_DOCUPATH, filesep, 'Mouse_workbooks', filesep , 'mouseDB.mat'], 'mdb')
    if ~suppressOutput
        fprintf('Initialization complete\n')
    end
    cd(presDir)
    
    
end

function out = build_DB_entry(fileInfo)
        
        % pull out the name and the date the file was modified.
        name = regexpi(fileInfo.name, '\.\w+', 'split');
        out.name = name{1};
        out.modDate = fileInfo.datenum;

        % pull out the histology information
        out.info = getGeneralInfo(fileInfo.name);
        out.phys = getPhysInfo(fileInfo.name);
        out.histo = getHistologyInfo(fileInfo.name);
        

end

function histo = getHistologyInfo(fName)
    
    % import the data from excel
    [~, ~, raw] = xlsread(fName, 'Histology', 'A1:B8');
    
    % force everything to be strings. Technically I could use the 'txt'
    % output of xlsread, but datestrings do not come across...
    raw = cellfun(@num2str, raw, 'uniformoutput', false);
    

    % define the keywords, and where the information is relative to the
    % keywords.
    % [<keyword>,<field name>,<data type>,<column for actual data>]
    keys = {'Notes on hist',    'notes',            'string',   2;...
            'Slice Thick',      'thickness',        'num',      2;...
            'Fix method',       'method',           'string',   2;...
            'Number of plates', 'numPlates',        'num',      2;...
            'Slices per plate', 'slicesPerPlate',   'num',      2;...
            'Location of ',     'plateLocation',    'string',   2};


    histo = assignFields(keys, raw);

end

function info = getGeneralInfo(fName)


    % import the data from excel and chuck the unused cells
    [~, ~, raw] = xlsread(fName, 'General Info', 'A1:B7');
    
    % convert everyting to text.
    raw = cellfun(@num2str, raw, 'uniformoutput', false); % force everything to be strings

    % define the keywords, and where the information is relative to the
    % keywords.
    % [<keyword>,<field name>,<data type>,<column for actual data>]
    keys = {'Mouse ID',     'mouseID',              'string', 2;...
            'Male/Female',  'sex',                  'string', 2;...
            'DOB',          'dob',                  'string', 2;...
            'Strain',       'strain',               'string', 2;...
            'Identifier',   'physicalIdentifier',   'string', 2;...
            'General',      'generalNotes',         'string', 2};

    info = assignFields(keys, raw);
    
    % make sure the dob is correct
    info.dob = convertDate(info.dob);
    

end

function phys = getPhysInfo(fName)
    
    % import the data from excel and chuck the unused cells
    [~, ~, raw] = xlsread(fName, 'Physiology', 'A1:H100');
    l_nan = cellfun(@(x) any(isnan(x)), raw(:,1));
    raw(l_nan,:) = []; % cut the empty rows
    
    % occasionally the number of columns is outrageous. Cut them all.
    if size(raw,2)>12
        raw(:,13:end) = [];
    end
    
    
    % convert everyting to text.
    raw = cellfun(@num2str, raw, 'uniformoutput', false); % force everything to be strings
    
    % Compile the data
    phys.int.type = raw{1,3};
    phys.int.osm = raw{1,8};
    phys.acsf.type = raw{2,3};
    phys.acsf.osm = raw{2,8};
    phys.temp = (raw{3,3});
    phys.notes = raw{4,3};
    
    % pull out the files and notes partaining to each cell/file
    tmp = regexpi({raw{:,1}}', 'Cell #');
    firstFileIdx = ~cellfun(@isempty, tmp); % now a logical vector (notice the 'not');
    firstFileIdx = find(firstFileIdx, 1, 'first') + 1;
    keys = raw(firstFileIdx-1,:);
    for a = 1:numel(keys) % remove the spaces
        keys{a}(keys{a}==' ')=[];
    end
    key_col = find(~strcmpi(keys, 'nan'));
    
    % assign the values for each cell to the mouse database
    neuronNums = cellfun(@str2double, raw(firstFileIdx:end, 1));
    uniqueNeurons = unique(neuronNums);
    for a = 1:numel(uniqueNeurons);
        l_cell = find(neuronNums == uniqueNeurons(a)) + firstFileIdx-1; % index to rows after firstFileIdx
        
        for i = 1:numel(l_cell)
            for k = 2:numel(key_col) % skip the first keyword
                phys.cell(a).file(i).(keys{k}) = raw{l_cell(i), key_col(k)};
            end
        end
    end
    
    
end

function out = assignFields(keys, xls)

    for a = 1:size(keys,1)
        % find the row where the keyword appears, and then grab the actual
        % data from the appropriate column
        row = regexpi({xls{:,1}}', keys{a,1}, 'ignorecase'); % returns a cell array
        row = ~cellfun(@isempty, row); % now a logical vector (notice the 'not')
        col = keys{a,4};
        
        
        % pull out the raw info, and assign it to the apprpriate field in
        % the structure
        tmp = xls{row, col};
        if strcmpi(keys{a,3}, 'num')
            tmp = str2num(tmp);
        end
        out.(keys{a,2}) = tmp;
        
    end

end

function newDate = convertDate(oldDate)

    % don't do anything if the date wasn't already specified.
    if any(strcmpi(oldDate, {'nan', 'mm/dd/yyyy'}))
        newDate = '';
        return
    end
    
    % deal with the date, which could be tricky depending on how flaky
    % Excel is and how the user inputs the date.
    seps = regexp(oldDate, '/|\.|-');
    if any(seps)
        
        % parse the unique parts of the date string.
        seps = [1, seps, numel(oldDate)];
        assert(~isempty(seps), 'Date string improperly specified')
        for a = 1:numel(seps)-1
            tmp{a} = oldDate(seps(a):seps(a+1));
            tmp{a} = regexprep(tmp{a}, '/|\.|-', '');
        end
        
        % figure out what to do with the parsed args
        if (str2double(tmp{1}) > 12) && (str2double(tmp{2}) <= 12)
            yy = tmp{1};
            mm = tmp{2};
            dd = tmp{3};
        elseif (str2double(tmp{1}) <= 12) && (str2double(tmp{3}) > 2012)
            mm = tmp{1};
            dd = tmp{2};
            yy = tmp{3};
        end
        
        % make sure YY is 4 digit
        if numel(yy)<4;
            yy = ['20', yy];
        end
        
        % store the value
        newDate = datestr([mm, '/', dd, '/', yy], 'mm/dd/yyyy');
        
        
    elseif numel(oldDate)<=5
        
        % excel on the PC specifies datenums from 12/30/1899. So add an
        % offset to bring things into alignment with the mac or linux.
        offset = datenum('12/30/1899');
        tmp = str2double(oldDate)+offset;
        newDate = datestr(tmp, 23); % for some reason, things seem to be better if I subtract one?!?!

    end
end


function out = build_DB_text(in)
    
    % concatenate the easy ones
    out = [in.name, ' ',...
           in.info.generalNotes, ' ',...
           in.info.sex, ' ',...
           in.info.strain, ' ',...
           in.histo.notes, ' ',...
           in.phys.notes, ' ',...
           in.phys.acsf.type, ' ',...
           in.phys.temp, ' ',...
           in.phys.int.type, ' ',...
           in.histo.notes, ' ',...
           in.histo.method, ' ',...
           ];
    
    % now loop over the phys files pulling out the notes and stuffing them
    % in with the others.
    if isfield(in.phys, 'cell')
        for c = 1:numel(in.phys.cell)
            for f = 1:numel(in.phys.cell(c).file)
                out = [out, ' ', in.phys.cell(c).file(f).Notes, ' ', in.phys.cell(c).file(f).FileName];
            end            
        end
    end
           
end


