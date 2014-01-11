%% experimenting with unpacking an abf file

fin

fileName = '13d20004';
mouse = 'EMX_3';
datPath = [GL_DATPATH, mouse, filesep, 'Physiology'];
filePath = findfile(fileName, datPath, '.abf');


% open the file

[dat, h] = my_abfload(filePath);

%
% plot all the sweeps, or the one sweep if it's gap free.
%
%%%%%%%%%%

% define the time base
sampRate = 1/(h.si.*10^(-6));
N = size(dat,1);
tt = [0:N-1]./sampRate;

sweep = 1;
figure
nplts = numel(h.recChNames);
for a = 1:nplts;
    subplot(nplts, 1, a)
    plot(tt, dat(:,a,sweep))
    ylabel(h.recChNames{a})
    if a == 1
        xlabel('Time (ms)')
    end
end


%% Make a database of excel files. 

fin

% tell the user what's happening
fprintf('Initializing the mouseView data base\n')

% grab all the excel files
cd(GL_DOCUPATH)
d = dir;

% determine if the database is already available
WBnames = {d(:).name};
DBavailable = any(strcmpi('mouseDB.mat', WBnames));
if DBavailable
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
    [num, txt, raw] = xlsread(d(a).name, 'Histology')
    
    
end

% save the database!
save([GL_DOCUPATH, 'mouseDB.mat'], 'mdb')
fprintf('Initialization complete\n')


%% NOTES

% ONE: 
% I'd like a way to add entries to the index only when the excel files have
% changed. If I edit a file on a PC, push the changes to GitHub, and then
% pull the updates onto my Mac, than the "date changed" field may differ
% between the Mac and PC.... I'll need to check this. Otherwise I may need
% to dip into the Git change log to determine when files change....














% seems like ABF load is applying "parameter" related header info based off
% of only one recorded channel. Is it possible to output the V/I clamp
% value for each analog input channel?

