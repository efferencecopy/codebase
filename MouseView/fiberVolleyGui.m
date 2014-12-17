function fiberVolleyGui


% load in the excel workbook
global GL_DOCUPATH

fname = [GL_DOCUPATH, 'Other_workbooks', filesep, 'fiberVolleyCellList.xlsx'];
[~,txt, raw] = xlsread(fname);
raw(size(txt,1)+1:end, :) = [];
raw(:,size(txt,2)+1:end) = [];
clear txt

% make the list
TFs = cellfun(@(x) num2str(x), raw(2:end,4), 'uniformoutput', false);
site = cellfun(@(x) num2str(x), raw(2:end, 2), 'uniformoutput', false);
list = cellfun(@(x,y,z) [x,'_site',y,'_',z], raw(2:end,1), site, TFs, 'uniformoutput', false);
[uniquelist, ~, udat.io.unique_idx] = unique(list, 'stable');


% initialize the gui

udat.h.main = figure;
set(udat.h.main, 'position', [1096 228 333 578]);

udat.h.list = uicontrol('style', 'listbox',...
                        'units', 'normalized',...
                        'position', [0.05, 0.1, 0.9, 0.8],...
                        'string', uniquelist,...
                        'enable', 'on',...
                        'fontsize', 14,...
                        'callback', @io_callAnalysisRoutine);

% package a few things away in the user data field.
udat.io.exptWorkbook = raw;
set(udat.h.main, 'userdata', udat);

end


function io_callAnalysisRoutine(varargin)

    persistent chk

    if isempty(chk)
        chk = 1;
        pause(0.3)
        if chk==1 % single click
            chk = [];
        end
        return
    else %double click
        chk = [];
    end

    % grab the user data
    udat = get(gcf, 'userdata');


    % determine the data files that need to be analyzed, and pass them onto the
    % downstream program.
    exptNum = get(udat.h.list, 'value');
    idx = udat.io.unique_idx == exptNum;
    idx = [false;idx]; % adding a leading zero to account for the header row.

    
    % call the main fiber volley analysis function
    fiberVolleyAnalysis(idx, udat.io.exptWorkbook);
end

