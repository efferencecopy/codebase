function out = strocat(varargin)

%   EXAMPLE:   stro = strocat(first, second, third, ...);
%              stro = strocat({first, second, third,...}); %input as a cell array
%
% concatinates STRO files after doing a small amount of error checking to
% make sure that the two files have the same experimental parameters.
% Frist, second, third are all themselves stro files (not nex files). I'll
% use the first file passed in as the standard with which to compare the
% subsequent files.
%
% CAH 9.10.08


if iscell(varargin{1})
    nFiles = numel(varargin{1});
else
    nFiles = length(varargin);
end

% use the first file as the standard for comparison
if iscell(varargin{1})
    out = varargin{1}{1};
else
    out = varargin{1};
end


%out.sum.concat = out.sum.fileName(end-14:end-4);
id = out.sum.paradigmID;
names = fieldnames(out.sum.exptParams);

%iterate over the subsequent files and concatenate
for a = 2:nFiles
    if iscell(varargin{1})
        tmp = varargin{1}{a};
    else
        tmp = varargin{a};
    end
    
    %break if the files don't have the same id
    if tmp.sum.paradigmID ~= id;
        error('One or more of the files came from different experimental paradigms');
    end
    
    %make sure that the experimental parameters are the same
    for b = 1:length(names);
        out_value  = getfield(out.sum.exptParams, names{b});
        tmp_value  = getfield(tmp.sum.exptParams, names{b});
        match = (out_value == tmp_value);
        if ~match & ~all(isnan([out_value; tmp_value]))  % in case they're both nans
            warning('\n*** File <%d> has an inconsistant <%s> parameter ***', a, names{b})
        end
    end
    %if you've made it to this point everything should check out, so
    %concatenate away.
    
    % First concatenating the filenames
    fname = integrateFilenames(out.sum.fileName, tmp.sum.fileName);
    out.sum.fileName = fname;
    out.trial = [out.trial; tmp.trial];
    if (size(out.ras,2) ~= size(tmp.ras,2))
        [out.ras, out.sum.rasterCells] = integrateRaster(out.ras, out.sum.rasterCells, tmp.ras, tmp.sum.rasterCells);
    else
        out.ras = [out.ras; tmp.ras];
    end
    out.other = [out.other; tmp.other];

    if any([210, 212] == tmp.sum.paradigmID)
        out.LMS = [out.LMS ;  tmp.LMS];
    end
end


% Nested functions

    % Concatenates filenames.  fname1 serves as the root and fname2 is used
    % to modify it.  First, we take fname2 and strip off everything before 
    % the final '\' and the suffix '.nex'.
    % Then we match succesively larger pieces of fname2 against fname1
    % until we can no longer find a match.
    % What is returned is: "<fname1>&<non-matching part of fname2>.nex"
    % so the character "&" denotes the concatenation.
    function fname = integrateFilenames(fname1, fname2)
        slashidx = find(fname2 == filesep, 1, 'last' );
        if ~isempty(slashidx)
            tmpname = fname2(slashidx+1:end);
        end
        % stripping off optional '.nex' from the name of the 2nd file.
        if (length(tmpname) > 4 && strcmp(tmpname(end-3:end),'.nex'))
            tmpname(end-3:end) = [];
        end
        
        for c = 1:length(tmpname)
            idx = strfind(fname1, tmpname(1:c));
            if isempty(idx)
                continue;
            end
        end
        bittoappend = tmpname(c:end);
        if (strcmp(fname1(end-3:end),'.nex'));
            fname = [fname1(1:end-4),'&',bittoappend,'.nex'];
        else
            fname = [fname1,'&',bittoappend];
        end
    end

    % Concatenates raster matrices in the event that they are unequally
    % sized.  This happens when a cell is picked up or lost between data
    % files.
    function [ras, cells] = integrateRaster(ras1, cells1, ras2, cells2)
        ras = {};
        cells = {};
        % First iterating over the columns of the first file
        for i = 1:length(cells1)
            cells(i) = cells1(i);
            idx = find(ismember(cells2,cells1(i)));
            if (~isempty(idx))
                ras(:,i) = [ras1(:,i);ras2(:,idx)]; 
            else
                ras(:,i) = [ras1(:,1);cell(size(ras2,1),1)];
            end
        end
        % Now taking care of any columns in the second file that didn't
        % appear in the first file
        for i = 1:length(cells2)
            idx = find(ismember(cells2(i),cells1), 1);
            if (isempty(idx))
                tmpras = ras(:,i:end);
                tmpcell = cells(i:end);
                ras(:,i) = [cell(size(ras1,1),1);ras2(:,i)];
                cells(i) = cells2(i);
                ras(:,i+1:i+size(tmpras,2)) = tmpras;
                cells(i+1:i+length(tmpcell)) = tmpcell;
            end
        end
    end
end
  