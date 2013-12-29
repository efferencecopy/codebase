function [fnames, spikeIdx] = fnamesFromTxt(txtFile, spikeID)

%
%   EXAMPLE: [fnames, spikeIdx] = fnamesFromTxt(txtFile, spikeID)
%
% opens a .txt file and returns each line of the text file to a row of a
% matrix.  This is helpful when batch procesing data from a list of data
% files contained in a text file. Each row should have the same number of
% character elements. 'txtFile' should be the entire path to the txt file.
% 'spikeID' is a flag, which when set, returns the variable 'spikeIdx'.
% spikeIdx is a vector of unit numbers to consider for each file.
%
% This program ignores any text following a percent sign so that the user
% can include comments in the text file. file names and spike IDs should be
% tab delimited in the text file.
%
% CAH 06/08

% GDLH added 5/1/09
if (nargin == 0 || isempty(txtFile))
    currentdir = pwd;
    cd('N:\NexFiles\nexfilelists');
    [filename, pathname] = uigetfile('*.txt', 'Select a list of files');
    txtFile = [pathname,'\',filename];
    cd(currentdir);
end

if ~strcmpi('.txt', txtFile(end-3:end))
    error('Not a valid .txt file')
end

% determine string values for this computer. This requires PTB (i think)
SPACE = KbName('space');
TAB = KbName('tab');

fnames = [];
fid = fopen(txtFile);
if fid > 0;
    while(1)
        line = fgetl(fid);
        if ~ischar(line);
            break
        end
        
        %skip over the commented lines and carrage returns
        if (numel(line) > 0)% && ~strcmpi(line(1), '%');
            
            %strip out the comments
            cmt = findstr(line, '%');
            if ~isempty(cmt)
                line(cmt:end) = [];  
            end
            
            %strip off any trailing space
            lastChar = find((abs(line)~=TAB) & (abs(line)~=SPACE), 1, 'last');
            if ~isempty(lastChar)
                line(lastChar+1:end) = []; 
            end
            
            %paste into the matrix of text lines
            try
                if ~isempty(line)
                    fnames = [fnames; line];
                end
            catch
                disp('Cannot create filelist.');
                keyboard
            end
        end
    end
    fclose(fid);
    if exist('spikeID', 'var') && spikeID
       [r,c] = find(fnames == char(TAB), 1);
       spikeIdx = str2num(fnames(:,c+1:end));
       fnames(:,c:end) = [];
    end
else
    error(sprintf('Unable to open the specified file: \n %s \n', txtFile));
end

