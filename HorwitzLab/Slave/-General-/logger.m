%LOGGER    Prints formatted data to the screen with timestamps.
%    LOGGER(FORMAT, ...) formats data and displays the results on the
%    screen.  Notice that the filename:linenumber entry is a hyperlink.
%
%    Example: Write a logging message to the screen from an mfile named
%    test.m on line 45 at 1:45:32 pm.
%
%        x = 12; y = 0.125;
%        logger('hello %s, x = %d, y = %0.3f\n', 'world', x, y);
%
%    Result: the formatted data prints to the screen
%
%    [13:45:32] <a href="matlab:">test.m:45</a>
%        hello world, x = 12, y = 0.125
%
%    See also FPRINTF.

% 2011/10/11 zlb

function logger(fmt, varargin)

if ~strcmp(fmt(end-1:end), '\n')
    fmt = [fmt '\n'];
end

relstr = version('-release');
if str2double(relstr(1:end-1)) > 2010
    useEditorAPI = true;
else
    useEditorAPI = false;
end

commandstr = '';
timestr = datestr(now,'HH:MM:SS');
stack = dbstack(1, '-completenames'); % omit 1 (this function's) frame

if ~isempty(stack)
    if useEditorAPI
        commandstr = sprintf( ...
            'matlab.desktop.editor.openAndGoToLine(''%s'',%d);', ...
            stack(1).file, stack(1).line);
    else
        commandstr = sprintf('opentoline(''%s'',%d,1);', ...
            stack(1).file, stack(1).line);
    end
else
    stack(1).name = 'cmdline'; stack(1).line = 1;
end

fprintf('[%s] <a href="matlab:%s">%s.m:%d</a>\n', ...
    timestr, commandstr, stack(1).name, stack(1).line);
fprintf(['\t' fmt], varargin{:});
