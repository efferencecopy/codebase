function wtf

% wtf a.k.a What The @#&!
%
%  rethrows the last error and prints the contents of the error stack. This
%  function won't do anything special (and may crash) for errors committed
%  from the command line, but should be useful for debugging errors that
%  occur in functions (and subfunctions of functions) called from the
%  command window.
%
%CAH 04/08


e = lasterror;

fprintf('\n\n                   ***** ERROR LOG *****\n');
fprintf('  MESSAGE:\n');
disp(e.message);

if numel(e.stack) > 0
    fprintf('\n IN CALLING FUNCTION:\n');
    disp(e.stack(end,1));

    %roll through all the subfunctions (if there are any)
    nSubFun = length(e.stack) - 1;
    for a = 1:nSubFun;
        fprintf('\n  IN SUBFUNCTION:\n');
        disp(e.stack(a,1));
        fprintf('\n');
    end
end