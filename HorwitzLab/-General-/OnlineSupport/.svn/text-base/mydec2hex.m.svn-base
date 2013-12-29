function out = mydec2hex(in)

% function out = mydec2hex(in)
%
% Takes a signed decimal integer and converts it to hex.
% If called with a matrix argument it returns a row vector
% containing the hex representation of the individual elements
% in the usual Matlab order (down the first column, the down
% the second, etc.).

in = in(:);
out = dec2hex(abs(in));
newcol = zeros(size(in,1),2);
if (any(in >= 0))
    newcol(logical(in>=0),:) = repmat('  ',sum(in>=0),1);
end
if (any(in < 0))
    newcol(logical(in<0),:) = repmat(' -',sum(in<0),1);
end
out = [newcol, out];
out = reshape(out',1,numel(out));