%%%%%%%%%%%%%%%%%%%%%%%%%
% Integrate this with GetVals?
function [out, p] = GetValsIfPossible(codes, p, type)

VALOFFSET = 4000;
out = [];

try
    if (length(codes) == 2 && codes(1) == codes(2))  % Bookending with double codes
        % GDLH 6/12/09 For backward compatibility assuming anything
        % bookended is double until otherwise specified.
        if (nargin < 3)
            type = 'double';
        end
        L = (p.events == codes(1));
        if (sum(L) < 2)
            return;
        else
            idxs = find(L);
            out = p.events(idxs(end-1)+1:idxs(end)-1)-VALOFFSET;
            % GDLH 3/12/12 commented out below line. We don't want to
            % eliminate negative values from the data stream - they might
            % be important.
            % out(out < 0) = [];
            % GDLH 6/12/09 commented out below line. uint2num is the more general function.
            % out = codes2num(out);     % assuming that all bookended codes are doubles
            out = uint2num(out,type);
            p.lastprocessed_t = p.times(max(idxs));
        end
    else   % Single codes or lists of codes
        if (nargin < 3)
            type = 'int';
        end
        [L, idxs] = ismember(codes, p.events);
        if (any(idxs == 0))
            return;
        end
        if (strcmpi(type,'int'))
            if (max(idxs)+1 <= length(p.events))
                out = p.events(idxs+1)-VALOFFSET;
                p.lastprocessed_t = p.times(max(idxs)+1);
            end
        elseif (strcmpi(type,'long'))
            if (max(idxs)+4 <= length(p.events))
                out = p.events(idxs+[1:4])-VALOFFSET;
                out = [2^0 2^8 2^16 2^24]*out;
                p.lastprocessed_t = p.times(max(idxs)+4);
            end
        elseif (strcmpi(type,'float'))
            if (max(idxs)+4 <= length(p.events))
                out = p.events(idxs+[1:4])-VALOFFSET;
                out = uint2num(out,'float');
                p.lastprocessed_t = p.times(max(idxs)+4);
            end
        elseif (strcmpi(type,'double'))
            if (max(idxs)+8 <= length(p.events))
                out = p.events(idxs+[1:8])-VALOFFSET;
                out = codes2num(out);
                p.lastprocessed_t = p.times(max(idxs)+8);
            end
        end
    end
catch ME
    keyboard
end
