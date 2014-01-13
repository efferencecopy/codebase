function out = GetVal(codes, events, type)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function out = GetVal(codes, events, <type>)
%
% Generic function for pulling the value of the ecode 
% that follows another specified ecode.
%
% First input argument should be one or more 8000-level codes.
% 1) If a single code is given, we find that code in the event stream, take
%   the code following it, and return it after having subtracted VALOFFSET.
% 2) If a vector of codes are given, we do this same thing for each element
%   of the vector.
% 3) If a two element vector is given and the two elements are the same, 
%   the codes are treated as bookends.  Two instances of the code are
%   found and the values between them are extracted and returned.
%
% the optional argument, Type, can be "long", "float", or "double".
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    VALOFFSET = 4000;
    
    if (nargin < 3)
       type = 'int'; 
    end
    if (length(codes) == 2 & codes(1) == codes(2))  % Bookending with double codes
        L = (events == codes(1));
        if (sum(L) < 2)
            error(['Found only ',num2str(sum(L)),' ',num2str(codes(1)),' codes']);
        else
            idxs = find(L);
            out = events(idxs(1)+1:idxs(2)-1)-VALOFFSET;
            out(out < 0) = [];
        end
    else   % Single codes or lists of codes
        [L, idxs] = ismember(codes, flipud(events));  % ismember returns the *last* idx if more than 1
        if (all(idxs ~= length(events)))
            if (strcmpi(type,'int'))
                out = events(length(events)-idxs+2)-VALOFFSET;
                % Convoluted logic, but this unflips the flipud and takes the
                % next element the sequence
            elseif (strcmpi(type,'long'))
                out = events(length(events)-idxs+[2:5])-VALOFFSET;
                out = [2^0 2^8 2^16 2^24]*out;
            elseif (strcmpi(type,'double'))
                out = events(length(events)-idxs+[2:9])-VALOFFSET;
                out = codes2num(out);
            elseif (strcmpi(type,'float'))
                out = events(length(events)-idxs+[2:5])-VALOFFSET;
                out = uint2num(out,'float');
            end
        else
            out = nan;
        end
    end
end