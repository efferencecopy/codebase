function out = linmodfiterr(colordirections, responses, mechvect)

    dotprods = colordirections*mechvect';
    out = norm(abs(dotprods)-responses);
 %   keyboard
end
