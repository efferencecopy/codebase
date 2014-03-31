function out = pixPerMicron(nrows, ncols)
    
% val = pixPerMicron(xdim, ydim)
%
% Returns the pixperdeg for a specific camera. THIS CODE ASSUMES THAT THE
% OBJECTIVE WAS 2x for the Nikon or Retiga cameras. The code assumes a 5x
% objective for the slice rig. Currently there is no error checking for
% this... Future versions should be more flexible in this regard.
%
% C.Hass 03/2014

    % convert ot char so 'switch' is happy
    switch char([nrows, ncols])
        case char([640, 512])
            %  For nikon camera on the nikon camera 2x set to quick focus, 1000 um = 95 pix
            out = 95 ./ 1000;
            
        case char([1280, 1024])
            %  For nikon camera on the nikon camera 2x set to 1280x1024 fine or quick:
            %  1000 um = 191 pix.
            out = 191 ./ 1000;
            
        case char([3840, 1024])
            % For nikon camera on the nikon camera 2x. 1000 um = 570 pix
            out = 570 ./ 1000;
            
        case char([1040 1392])
            % for the retiga camera on the nikon scope 2x, 1000 um = 235 pix
            out = 235 / 1000;
            
        case char([768 1024])
            % for the slice rig 5x, 131.746 pix = 329.6195 um
            out = 131.746 ./ 329.6195;
            
        otherwise
            out = 191 ./ 1000;
            warning('unknown camera')
    end



end