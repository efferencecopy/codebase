function out = pixPerMicron(xdim, ydim)
    

    % convert ot char so 'switch' is happy
    switch char([xdim, ydim])
        case char([640, 512])
            % For the nikon camera set to quick focus, 1000 um = 95 pix
            out = 95 ./ 1000;
            
        case char([1280, 1024])
            % For the nikon camera set to 1280x1024 fine or quick:
            %  1000 um = 191 pix.
            out = 191 ./ 1000;
            
        case char([3840, 1024])
            % For nikon camera. 1000 um = 570 pix
            out = 570 ./ 1000;
            
        case char([1040 1392])
            error('Regita calibration unknown');            
    end



end