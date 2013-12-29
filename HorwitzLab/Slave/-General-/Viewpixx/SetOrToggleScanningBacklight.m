% If there are no input arguments, this function will toggle the scanning
% backlight mode. If you supply an argument (-1, 0, or 1), this function
% will toggle, turn off, or turn on the scanning backlight mode.
function scanningbl_on = SetOrToggleScanningBacklight(force_to)
if ~nargin || isempty(force_to)
    force_to = -1;
elseif nargin == 1
    force_to = min(max(force_to, -1), 1);
else
    scanningbl_on = -1;
    return
end

try
    PsychDataPixx('Open');
catch %#ok<CTCH>
    scanningbl_on = -1;
    return
end

status = Datapixx('GetVideoStatus');
switch force_to
    case -1
        if status.scanningBacklight
            PsychDataPixx('DisableVideoScanningBacklight');
        else
            PsychDataPixx('EnableVideoScanningBacklight');
        end
        scanningbl_on = ~status.scanningBacklight;
    case 0
        PsychDataPixx('DisableVideoScanningBacklight');
        scanningbl_on = 0;
    case 1
        PsychDataPixx('EnableVideoScanningBacklight');
        scanningbl_on = 1;
end
PsychDataPixx('Close');
