function spds = DarkAdaptTest(device_name)
% Numpad 7 and 8: red-   & red+
% Numpad 4 and 5: green- & green+
% Numpad 1 and 2: blue-  & blue+
% Numpad ENTER cycles through color set for spatial measurements
% Escape quits the script

switch upper(device_name)
    case 'PR705'
        use_device = 2;
    case 'PR650'
        use_device = 1;
    case 'NONE'
        use_device = 0;
    otherwise
        error('Unrecognized device name');
end

time0 = GetSecs(); % load the mex file into memory

hport = -1;
device_handler('open');

use_viewpixx = PromptScanningBacklightMode(0);

if use_viewpixx
    PsychImaging('PrepareConfiguration');
    PsychImaging('AddTask', 'General', 'EnableDataPixxC48Output', 0);
    [window,screenRect] = PsychImaging('OpenWindow', 0, [0 0 0]);
else
    [window,screenRect] = Screen(0, 'OpenWindow', [0 0 0]);
    Screen('LoadNormalizedGammaTable', window, ...
        linspace(0, 1, 256)' * ones(1, 3));
end
Screen('TextSize', window, 20);

spds = {};

screenwidth = 57; % cm
mondist = 100; % cm
box_width_deg = 2.5; % deg
% pixpercm = sqrt(screenRect(3)^2 + screenRect(4)^2) / 23.6 / 2.54;
% 23.6" diagonal (advertised), 2.54 cm per inch.

theta = atand(screenwidth/2/mondist);
pixperdeg = screenRect(3)/2/theta;

box_template = box_width_deg * pixperdeg * [-1 -2 1 2] / 2;

box_ = [CenterRectOnPoint(box_template, 3.75*pixperdeg, 2*3.25*pixperdeg);
    CenterRectOnPoint(box_template, 3.75*pixperdeg, 2*8.75*pixperdeg);
    CenterRectOnPoint(box_template, 3.75*pixperdeg, 2*14.25*pixperdeg);
    CenterRectOnPoint(box_template, 15.75*pixperdeg, 2*3.5*pixperdeg);
    CenterRectOnPoint(box_template, 16*pixperdeg, 2*8.75*pixperdeg);
    CenterRectOnPoint(box_template, 15.75*pixperdeg, 2*14.25*pixperdeg);
    CenterRectOnPoint(box_template, 28.5*pixperdeg, 2*3.75*pixperdeg);
    CenterRectOnPoint(box_template, 28.5*pixperdeg, 2*8.75*pixperdeg);
    CenterRectOnPoint(box_template, 28.25*pixperdeg, 2*14.25*pixperdeg)]';

spatial_meas = size(box_, 2) >= 9; % are we doing the 9 spatial homogeneity measurements?

color_queue = .15 * [1 1 1; 1 0 0; 0 1 0; 0 0 1; 0 0 0];
box_rgb = color_queue(1,:);

do_refresh = true;
leave = false;

ListenChar(2);

while true
    [is_keydown,null,keycodes] = KbCheck();
    if is_keydown
        parse_key(find(keycodes, 1));
    end
    if do_refresh
        refresh_box();
        draw_text();
        do_refresh = false;
        Screen('Flip', window);
    end
    
    if leave, break; end
end

device_handler('save');
device_handler('close');
ListenChar(0);

if use_viewpixx
    CloseViewpixx();
end
sca();

    function parse_key(keycode)
        KbReleaseWait();
        if keycode == 41 % escape
            leave = true;
        elseif keycode == 40 % get spectra
            pause(.1);
            try
                time0 = GetSecs();
                spd = device_handler('measure');
                spds = [spds; {GetSecs()-time0 box_rgb spd}];
            catch Ex
                device_handler('close')
                ListenChar(0);
                if use_viewpixx
                    CloseViewpixx();
                end
                sca();
                rethrow(Ex);
            end
        else
            update_box_rgb(keycode);
        end
    end

    function update_box_rgb(keycode)
        GUNSTEP = 5;
        do_refresh = true;
        switch keycode
            case 95 % red-
                box_rgb(1) = max([0 box_rgb(1) - GUNSTEP]);
            case 92 % green-
                box_rgb(2) = max([0 box_rgb(2) - GUNSTEP]);
            case 89 % blue-
                box_rgb(3) = max([0 box_rgb(3) - GUNSTEP]);
            case 96 % red+
                box_rgb(1) = min([255 box_rgb(1) + GUNSTEP]);
            case 93 % green+
                box_rgb(2) = min([255 box_rgb(2) + GUNSTEP]);
            case 90 % blue+
                box_rgb(3) = min([255 box_rgb(3) + GUNSTEP]);
            case 88 % numpad enter: cycle through colors
                color_queue = color_queue([2:end 1], :);
                box_rgb = color_queue(1,:);
            otherwise
                do_refresh = false;
        end
    end

    function spd = device_handler(request)
        spd = [];
        if use_device
            switch request
                case 'open'
                    if use_device == 2
                        hport = OpenPortUsingIO(1000, 5);
                        % 1 sec exposure, 5 measurements
                    else
                        CMCheckInit();
                    end
                case 'close'
                    if use_device == 2
                        ClosePortUsingIO(hport);
                    else
                        PR650close();
                    end
                case 'measure'
                    if use_device == 2
                        [null, spd] = ...
                            ParseSpectrum705(MeasureSpectraUsingIO(hport)); %#ok<SETNU>
                    else
                        spd = PR650measspd([380 4 101]);
                    end
                case 'save'
                    save(sprintf('%s-DarkAdaptTest.mat', datestr(now, ...
                        'yyyymmddTHHMMSS')), 'spds');
            end
        end
    end

    function refresh_box()
        Screen('FillRect', window, box_rgb, box_);
    end

    function draw_text()
        if ~spatial_meas
            DrawFormattedText(window, sprintf('[%d %d %d]', box_rgb), ...
                screenRect(3) * 3 / 4, screenRect(4) / 4, [255 255 255]);
        end
    end
end % end main function