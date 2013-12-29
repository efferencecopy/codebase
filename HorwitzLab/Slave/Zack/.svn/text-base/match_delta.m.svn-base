function match_delta()
global gl
KbName('UnifyKeyNames');
keys.quit = KbName('ESCAPE');
keys.inc_gun = KbName('UpArrow');
keys.dec_gun = KbName('DownArrow');
keys.step_finely = KbName('space'); % cycles through predefined gun steps
keys.togglegetimagelock = KbName('RightShift');
keys.getimage = KbName('Return');
keys.getimage = keys.getimage(1);

RestrictKeysForKbCheck(cell2mat(struct2cell(keys))');
ListenChar(2);

calfilename = 'ViewPixx.mat';
fundfilename = 'T_cones_smj.mat';
mondist = 100;
screenwidth = 57;
gl.windowPtr = 0;
gl.fliprequest = 0;

GUN = 2;

gl.bar.height = 1.5; % deg
gl.bar.width = 1.5;

gl.getnextframe = 0;
gl.saved_frames = [];
gl.saved_framenumber = [];
gl.allowframecapture = 1;

timetoleave = false;

gl.vpixx = ConnectedToViewPixx();

% gives about a 2% change in gun intensity per step
if gl.vpixx % 10-bit color depth
    init_step = 5 * 2^2 / (2^10-1);
else % Bits++: 14-bit color depth
    init_step = 5 * 2^6 / (2^14-1);
end
   
gl.GUNSTEP = linspace(0.001, init_step, 6);
gl.GUNSTEP = fliplr(gl.GUNSTEP);
  
InitDisplay(mondist, screenwidth, calfilename, fundfilename);

gl.peak_rgb = gl.bkgndrgb;
% gl.peak_rgb(GUN) = gl.bkgndrgb(GUN) + .2*rand;
gl.peak_rgb(GUN) = 0.918141;
gl.trough_rgb = gl.bkgndrgb;
gl.trough_rgb(GUN) = 0;

gun_names = {'red' 'green' 'blue'};

framecount = uint32(0);
while true
    [keyisdown,~,keycodes] = KbCheck();
    if keyisdown
        if keycodes(keys.quit)
            timetoleave = true;
        else
            InterpretKey(keys, keycodes, GUN);
        end
    end
    
    if ~mod(framecount, 2)
        Screen('FillRect', gl.windowPtr, gl.peak_rgb, gl.drawrect);
    else
        Screen('FillRect', gl.windowPtr, gl.trough_rgb, gl.drawrect);
    end
    
    DrawFormattedText(gl.windowPtr, sprintf('gun step: %0.2f%%\npeak %s: %g\ntrough %s: %g', ... 
        100*gl.GUNSTEP(end), gun_names{GUN}, gl.peak_rgb(GUN), gun_names{GUN}, ...
        gl.trough_rgb(GUN)), 'center', ceil(.75*gl.screenHeightpix), [0 0 0]);
    
    Screen('Flip', gl.windowPtr);
    framecount = framecount + 1;
    
    if gl.allowframecapture && gl.getnextframe
        img = Screen('GetImage', gl.windowPtr, gl.drawrect, [], 1);
        gl.saved_frames = [gl.saved_frames {img}];
        gl.saved_framenumber = [gl.saved_framenumber framecount];
        gl.getnextframe = gl.getnextframe - 1;
        if gl.getnextframe == 0
            gl.allowframecapture = 0;
        end
    end
    
    if timetoleave
        break
    end
end

sca();
ListenChar(0);
RestrictKeysForKbCheck([]);
fprintf('final peak rgb:   [%g %g %g]\n', gl.peak_rgb);
fprintf('final trough rgb: [%g %g %g]\n', gl.trough_rgb);
fprintf('background:       [%g %g %g]\n', gl.bkgndrgb);
fprintf('%s channel mean: %g\n', gun_names{GUN}, (gl.peak_rgb(GUN) + gl.trough_rgb(GUN))/2);

function InterpretKey(keys, keycodes, GUN)
global gl
if keycodes(keys.inc_gun)
    gl.trough_rgb(GUN) = min(gl.bkgndrgb(GUN), gl.trough_rgb(GUN) + gl.GUNSTEP(end));
elseif keycodes(keys.dec_gun)
    gl.trough_rgb(GUN) = max(0, gl.trough_rgb(GUN) - gl.GUNSTEP(end));
elseif keycodes(keys.step_finely)
    gl.GUNSTEP = [gl.GUNSTEP(end) gl.GUNSTEP(1:end-1)];
elseif keycodes(keys.togglegetimagelock)
    gl.allowframecapture = 1;
elseif gl.allowframecapture && gl.getnextframe == 0 && keycodes(keys.getimage)
    gl.getnextframe = 2; % save the next two frames
end

function InitDisplay(mondist, screenwidth, varargin) % varargin{1} and {2} are defined to be the calibration and fundamentals
global gl
calfilename = varargin{1};
load(calfilename);
cal = cals{end}; %#ok<USENS>
gl.mondistcm = mondist;
gl.screenWidthcm = screenwidth;
gl.bkgndrgb = cal.bgColor';
gl.cal.gammaTable = cal.gammaTable;
gl.cal.monSpd = cal.P_device;
gl.cal.invGamma = InvertGammaTable(cal.gammaInput, cal.gammaTable, 2^16);

if nargin > 3 && ~isempty(varargin{2})
    fundfilename = varargin{2};
    s = load(fundfilename);
    fns = fieldnames(s);
    gl.cal.fundamentals = s.(fns{1})';
    wavelength_spacing = s.(fns{2});
    P_device = SplineSpd(SToWls(cal.S_device), cal.P_device, SToWls(wavelength_spacing));
    gl.cal.M = gl.cal.fundamentals'*P_device;
    gl.cal.invM = inv(gl.cal.M);
end

% start up the imaging pipeline
if ~isempty(Screen('Windows'))
    gl.windowPtr = max(Screen('Windows'));
    Screen('FillRect', gl.windowPtr, cal.bgColor);
else
    PsychImaging('PrepareConfiguration');
    if gl.vpixx
        PsychImaging('AddTask', 'General', 'EnableDataPixxC48Output', 1); % in mode '1': every 2nd column of pixels is ignored
    else
        PsychImaging('AddTask', 'General', 'EnableBits++Color++Output', 1);
    end
    gl.windowPtr = PsychImaging('OpenWindow', 0, cal.bgColor);
end

gl.framerate = Screen('NominalFrameRate', gl.windowPtr, 1);
[screenwidthpix, screenheightpix]  = Screen('WindowSize', gl.windowPtr);
gl.screenWidthpix = screenwidthpix; % using Bits++ in Colour mode each pixel has a 1x2 aspect ratio
gl.screenHeightpix = screenheightpix;  % but Psychophysicstoolbox doesn't (need to) know about this
gl.screenCenterXpix = screenwidthpix/2;
gl.screenCenterYpix = screenheightpix/2;

pixpercm = gl.screenWidthpix/gl.screenWidthcm;
theta = atan2(gl.screenWidthcm/2, gl.mondistcm)*180/pi;
cmperdeg = gl.screenWidthcm/(2*theta);
gl.pixperdeg = pixpercm*cmperdeg;

gl.barheightinpix = round(gl.bar.height*gl.pixperdeg);
gl.barwidthinpix = round(gl.bar.width*gl.pixperdeg);

gl.drawrect = [gl.screenCenterXpix-floor(gl.barwidthinpix/2)...
    gl.screenCenterYpix-floor(gl.barheightinpix/2)...
    gl.screenCenterXpix+ceil(gl.barwidthinpix/2)...
    gl.screenCenterYpix+ceil(gl.barheightinpix/2)];

Screen('TextSize', gl.windowPtr, 20);

HideCursor();
