% Test code working with the PR650
% This stuff only works in Windows Matlab
%%% Initialization stuff
EOM = char(10);
BYTESPERLINE = 17;
s = serial('COM2');
s.Inputbuffersize = 2000;
fopen(s);
fwrite(s,['s,,,,,,1,',EOM])
% Clearing the buffer
fread(s,s.BytesAvailable,'uchar');

fwrite(s,['m0' EOM]);
while (s.BytesAvailable == 0);
end

% Clearing the buffer
fread(s,s.BytesAvailable,'uchar');

fwrite(s,['d5', EOM]);
% Line below: 101 wavelengths + a header of the same length
[response] = fread(s,(101+1)*BYTESPERLINE,'uchar');
char(response')

response(1:17) = [];
data = [];
for i = 1:101
   tmp = response((i-1)*BYTESPERLINE+(2:4));
   wavelength = str2num(char(tmp'));
   tmp = response((i-1)*BYTESPERLINE+(7:15));
   power = str2num(char(tmp'));
   data = [data; wavelength, power];
end
plot(data(:,1),data(:,2),'k-');


%% Making the measurement
SerialComm('write', portNumber, ['m', EOM]);

tmpdata = '';
tic;
while (isempty(tmpdata) & toc < MAXWAITSEC)
    WaitSecs(0.05);
    tmpdata = char(SerialComm('readl', portNumber));  % Return from 'M' command
end
if (isempty(tmpdata))
    error('Measurement timed out.');
end

%% Retrieving the data
SerialComm('write', portNumber, ['d5', EOM]);
tmpdata = '';
tic;
while (isempty(tmpdata) & toc < MAXWAITSEC)
    tmpdata = char(SerialComm('readl', portNumber))  % return from 'D' command
    WaitSecs(0.05);
end
if (isempty(tmpdata))
    error('return timed out.');
else
    WaitSecs(5);  % Give the PR650 a chance to send the data.  It gets buffered.
    % A better way to do this of course is to wait until the next line
    % comes in.  2 seconds of wait here isn't enough.
    tmpdata = char(SerialComm('readl', portNumber))  % Header
end

data = [];
while (~isempty(tmpdata))
    tmpdata = char(SerialComm('readl', portNumber));
    commaPosition = find(tmpdata == ',');
    wavelength = str2num(tmpdata(1:commaPosition-1));
    power = str2num(tmpdata(commaPosition+1:end));
    data = [data; wavelength power];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Trying to use the psychtoolbox routines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

portNumber = 1;
PR650init(portNumber)
[spd, qual] = PR650measspd

%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
%% Code for finding a good background color
%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%

KEY.ESC = 41;
KEY.TILDE = 53;
KEY.ONE = 30;
KEY.TWO = 31;
KEY.THREE = 32;
KEY.Q = 20;
KEY.W = 26;
KEY.A = 4;
KEY.S = 22;
KEY.Z = 29;
KEY.X = 27;

rgb = [128 128 128];

[width, height]=Screen('DisplaySize', 0);
[window,screenRect] = Screen(0,'OpenWindow',rgb,[0 0 width height]);

oldClut = Screen('ReadNormalizedGammaTable', window);
Screen('LoadNormalizedGammaTable',window,linspace(0,1,256)'*ones(1,3));

Done = 0;
while (~Done)
   [keyisdown,secs,keycode] = KbCheck();
   if (keyisdown)
       if(find(keycode) == KEY.Q)
           rgb(1) = max([0 rgb(1)-1]);
       elseif (find(keycode) == KEY.W)
           rgb(1) = min([255 rgb(1)+1]);
       elseif (find(keycode) == KEY.A)
           rgb(2) = max([0 rgb(2)-1]);
       elseif (find(keycode) == KEY.S)
           rgb(2) = min([255 rgb(2)+1]);
       elseif (find(keycode) == KEY.Z)
           rgb(3) = max([0 rgb(3)-1]);
       elseif (find(keycode) == KEY.X)
           rgb(3) = min([255 rgb(3)+1]);
       elseif (find(keycode) == KEY.ESC)
           Done = 1;
       end
      Screen('FillRect',window,rgb);
      Screen('Flip', window);
   end
end
Screen('LoadNormalizedGammaTable', window, oldClut);
Screen('Close', window);
disp(rgb./[255 255 255]);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Looking for quantal steps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

portNumber = 1;
PR650init(portNumber);


[width, height]=Screen('DisplaySize', 0);
[window,screenRect] = Screen(0,'OpenWindow',[0 0 0],[0 0 width height]);

oldClut = Screen('ReadNormalizedGammaTable', window);
Screen('LoadNormalizedGammaTable',window,linspace(0,1,256)'*ones(1,3));

RGBs = [0 253 0 ; 0 254 0 ; 0 255 0];
spds = [];
nreps = 10;
for i = 1:nreps
    for j = 1:size(RGBs,1)
        Screen('FillRect', window, RGBs(j,:));
        Screen('Flip', window);
        [spd, qual] = PR650measspd;
        spds = [spds, spd];
    end
end

Screen('LoadNormalizedGammaTable', window, oldClut);
Screen('Close', window);

% And now, a little analysis...
[u,s,v] = svd(spds);
figure;
plot(v(:,1),'k.-');

mat = reshape(v(:,1),size(RGBs,1),nreps);
figure;
hist(mat');

% Coefficient of variation is probably a good measure of spread
std(v(:,1))/abs(mean(v(:,1)))

%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%
%% Code for testing to see whether the Bits++ box is actually increasing the
% color depth.
%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%

portNumber = 1;
PR650init(portNumber);
spds = [];
for i = 11:30
   disp(['Current level: ',num2str(mod(i-1,3)+1)]);
   input('Hit return to make a measurement');
   [spd, qual] = PR650measspd;
   spds = [spds, spd];
end

% The upshot from the above two cells of code is that the Mac Pro using the
% Ati Radeon card stepping the green gun from 253, 254, 255 leads to PR650
% measurements that are resolvably different.  Doing the same procedure
% using Charlie's G4 (unknown graphics card) and the Bits++ box, but the
% same monitor stepping the green gun from 65021, 65022, 65023 (the lowest
% value, 65021, should be driving the green gun at 99.22% of max or
% 253/255, but we forgot to check this) did not lead to resolvable values
% (all the values were very close together).  Seems like the Bits++ box is
% working correctly.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Testing whether LoadNormalizedGammaTable is
%% working the way I think it should.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

portNumber = 1;
PR650init(portNumber);

[window,screenRect] = Screen('OpenWindow',1,[0 0 0]);

Screen('LoadNormalizedGammaTable',window,linspace(0,1,256)'*ones(1,3));
nvals = 15;
vals = linspace(0,255,nvals);

spds = [];
for i = 1:nvals
    Screen('FillRect', window, [0 vals(i) 0]);
    Screen('Flip', window);
    [spd, qual] = PR650measspd;
    spds(:,i) = spd;
end

[u,s,v] = svd(spds);

if(max(abs(u(:,1))) ~= max(u(:,1)))
   v = -v;
end
gamma = v(:,1);
% ignoring blacklight for now
% Scaling so we're between [0, 1]
gamma = (gamma-min(gamma))./(max(gamma)-min(gamma))

invgamma = interp1(gamma, vals/255, linspace(0,1,256));

% OK, now trying it again to see if LoadNormalizedGammaTable
% is doing the trick.

Screen('LoadNormalizedGammaTable',window,repmat(invgamma',1,3));

newspds = [];
for i = 1:nvals
    i
    Screen('FillRect', window, [0 vals(i) 0]);
    Screen('Flip', window);
    [spd, qual] = PR650measspd;
    newspds(:,i) = spd;
end

Screen('Close',window);

[u,s,v] = svd(newspds);
plot(v(:,1),'k.')
plot(diff(v(:,1)),'k.')

% Cool!  This all works with the Bits++ in Colour mode.  The reason is
% probably that when we measure the gamma curves we make linear increases 
% in both the high and low order bytes, with the result that we're linearly
% increasing the voltage to the electron gun.  So we're going from [0, 0]
% to [255, 255] in even increments.  So I can probably use all the standard
% calibration stuff with Bits++ in colour mode but when I draw things to
% the screen I should set the gamma functions on the graphics card to flat
% and deal with the gamma stuff in software - the hardware gamma table
% treats the high and low order bytes identically which is not what we
% want.

%%
% Checking phosphor spectra with different profiles, different heads,
% different monitor configurations, etc.

portNumber = 1;
PR650init(portNumber);

[width, height]=Screen('DisplaySize', 0);
[window,screenRect] = Screen('OpenWindow',0,[0 0 0],[0 0 width height]);

Screen('LoadNormalizedGammaTable',window,linspace(0,1,256)'*ones(1,3));
rgbs = [190 0 0; 0 190 0; 0 0 190];
nreps = 5;
spds = zeros(81,3,nreps);
for i = 1:nreps
    for j = 1:size(rgbs,1)
        Screen('FillRect', window, rgbs(j,:));
        Screen('Flip', window);
        [spds(:,j,i), qual] = PR650measspd;
    end
end

plot(reshape(spds,[81,3*nreps]),'k-');

Screen('CloseAll');
%%
% Just changing the normalized lookup table in the two monitor
% configuration to see if it affects both monitors.
% The upshot of this experiment is that using both of the graphics card
% heads is a bad idea.  Use one with a splitter - it's safer.  The reason
% is that "LoadNormalizedGammaTable" only knows about one of the displays
% and it's not always obvious which one this is!

[width, height]=Screen('DisplaySize', 0);
[window,screenRect] = Screen('OpenWindow',0,[0 0 0],[0 0 width height]);
for i = 1:10 
    Screen('LoadNormalizedGammaTable',window,unifrnd(0,1,256,3));
    pause(.25);
end
Screen('LoadNormalizedGammaTable',window,linspace(0,1,256)'*ones(1,3));
Screen('CloseAll');
%%
% Changing 
%