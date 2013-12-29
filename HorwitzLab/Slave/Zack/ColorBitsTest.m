% assumes you're hooked up to the PR-705

boxSize = 200;
Xoffset = 0; Yoffset = 0;
screenid = 0;

PsychImaging('PrepareConfiguration');
PsychImaging('AddTask', 'General', 'EnableBits++Color++Output', 1);
[wPtr,screenRect] = PsychImaging('OpenWindow', screenid, .5);

CMCheckInit(6);

boxRect = [0 0 boxSize boxSize];
boxRect = CenterRect(boxRect, screenRect);
boxRect = boxRect + [Xoffset -Yoffset Xoffset -Yoffset];

boxMat = zeros(boxSize, boxSize, 3);
boxMat(:,:,2) = 1;
counter = 1;
whichlevels = (64000+(0:1:20))/2^16;
spdMat = [];
S = [380 2 201];
input('');
try
    for mult = whichlevels
        tex = Screen('MakeTexture', wPtr, boxMat * mult, [], [], 2);
        Screen('DrawTexture', wPtr, tex, [], boxRect, [], 0);
        Screen('Flip', wPtr);
        spd = PR705measspd(S);
        if isempty(spd), break; end
        spdMat(counter,:) = spd(:);
        counter = counter + 1;
    end
    Screen('Close', tex);
catch ME
    CMClose(6);
    sca();
    disp(getReport(ME));
    rethrow(ME);
end
save tmpdata
pause(2);
CMClose(6);
sca();

[u,s,v] = svd(spdMat');
figure(); plot(v(:,1),'ko'); title('new method');
