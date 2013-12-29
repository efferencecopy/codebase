function [stimIntensityStair, oog] = stairCase2(currentLMS,numReversals,thresh,stepSize,scale,originalquadparams)

%% Variables
%currentLMS
stimIntensity=1;
stepSizeIndx=[];
reversal=0;
responseIndx=[];
stimIntensityStair=[];
respAtRev=[];
iteration=0;
iterationAtRev=[];

%% Test initial color directions at intensity=1
oog=isItOog(currentLMS,stimIntensity);
if oog==0
    response=neuron(currentLMS,stimIntensity,originalquadparams,thresh);
end


%% Staircase
while reversal<numReversals && oog==0  
    if response<=thresh
        norm=1+stepSize;
    else
        norm=1-stepSize;
    end
    iteration=iteration+1; %Keeps a running total of the total loop iterations
    stimIntensity=stimIntensity*norm;
    stimIntensityStair=[stimIntensityStair; stimIntensity]; %#ok<AGROW>  Keeps running total of cone stimulation
    oog=isItOog(currentLMS,stimIntensity);
    if oog==0
        response=neuron(currentLMS,stimIntensity,originalquadparams,thresh); %Response of cell to increased stim + noise
        responseIndx=[responseIndx; response]; %#ok<AGROW> Keeps running total of cell responses
    end
    if length(responseIndx)>1 && sign(thresh-responseIndx(end))~=sign(thresh-responseIndx(end-1)) && oog==0
        reversal=reversal+1; %Keeps track of reversals
        respAtRev=[respAtRev; responseIndx(iteration)]; %#ok<AGROW> Keeps running total of cell's response at reversal
        iterationAtRev=[iterationAtRev; iteration]; %#ok<AGROW>
        stepSizeIndx=[stepSizeIndx; stepSize]; %#ok<AGROW> Sanity check to see that step size is decreasing
        stepSize=stepSize*scale; %Reduces step size after each reversal
    end
end

%% Plot the staircase results
% if oog==0
%     figure(1); 
%     subplot(2,1,1) 
%     title('NeuroThresh Model')
%     hold on; grid on;
%     ylabel('Spike Rate')
%     plot(iterationAtRev, respAtRev, 'b o')
%     plot(responseIndx, 'b')
%     legend('Response at Reversal',4)
%     subplot(2,1,2); 
%     plot(stimIntensityStair, '--r')
%     hold on; grid on;
%     xlabel('Iteration of Stimulation')
%     ylabel('Stimulation Intensity')
% end

