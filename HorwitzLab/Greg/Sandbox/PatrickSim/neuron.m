function [response] = neuron(currentLMS,stimIntensity,quadparams,thresh)

% Isoresponse surface for K082609010
% threshold = 25.75 sp/sec
%quadparams = [536.1387 456.5097 27.2455 -21.9671 8.4686 -67.8732];
lms = currentLMS*stimIntensity;
variables = [lms(1)^2 lms(2)^2 lms(3)^2 2*lms(1).*lms(2) 2*lms(1).*lms(3) 2*lms(2).*lms(3)];
% Linearizing the contrast response function
response=poissrnd(sqrt(thresh.^2*(variables*quadparams')));
