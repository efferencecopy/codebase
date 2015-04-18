% game plan
%
%
%  use euler's method to determine the predicted response for a particular
%  set of depression and facilitation constants.
%
%  Technically, I only need to calculate the predictions once for each
%  unique spike train condition, and then calculate the residules for each
%  measured condition. 
%
%  I don't know if there's a faster way to fit a dynamical system. If there
%  is, that would be cool!
%
%  Varela 1997 used 2 depression factors and 1 facilitation factor. Each
%  factor consisted of a constant that was applied to each action
%  potential, and a recovery time constant that forced the factor back to
%  its resting state (a value of 1). 


% Step 1: import a file that's been pre-processed using the
% 'anlyMod_avgOuterleave' module. Make sure that it has the correct field.
% For example 'recovery' group is permissible. Bottom line, you're looking
% for a tdict structure.

% Step 2: iterate over channels. It would be nice if I can figure out a
% clever way of excluding junk channels.

% Step 2.1: For each trial type in the tdict, determine what the predicted
% pulse amplitude would be given a set of parameters. Do fminsearch to find
% the best fitting paramters.