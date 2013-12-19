
function [PHYSPATH, IMGPATH] = defineLocalPaths

%
% Defines the data paths to Glickfeld lab data files. These paths are
% declared as GLOBAL by startup.m
%

if ismac && strcmpi(license, '359028') %charlie's macbook pro
    PHYSPATH = '~/Dropbox/Duke/Data/Physiology/';
    IMGPATH = '~/Dropbox/Duke/Data/Imaging/';
    defaultpath = '~/LabStuff/';
    
elseif isunix % charlie's desktop linux machine
    PHYSPATH = '~/Crash/Data/Mouse_Phys/';
    IMGPATH = '~/Crash/Data/Imaging/';
    defaultpath = '~/Crash/';
    
elseif ispc
    PHYSPATH  = '\\crash.dhe.duke.edu\charlie\Data\Mouse_Phys\';
    IMGPATH =  '\\crash.dhe.duke.edu\charlie\Imaging\';
    defaultpath = '\\crash.dhe.duke.edu\charlie\';
end

cd(defaultpath)