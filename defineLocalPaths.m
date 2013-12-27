
function [GL_PHYSPATH, GL_IMGPATH] = defineLocalPaths

%
% Defines the data paths to Glickfeld lab data files. These paths are
% declared as GLOBAL by startup.m
%

if ismac && strcmpi(license, '359028') %charlie's macbook pro
    GL_PHYSPATH = '~/Dropbox/Duke/Data/Physiology/';
    GL_IMGPATH = '~/Dropbox/Duke/Data/Imaging/';
    defaultpath = '~/LabStuff/';
    
elseif isunix % charlie's desktop linux machine
    GL_PHYSPATH = '~/Crash/Data/Mouse_Phys/';
    GL_IMGPATH = '~/Crash/Data/Imaging/';
    defaultpath = '~/Crash/';
    
elseif ispc
    GL_PHYSPATH  = '\\crash.dhe.duke.edu\charlie\Data\Mouse_Phys\';
    GL_IMGPATH =  '\\crash.dhe.duke.edu\charlie\Imaging\';
    defaultpath = '\\crash.dhe.duke.edu\charlie\';
end

cd(defaultpath)