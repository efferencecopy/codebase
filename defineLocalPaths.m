function [GL_PHYSPATH, GL_IMGPATH, GL_DOCUPATH] = defineLocalPaths

%
% Defines the data paths to Glickfeld lab data files. These paths are
% declared as GLOBAL by startup.m
%

if ismac && strcmpi(license, '359028') %charlie's macbook pro
    GL_PHYSPATH = '~/Dropbox/Duke/Data/Physiology/';
    GL_IMGPATH = '~/Dropbox/Duke/Data/Imaging/';
    GL_DOCUPATH = '/LabStuff/DocuBase/';
    defaultpath = '~/LabStuff/';
    
elseif isunix % charlie's desktop linux machine
    GL_PHYSPATH = '~/Crash/Physiology/';
    GL_IMGPATH = '~/Crash/Imaging/';
    GL_DOCUPATH = '/DocuBase/';
    defaultpath = '~/Crash/';
    
elseif ispc
    GL_PHYSPATH  = '\\crash.dhe.duke.edu\charlie\Physiology\';
    GL_IMGPATH =  '\\crash.dhe.duke.edu\charlie\Imaging\';
    defaultpath = '\\crash.dhe.duke.edu\charlie\';
end

cd(defaultpath)