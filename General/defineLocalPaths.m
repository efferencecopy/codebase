function [GL_DATPATH, GL_DOCUPATH] = defineLocalPaths

%
% Defines the data paths to Glickfeld lab data files. These paths are
% declared as GLOBAL by startup.m
%

if ismac && strcmpi(license, '359028') %charlie's macbook pro
    GL_DATPATH = '~/LabStuff/Data/Mice/';
    GL_DOCUPATH = '/LabStuff/DocuBase/';
    defaultpath = '~/LabStuff/';
    
elseif isunix % charlie's desktop linux machine
    GL_DATPATH = '~/Crash/Data/Mice/';
    GL_DOCUPATH = '/DocuBase/';
    defaultpath = '~/Crash/';
    
elseif ispc
    GL_DATPATH =  '\\crash.dhe.duke.edu\charlie\Data\Mice\';
    GL_DOCUPATH = '\\crash.dhe.duke.edu\charlie\DocuBase';
    defaultpath = '\\crash.dhe.duke.edu\charlie\';
end

cd(defaultpath)