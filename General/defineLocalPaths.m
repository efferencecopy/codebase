function [GL_DATDIR, GL_DOCUPATH] = defineLocalPaths

%
% Defines the data paths to Glickfeld lab data files. These paths are
% declared as GLOBAL by startup.m
%

if ismac && strcmpi(license, '359028') %charlie's macbook pro
    GL_DATDIR = '~/Dropbox/Duke/Data/Mice/';
    GL_DOCUPATH = '/LabStuff/DocuBase/';
    defaultpath = '~/LabStuff/';
    
elseif isunix % charlie's desktop linux machine
    GL_DATDIR = '~/Crash/Data/Mice/';
    GL_DOCUPATH = '/DocuBase/';
    defaultpath = '~/Crash/';
    
elseif ispc
    GL_DATDIR =  '\\crash.dhe.duke.edu\charlie\Data\Mice\';
    GL_DOCUPATH = '\\crash.dhe.duke.edu\charlie\DocuBase';
    defaultpath = '\\crash.dhe.duke.edu\charlie\';
end

cd(defaultpath)