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
    GL_DOCUPATH = 'C:\Users\glickfeld_lab\Desktop\CharlieDocs';
    defaultpath = 'C:\Users\glickfeld_lab\Documents\MATLAB';
    GL_DATPATH =  '\\crash.dhe.duke.edu\charlie\Data\Mice\';
end

cd(defaultpath)