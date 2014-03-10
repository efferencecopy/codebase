function [GL_DATPATH, GL_DOCUPATH] = defineLocalPaths

%
% Defines the data paths to Glickfeld lab data files. These paths are
% declared as GLOBAL by startup.m. This function also changes directory to
% something more reasonable than the matlab startup dir.
%
% C.HASS 12/2013

switch whoami
    case 'hass_mbp'
        if exist('/Volumes/Charlie HD', 'dir') % external HD is plugged in
            GL_DATPATH = '/Volumes/Charlie HD/Crash HD/Mice/';
            GL_DOCUPATH = '/Volumes/Charlie HD/Rig DocuBase/';
        else % grab things from the local internal HD
            GL_DATPATH = '~/LabStuff/Data/Mice/';
            GL_DOCUPATH = '~/LabStuff/DocuBase/';
        end
        
        defaultpath = '~/LabStuff/';
        
    case 'hass_linux'
        GL_DATPATH = '~/Crash/Data/Mice/';
        GL_DOCUPATH = '~/DocuBase/';
        defaultpath = '~/Crash/';
        
    case 'glick_rig1'
        GL_DOCUPATH = 'C:\Users\glickfeld_lab\Documents\docubase\';
        GL_DATPATH =  '\\crash.dhe.duke.edu\charlie\Data\Mice\';
        defaultpath = 'C:\Users\glickfeld_lab\Documents\MATLAB';
        
    case 'nuke'
        GL_DOCUPATH = 'C:\Users\charlie\Documents\GitHub\docubase\';
        GL_DATPATH =  '\\crash.dhe.duke.edu\charlie\Data\Mice\';
        defaultpath = 'C:\Users\charlie\Documents\MATLAB';
end

cd(defaultpath)