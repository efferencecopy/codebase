function [GL_DATPATH, GL_DOCUPATH, GL_POPDATPATH, GL_ALLENPATH] = defineLocalPaths

%
% Defines the data paths to Glickfeld lab data files. These paths are
% declared as GLOBAL by startup.m. This function also changes directory to
% something more reasonable than the matlab startup dir.
%
% C.HASS 12/2013

switch whoami
    case 'hass_mbp'
        if exist('/Volumes/MobileCrash', 'dir') % external HD is plugged in
            GL_DATPATH = '/Volumes/MobileCrash/MobileCrash/Data/Mice/';
            GL_POPDATPATH = '/Volumes/MobileCrash/MobileCrash/Data/population_analysis_DBs/';
            GL_ALLENPATH = '/Volumes/MobileCrash/MobileCrash/Allen Images/';
        else % grab things from the local internal HD
            GL_DATPATH = '~/LabStuff/Data/Mice/';
            GL_POPDATPATH = ''; % not defined
            GL_ALLENPATH = ''; % not defined
        end
        
        GL_DOCUPATH = '~/github/DocuBase/';
        defaultpath = '~/LabStuff/';
        
    case 'hass_linux'
        GL_DATPATH = '~/Crash/Data/Mice/';
        GL_POPDATPATH = '~/Crash/Data/population_analysis_DBs/';
        GL_DOCUPATH = '~/Documents/CodeRepos/docubase/';
        GL_ALLENPATH = ''; % not defined
        defaultpath = '~/Crash/';
        
    case 'glick_rig1'
        GL_DOCUPATH = 'C:\Users\glickfeld_lab\Documents\docubase\';
        GL_POPDATPATH = '\\crash.dhe.duke.edu\charlie\Data\population_analysis_DBs\';
        GL_DATPATH =  '\\crash.dhe.duke.edu\charlie\Data\Mice\';
        GL_ALLENPATH = '\\crash.dhe.duke.edu\charlie\Allen Images\';
        defaultpath = 'C:\Users\glickfeld_lab\Documents\MATLAB';
        
    case 'nuke'
        GL_DOCUPATH = 'C:\Users\charlie\Documents\SourceTree_local\docubase\';
        GL_DATPATH =  '\\crash.dhe.duke.edu\charlie\Data\Mice\';
        GL_POPDATPATH = '\\crash.dhe.duke.edu\charlie\Data\population_analysis_DBs\';
        GL_ALLENPATH = '\\crash.dhe.duke.edu\charlie\Allen Images\';
        defaultpath = 'C:\Users\charlie\Documents\MATLAB';
end

cd(defaultpath)

% modest error checking
assert(~isempty(dir(GL_DATPATH)), '\n*********\n ERROR: %s is poorly defined. Sign into Crash, and re-run startup.m \n*********\n', 'GL_DATPATH')



