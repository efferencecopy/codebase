classdef blkobj
    % 
    
    
    properties
        sum = [];
        trial = [];
        ras = {};
        other = {};
    end
    
    methods
        function obj = blkobj(fileName, exptParamsDefs, mdb)
            %
            % find and import the data. mdb is an optional argument. If
            % it's not supplied, than the mdb will be generated de novo.
            %
            global GL_DATPATH
            
            if ~exist('fileName', 'var') || isempty(fileName) % no file name supplied
                currentDir = pwd;
                cd(GL_DATPATH)
                [fileName,fpath] = uigetfile({'*.nev'; '*.ns1'; '*.ns2'; '*.ns3'; '*.ns4'; '*.ns5'; '*.ns6'});
                fpath = [fpath,fileName];
                cd(currentDir);
                
            elseif any(regexpi(fileName, filesep))
                
                fpath = fileName; % the user supplied a fully qualified path
                
            else
                % narrow down the search for findfile.m
                if ~exist('mdb', 'var')
                    suppressVerbose = true;
                    mdb = initMouseDB('update', suppressVerbose);
                end
                mouseName = mdb.search(fileName(1:10)); %ignore the exact file name and just look at the date
                assert(numel(mouseName)==1, 'BLKOBJ ERROR: found too many mice with file name <%s>', fileName);
                
                % use find file to locate any potential blackrock files
                extensions = {'.nev', '.ns1', '.ns2', '.ns3', '.ns4', '.ns5', '.ns6'};
                for i_ext = 1:numel(extensions)
                    fpath{i_ext} = findfile(fileName, [GL_DATPATH, mouseName{1}], extensions{i_ext});
                end
                
                validPaths_idx = cellfun(@(x) ~isempty(x), fpath);
                assert(any(validPaths_idx), 'BLKOBJ ERROR: could not find a valid path')
            end
            
            % convert the blackrock files to a stro structure.
            if ~exist('exptParamsDefs', 'var')
                exdef = [];
            else
                exdef = exptParamsDefs; % execute the script of trial definitions
            end
            
            % call blk2stro to unpack the data
            stro = blk2stro('exdef', exdef,...
                           'nev', fpath{1},...
                           'ns1', fpath{2},...
                           'ns2', fpath{3},...
                           'ns3', fpath{4},...
                           'ns4', fpath{5},...
                           'ns5', fpath{6},...
                           'ns6', fpath{7});
            
            % package the output of the class 'blkobj'. This can't be
            % defined as the output of blk2stro b/c that is a structure and
            % not a blkobj class...
            obj.sum = stro.sum;
            obj.trial = stro.trial;
            obj.ras = stro.ras;
            obj.other = stro.other;
            
        end
        
    end %methods
    
end