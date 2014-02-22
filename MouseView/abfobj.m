classdef abfobj
    % creates a DT object. Basically the same as a STRO structure except
    % that a few other things are defined here. Indicies to trial events
    % are defined here too.
    
    properties
        head = [];
        dat = [];
        idx = [];
        name = [];
    end
    
    methods
        function obj = abfobj(fileName, mdb)
            %
            % find and import the data. mdb is an optional argument. If
            % it's not supplied, than the mdb will be generated de novo.
            %
            global GL_DATPATH
            
            if ~exist('fileName', 'var') || isempty(fileName) % no file name supplied
                currentDir = pwd;
                cd(GL_DATPATH)
                [fileName,fpath] = uigetfile('*.nex');
                fpath = [fpath,fileName];
                cd(currentDir);
            else
                % narrow down the search for findfile.m
                if ~exist('mdb', 'var')
                    mdb = initMouseDB();
                end
                [~, mouseName] = tsearch(mdb.searchText, fileName, {mdb.mice.name}');
                assert(~isempty(mouseName), 'ABFOBJ ERROR: could not find data file <%s>', fileName);
                assert(numel(mouseName)==1, 'ABFOBJ ERROR: too many data files match the input <%s>', fileName);
                
                % use find file to locate the .abf file
                fpath = findfile(fileName, [GL_DATPATH, mouseName{1}], '.abf');
                assert(~isempty(fpath), 'ABFOBJ ERROR: could not locate <%s> in directory <%s>', fileName, mouseName{1})
                
            end
            
            % convert the .abf file to a matlab structure. Define a few
            % other things.
            [obj.dat, obj.head] = my_abfload(fpath);
            obj.name = fileName;
            
            
        end
        
        function obj = removeSweeps(obj, idx)
            
            assert(size(obj.dat,3)>1, 'ABFOBJ ERROR: only one sweep present, nothing to delete')
            obj.dat(:,:,idx) = [];
        end
        
        function out = fitexp(raw)
            % do this as a cell fun type of thing?
            % Operate separately for cell inputs (cellfun) vs simple
            % inputs?
        end
        
        function out = grabData(obj, ch, timeStart, timeEnd)
        end
        
    end %methods
    
end