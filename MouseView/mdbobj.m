classdef mdbobj
    % Mouse Database object
    
    properties
        mice = {};
        searchText = {}
    end
    
    methods
        function [match, idx] = search(obj, searchString)
            names = obj.mouseNames;
            [idx, match] = tsearch(obj.searchText, searchString, names);
            
            %return columns
            idx = idx(:);
            match = match(:);
        end
        
        function names = mouseNames(obj)
            names = cellfun(@(x) x.name, obj.mice, 'uniformoutput', false)';
        end
        
    end
    
end

