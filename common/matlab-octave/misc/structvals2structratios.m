function [s, fnames] = structvals2structratios(s, fields1, fields2, splitchar)
% creates ratios of two sets of fields in a structure in appropriately
% named fields in the same structure
%
% Syntax
% 
% s = structvals2structratios(s, fields1, fields2)
% s = structvals2structratios(s, fields1, fields2, splitchar)
% [s, fnames] = structvals2structratios(...)
%
% Description
%
% s is a structure containing the fields specified in fields1 and fields2
%
% s = structvals2structratios(s, fields1, fields2) adds fields to the
% structure 's' which contain the ratios of the values contained in the
% fields specified by the strings, or cell arrays of strings in fields1 and
% fields2. 
%
% By default the new ratio fields are given a name constructed by joining
% the two specified filed names with a capital 'V' character between them.
% An alternative character or string can be supplied with the optional
% argument 'splitchar'. splitchar can either be a single string, used for
% all supplied pairs of field names, or a cell array of strings of the same
% size as fields1 and fields2, in which case the corresponding value of
% splitchar is inserted between the field name pair.
%
% structvals2structratios divides the fields using the ./ operator and so
% will also work on fields containing matrices provided the normal size
% agreement rules are observed.
%
% The new structure names that were created in the process are returned in
% fnames, a cell string array.
%
% Example 1
%
% s.a = 1;
% s.b = 2;
% s = structvals2structratios(s, 'a', 'b')
%     
%     s = 
%     
%           a: 1
%           b: 2
%         aVb: 0.5000
%
% Example 2
%
% s.a = 1;
% s.b = 2;
% s.c = 2;
% s.d = 3;
% s = structvals2structratios(s, {'a', 'c'}, {'b', 'd'}, 'RATIO')
% 
%     s = 
% 
%               a: 1
%               b: 2
%               c: 2
%               d: 3
%         aRATIOb: 0.5000
%         cRATIOd: 0.6667
%
%

    % do some error checking on the input
    if ~(iscellstr(fields1) && iscellstr(fields2))
        if ischar(fields1) && ischar(fields2)
            fields1 = {fields1};
            fields2 = {fields2};
        else
            error('fields1 and fields2 must be strings or cell arrays of strings')
        end
    end
    
    if ~samesize(fields1, fields2)
        error('value field string arrays must be the same size.');
    end
    
    if nargin < 4
        splitchar = repmat({'V'}, size(fields1));
    else
        if ischar(splitchar)
            splitchar = repmat({splitchar}, size(fields1));
        elseif ~iscellstr(splitchar)
            error('ratiochar must be a string or cell array of strings')
        end
    end
    
    if ~every(isfield(s, [fields1, fields2]))
        error('Not all value fields are present in structure');
    end
    
    % now loop through the pairs of input fields and create new fields
    % based on their names and the supplied (or default) splitting string
    % with contents the division of the existing fields
    
    fnames = cell(size(fields1));
    
    for i = 1:numel(fields1)
    
        fnames{i} = [fields1{i}, splitchar{i}, fields2{i}];
        
        s.(fnames{i}) = s.(fields1{i}) ./ s.(fields2{i});

    end
    
end