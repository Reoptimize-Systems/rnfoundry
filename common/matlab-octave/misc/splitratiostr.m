function [numerstr, denomstr] = splitratiostr(ratiofield, splitchar)
% attempts to split a series of strings representing ratios into two sets
% of strings representing the numerator and denominator of the ratio in the
% string
%
% Syntax
% 
% [numerstr, denomstr] = splitratiostr(ratiofield, splitchar)
%
% Description
%
% ratiofield is a cell array of strings contianing ratios, these are
% described as two strings separated by a second string denoting the split
% point. splitchar is a cell array of strings of the same size as
% ratiofield which contains the strings separating the numerator and
% denominator in each corresponding ratio in ratiofield
%
% Output
%
% numerstr and denomstr are cell arrays of strings of the same size. numerstr
% contains the numerator of each of the ratios and denomstr the denominator
%
% See also STRUCTRATIOS2STRUCTVALS

    k = strfind(ratiofield, splitchar);
    
    numerstr = cell(size(ratiofield));
    denomstr = cell(size(ratiofield));
    
    for i = 1:numel(k)
        
        if ~isempty(k{i})
            
            if k{i}(1) == 1
                
                error('STRUCTRATIOS2STRUCTVALS:splitratiofield:firstchars', ...
                    'Split character or string occurs in first characters of ratio name %s.', ratiofield{i})
                
            elseif k{i}(1) == length(ratiofield{i})
                
                error('STRUCTRATIOS2STRUCTVALS:splitratiofield:lastchars', ...
                    'Split character or string occurs in last characters of ratio name %s.', ratiofield{i})
                
            else
                
                numerstr{i} = ratiofield{i}(1:k{i}-1);
                
                denomstr{i} = ratiofield{i}(k{i} + length(splitchar):end);
                
            end
            
            
        else
            error('STRUCTRATIOS2STRUCTVALS:splitratiofield:nosplitchar', ...
                'Could not find split character in ratio field name %s.', ratiofield{i})
        end
        
    end

end