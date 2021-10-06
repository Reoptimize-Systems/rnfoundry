function s = structratios2structvals(s, ratiofields, base, splitchar)
% takes a structure containing fields representing ratios of two variables,
% and a based value field, and converts the ratios to sets of values which
% are added to the same structure.
% 
% Syntax
% 
% s = structratios2structvals(s, ratiofields, base, splitchar)
%
% Description
%
% structratios2structvals takes in a structure containing fields
% representing ratios, and a base value, and converts the ratios to values
% in appropriately named fields in the same structure. By default the
% ratios are supplied in the form n, where the 'V' character splits
% the quantities making up the ration. The conversion is done recursively
% following the chain of ratios.
%
% For example:
%
% s.bVa = 1;
% s.cVb = 2;
% s.dVb = 3;
% s.eVd = 10;
% s.a = 5;
% 
% ratiofields = {'bVa', 'cVb', 'dVb', 'eVd'};
% base = 'a';
% s = structratios2structvals(s, ratiofields, base)
%
% s = 
% 
%     bVa: 1.0000e+000
%     cVb: 2.0000e+000
%     dVb: 3.0000e+000
%     eVd: 10.0000e+000
%       a: 5.0000e+000
%       b: 5.0000e+000
%       c: 10.0000e+000
%       d: 15.0000e+000
%       e: 150.0000e+000
% 
% The splitting character can be any character or string, if the
% 'splitchar' optional input is supplied. e.g.
%
% s.bVERSUSa = 1;
% s.cVERSUSb = 2;
% s.dVERSUSb = 3;
% s.eVERSUSd = 10;
% s.a = 5;
% 
% ratiofields = {'bVERSUSa', 'cVERSUSb', 'dVERSUSb', 'eVERSUSd'};
% base = 'a';
% splitchar = 'VERSUS';
% s = structratios2structvals(s, ratiofields, base, splitchar)
% 
% s = 
% 
%       bVERSUSa: 1.0000e+000
%       cVERSUSb: 2.0000e+000
%       dVERSUSb: 3.0000e+000
%       eVERSUSd: 10.0000e+000
%              a: 5.0000e+000
%      bbVERSUSa: 1.0000e+000
%     ccVERSUSbb: 2.0000e+000
%     ddVERSUSbb: 3.0000e+000
%     eeVERSUSdd: 10.0000e+000
%              b: 5.0000e+000
%              c: 10.0000e+000
%              d: 15.0000e+000
%              e: 150.0000e+000
%
% Input
%
%  s - structure containing the fields to be evaluated
%
%  ratiofields - cell array of strings containing the names of the fields
%    to be converted from ratios to quantities
%
%  base - string containing field name containing base value for ratios
%
%  splitchar - optional string containing the character or string splitting
%    the numerator and denominator in the ratiofields
%
% Output
%
%  s - the input structure, but with the calculated output fields added
%

    % do some error checking on the input
    if ischar(ratiofields)
        ratiofields = {ratiofields};
    elseif ~iscellstr(ratiofields)
       error('STRUCTRATIOS2STRUCTVALS:badratiofields', ...
           '''ratiofields'' must be a string or cell array of strings.') 
    end

    if nargin < 4
        splitchar = repmat({'V'}, size(ratiofields));
    elseif iscellstr(splitchar) && numel(splitchar) > 1 && ~samesize(splitchar, ratiofields)
        error('STRUCTRATIOS2STRUCTVALS:badsplitchar', ...
            'If supplied splitchar must be a string, or single element cell string array, or a cell array of strings of the same size as ratiofields');
    elseif ischar(splitchar)
        splitchar = repmat({splitchar}, size(ratiofields));
    end
    
    if anyalldims(cellfun(@(x,y) length(x) < (length(y) + 2), ratiofields, splitchar))
       error('STRUCTRATIOS2STRUCTVALS:shortratio', ...
           'All ratio fields must contain at least two more characters than the corresponding split string') 
    end
    
    % check for duplicate ratios
    if numel(unique(ratiofields)) ~= numel(ratiofields)
        
        ratiofields = unique(ratiofields);
        
        warning('STRUCTRATIOS2STRUCTVALS:duplicates', ...
            '''ratiofields'' contained duplicates, using only unique fields.')
        
    end    
    
    % Now for each ratio field attempt to split it at the split character
    % to get the relevent numerator and denominator strings for the ratio
    numerstr = cell(size(ratiofields));
    denomstr = cell(size(ratiofields));
    for i = 1:numel(ratiofields)
        [numerstr(i), denomstr(i)] = splitratiostr(ratiofields(i), splitchar{i});
    end
    
    % check none of the denominators are the same as the numerators
    if anyalldims(strcmp(numerstr, denomstr))
        error('STRUCTRATIOS2STRUCTVALS:denomisnumerator', ...
            'Ratios cannot have the same numerator as denominator')
    end
    
    % recursively calculate the ratios in the structure, starting from
    % those with a give base name
    [s, numerstr] = recursivevalcalc(s, numerstr, denomstr, base, splitchar);
    
    if ~isempty(numerstr)
        warning('STRUCTRATIOS2STRUCTVALS:conversionfailure', ...
            'Not all ratios were converted to dimensions');
    end
    
end


function [s, newbases] = valcalc(s, numerstr, basestr, splitchar)
% valcalc calculates the value of the numerator of a ratio in a structure
% given the denominator, and puts the numerator in the structure as a new
% field
%
% Syntax
%
% [s, newbases] = valcalc(s, numerstr, basestr, splitchar)
%
% Input
%
% s - input structure containing the ratio fields constructed using
%   the numerstr, basestr and splitchar input (see below)
% 
% numerstr - A cell array of strings, where each string is the
%   numerator of a ratio field in the structure 's'.
% 
% basestr - A string containing the denominator of all the ratio
%   fields to be evaluated in 's'.
%
% splitchar - Cell array of strings or characters denoting the ratio
%   splitting string in the structure. 
%
% Description
%
% The ratio fields in the structure are constructed by joining each member
% of numerstr with the corresponding member of splitchar and the string
% basestr. The value fields are then made from the members of numerstr and
% are the constructed ratio fields multiplied by the base field
%
% See also STRUCTRATIOS2STRUCTVALS

    newbases = cell(size(numerstr));
    
    for i = 1:numel(numerstr)
        
        rationame = [numerstr{i},  splitchar{i}, basestr];
        
        s.(numerstr{i}) = s.(rationame) .* s.(basestr);
        
        newbases{i} = numerstr{i};
        
    end

end


function [s, numerstr, denomstr, splitchar] = recursivevalcalc(s, numerstr, denomstr, base, splitchar)
% recursively calculates the values in the structure based on a starting
% base value field and sets of field pairs
%
% Syntax
% 
% [s, numerstr, denomstr, splitchar] = recursivevalcalc(s, numerstr, denomstr, base, splitchar)
%
% Description
%
% recursivevalcalc recursively converts the ratio field pairs supplied in
% numerstr and denomstr to actual values starting from the base value in the
% structure field in 'base'. The first numerators created by the use of the
% starting base value are then used as new based values for any ratios not
% base on the original base value, and so on recursively until no ratio
% fields, or new base values remain to be evaluated.
%

    % Look for the base value in the ratios, we must start with ratios
    % denominated by this value
    basematch = strcmp(base, denomstr);

    % return if we can't find base in the list of denominators
    if ~anyalldims(basematch)
        return;
    end

    [s, newbases] = valcalc(s, numerstr(basematch), base, splitchar(basematch));

    % remove the completed ratios
    numerstr(basematch) = [];
    denomstr(basematch) = [];
    splitchar(basematch) = [];
    
    for i = 1:numel(newbases)
        % Recursion will end when their are no pairs left
        if ~isempty(numerstr)
            [s, numerstr, denomstr, splitchar] = recursivevalcalc(s, numerstr, denomstr, newbases{i}, splitchar);
        end
    end
        
end