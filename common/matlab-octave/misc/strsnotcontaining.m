function outstrs = strsnotcontaining(instrs, pattern)
% returns the subset of strings from a cell array of strings which do not
% contain the substring 'pattern'
%
% Syntax
%
% outstrs = strsnotcontaining(instrs, pattern)
%

% Copyright Richard Crozier 2012

    if ~iscellstr(instrs)
        error('instrs must be a cell array of strings');
    end
    
    if ~ischar(pattern)
        error('pattern must be a string');
    end

    nosubsstr = cellfun(@(x) isempty(strfind(x, pattern)), instrs, 'UniformOutput', true);

    outstrs = instrs(nosubsstr);
    
end