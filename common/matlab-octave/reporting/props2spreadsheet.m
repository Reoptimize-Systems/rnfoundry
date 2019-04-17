function props2spreadsheet (obj, filename, varargin)
% output numeric scalar and vector object properties to a spreadsheet
%
%

    options.WriteFcn = 'xlwrite';
    options.SheetName = 'props2spreadsheet Output';
    options.NameReplacements = {'', ''};
    options.FirstColumn = 1;
    
    options = parse_pv_pairs (options, varargin);
    
    writefcn = str2func (options.WriteFcn);

    propnames = properties (obj);
    
    for ind = 1:numel (propnames)
        
        varname = propnames{ind};
        
        repstrind = find (strcmp(options.NameReplacements(:,1), varname));
        
        if ~isempty (repstrind)
            varname = options.NameReplacements(repstrind, 2);
        end
        
        var = obj.(propnames{ind});
        
        if isvector (var) && isnumeric (var)
            
            writefcn (filename, {varname}, options.SheetName, xlsrange (ind+2, options.FirstColumn));
            writefcn (filename, var(:)', options.SheetName, xlsrange (ind+2, options.FirstColumn+1));
            
        end
        
    end

end