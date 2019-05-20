function props2spreadsheet (obj, filename, varargin)
% output numeric scalar and vector object properties to a spreadsheet
%
% Syntax
%
% props2spreadsheet (obj, filename)
% props2spreadsheet (..., 'Parameter', Value)
%
% Description
%
% output numeric scalar and vector object properties to a spreadsheet
%
% Input
%
%  obj - matlab object with properties to be exported
%
%  filename - name of the output spreadsheet file
%
% Addtional arguments may be supplied as parameter-value pairs. The available options are:
%
%  'WriteFcn' - optional character vector containing the name of the 
%    function to be used to write the output spreadsheet. Default is
%    'xlwrite'.
%
%  'SheetName' - optional name of the sheet to which to export the data
%
%  'NameReplacements' - optional (n x 2) cell array containing character
%    vectors to be swapped to alternative names.
%
%  'FirstColumn' - optional scalar containing the number of the column in
%    which the property names (or their replacemnets) will be written
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