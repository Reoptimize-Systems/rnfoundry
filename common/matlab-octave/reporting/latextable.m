function latextable (data,varargin)
% Prints formatted matrix of numerical data with headings
% 
% Syntax
% 
% latextable (data, 'Parameter', Value)
% 
% Input
% 
% data - a matrix or cell array, containing the data to be put in the
%   table. If a matrix the numbers in the table will be printed using the
%   default format specifier (f), or the specifier(s) supplied in fms. data
%   can also be a cell array containing a mixture of strings and numbers.
%   each cell in this case can contain only a single scalar value. In this
%   case numbers are printed as in the matrix case, while strings are
%   printed up to the maximum width of the column.
% 
% Following this further options can be supplied as Parameter-Value pairs.
% The available options are:
%
% 'ColumnHeadings' 
%
%   a cell array of strings for the headings of each column. If no column
%   widths are supplied (see below), the columns will be made wide enough
%   to accomodate the headings.
% 
% 'NumberWidth' 
%
%   scalar or vector of column widths to use for the table. If scalar,
%   every column will have the same width. If a vector it must be of the
%   same length as the number of columns of data. If not supplied, and
%   column headers are supplied, the width of the widest column header will
%   be used. If not supplied and column headers are not supplied, a default
%   with of 16 characters is used.
% 
% 'FormatSpec'
%
%   a string, or cell array of strings containing format specifiers for
%   formatting the numerical output in each column. If a single string, the
%   same specifier is used for every column. If not supplied, the 'g'
%   specifier is used for every column.
% 
% 'RowHeadings'
%
%   a cell array of strings for the start of each row. If row headings are
%   supplied, the first column will be made wide enough to accomodate all
%   the headings.
% 
% 'FileID'
%
%   the file id to print to. Use 1 for stdout (to print to the command
%   line). Defaults to the command line if not supplied (i.e. same
%   behaviour as supplying the value 1).
% 
% 'BookTabs'
%
%   Flag indicating whether to include toprule midrule and bottomrule, as
%   provided by the latex booktabs package for nicely formatted tables.
%
% Example 
% 
% colheadings = {'number of projects','sales','profit'};
% rowheadings = {'Jimmy Slick', 'Norman Noob'}
% data = [3 rand(1) rand(1); 1 rand(1) rand(1)];
% 
% wid = 16;
% fms = {'d','.4f','.5E'};
% 
% fileID = 1;
% 
% latextable ( data, ...
%              'ColumnHeadings', colheadings, ...
%              'FormatSpec', fms, ...
%              'RowHeadings', rowheadings, ...
%              'BookTabs', true);
%
% See Also: DISPLAYTABLE, DISP, FPRINTF

% Created by Richard Crozier 2015

    options.NumberWidth = [];
    options.FormatSpec = 'g';
    options.FileID = 1;
    options.RowHeadings = {};
    options.ColumnHeadings = {};
    options.BookTabs = false;
    
    options = parse_pv_pairs (options, varargin);
        
	% use appropriate column separators and row endings to create a LaTeX
	% table
    colsep = ' & ';
    rowending = ' \\';
    
    if options.BookTabs
        fprintf (options.FileID, '\\toprule\n');
    end
    
    % loop through the column headings printing them out with the
    % desired column separator
    for i = 1:numel (options.ColumnHeadings)

        str = sprintf ('%s', options.ColumnHeadings{i});

        if i == numel (options.ColumnHeadings)
            % If we are at the end of a row, don't print the column
            % separator
            fprintf (options.FileID, '%s \\\\ \n', str);
        else
            % print the column header and column separator
            fprintf (options.FileID, ['%s',colsep], str);
        end

    end
    
    if options.BookTabs
        fprintf (options.FileID, '\\midrule\n');
    end
    
    displaytable ( data, ...
                   'ColHeadings', {}, ...
                   'ColWidth', options.NumberWidth, ...
                   'Format', options.FormatSpec, ...
                   'RowHeadings', options.RowHeadings, ...
                   'FileID', options.FileID, ...
                   'ColSep', colsep, ...
                   'RowEnding', rowending );

    if options.BookTabs
       fprintf (options.FileID, '\\bottomrule\n');
    end

end