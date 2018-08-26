function displaytable(data, colheadings, wid, fms, rowheadings, fid, colsep, rowending)
% Prints formatted matrix of numerical data with headings
% 
% Syntax
% 
% displaytable(data, colheadings, wid, fms, rowheadings, fid, colsep, rowending)
% 
% Input
% 
%  data - a matrix or cell array, containing the data to be put in the
%   table. If a matrix the numbers in the table will be printed using the
%   default format specifier (f), or the specifier(s) supplied in fms. data
%   can also be a cell array containing a mixture of strings and numbers.
%   each cell in this case can contain only a single scalar value. In this
%   case numbers are printed as in the matrix case, while strings are
%   printed up to the maximum width of the column.
% 
%  colheadings - a cell array of strings for the headings of each column.
%   Can be an empty cell array if no headings are required. If no column
%   widths are supplied (see below), the columns will be made wide enough
%   to accomdate the headings.
% 
%  wid - (optional) scalar or vector of column widths to use for the table.
%   If scalar, every column will have the same width. If a vector it must
%   be of the same length as the number of columns of data. If not
%   supplied, and column headers are supplied, the width of the widest
%   column header will be used. If not supplied and column headers are not
%   supplied, a default with of 16 characters is used.
% 
%  fms - (optional) a string, or cell array of strings containing format
%   specifiers for formatting the numerical output in each column. If a
%   single string, the same specifier is used for every column. If not
%   supplied, the 'g' specifier is used for every column.
% 
%  rowheadings - (optional) a cell array of strings for the start of each
%   row. Can be an empty cell array if no row headings are required. If row
%   headings are supplied, the first column will be made wide enough to
%   accomodate all the headings.
% 
%  fid - (optional) the file id to print to. Use 1 for stdout (to print to
%   the command line).
% 
%  colsep - (optional) A string or character to insert between every column.
%   The default separation string is ' | ', i.e. a space followed by a
%   vertical bar, followed by a space. A table suitible for inclusion in a
%   LaTeX document can be created using the ' & ' string, for example.
% 
%  rowending - (optional) An optional string or character to be appended at
%   the end of every row. Default is an empty string.
%
% 
% Examples
%
% Example 1 - Basic useage
% 
% colheadings = {'number of projects','sales','profit'};
% rowheadings = {'Jimmy Slick', 'Norman Noob'}
% data = [3 rand(1) rand(1); 1 rand(1) rand(1)];
% 
% To format the first number in each row as a decimal (%d), the second
% number %16.4f, and the third as %16.5E do the following:
% 
% wid = 16;
% fms = {'d','.4f','.5E'};
% 
% In this case 16 will be the field width, and '.5E' is what to use for the
% fms argument
% 
% fileID = 1;
% 
% >> displaytable(data,colheadings,wid,fms,rowheadings,fileID);
%             |number of projec |           sales |          profit 
% Jimmy Slick |               3 |          0.4502 |    5.22908E-001 
% Norman Noob |               1 |          0.9972 |    2.78606E-002 
%
% Example 2 - Produce a latex table
% 
% colheadings = {'number of projects','sales','profit'};
% rowheadings = {'Jimmy Slick', 'Norman Noob'};
% data = [3 rand(1) rand(1); 1 rand(1) rand(1)];
% wid = 16;
% fms = {'d'};
% 
% colsep = ' & ';
% rowending = ' \\';
% 
% fileID = 1;
% 
% >> displaytable(data,colheadings,wid,fms,rowheadings,fileID,colsep,rowending);
%
%             & number of projec &            sales &           profit \\
% Jimmy Slick &                3 &    6.948286e-001 &    3.170995e-001 \\
% Norman Noob &                1 &    9.502220e-001 &    3.444608e-002 \\
% 
% Example 3 - Mixed numeric and strings
%
% colheadings = {'number of projects','sales','profit'};
% rowheadings = {'Jimmy Slick', 'Norman Noob'};
% 
% data = {3, rand(1), rand(1); 
%         1, 'none', 0};
%     
% wid = 16;
% fms = {'d','.4f','.5E'};
% 
% fileID = 1;
% 
% >> displaytable(data,colheadings,wid,fms,rowheadings,fileID);
%             | number of projec |            sales |           profit
% Jimmy Slick |                3 |           0.4854 |     8.00280E-001
% Norman Noob |                1 |             none |     0.00000E+000
%
% See Also: DISP, FPRINTF

% Created by Richard Crozier 2012

    if nargin < 8 || isempty(rowending)
        rowending = '';
    end

    if nargin < 7 || isempty(colsep)
        colsep = ' | ';
    end

    % do some basic checking of input
    if nargin < 6 || isempty(fid)
        % print to the command line
        fid = 1;
    end
    
    if nargin < 5 || isempty(rowheadings)
        % no row headings supplied, use empty cell array
        rowheadings = {};
    end
    
    if nargin < 4 || isempty(fms)
        % no format specifiers supplied, use 'g' for all columns
        fms = 'g';
    end
    
    if nargin < 3 || isempty(wid)
        % default width is 10, this will be modified if column headers are
        % supplied
        wid = 10;
    end
    
    if nargin < 2
        colheadings = {};
    end

    if nargin >= 6 ...
            && ~isempty (rowheadings) ...
            && (~iscellstr(rowheadings) || ~isvector(rowheadings))
        error ('row headings must be vector cell array of strings');
    end
    
    tablestr = maketablestr (data, colheadings, wid, fms, rowheadings, colsep, rowending);
        
    % print to fid
    fprintf (fid, '%s', tablestr);
        
end