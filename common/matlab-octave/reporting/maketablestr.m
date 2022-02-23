function tablestr = maketablestr (data, varargin)
% makes a char vector with formatted matrix of numerical data with headings
% 
% Syntax
% 
% tablestr = maketablestr (data)
% tablestr = maketablestr (..., 'Parameter', value)
% 
% Input
% 
%  data - a matrix or cell array, containing the data to be put in the
%   table. If a matrix the numbers in the table will be printed using the
%   default format specifier (f), or the specifier(s) supplied in 'Format'
%   option. data can also be a cell array containing a mixture of strings
%   and numbers. each cell in this case can contain only a single scalar
%   value. In this case numbers are printed as in the matrix case, while
%   strings are printed up to the maximum width of the column.
% 
% Additional optional arguments can be supplied using parameter-value
% pairs. The available options are:
%
%  'ColHeadings - a cell array of strings for the headings of each column.
%    Can be an empty cell array if no headings are required. If no column
%    widths are supplied (see below), the columns will be made wide enough
%    to accomdate the headings.
% 
%  'ColWidth - (optional) scalar or vector of column Widths to use for the
%    table. If scalar, every column will have the same width. If a vector
%    it must be of the same length as the number of columns of data. If not
%    supplied, and column headers are supplied, the width of the column
%    header will be used. If not supplied and column headers are not
%    supplied, a default with of 16 characters is used.
% 
%  'Format - (optional) a string, or cell array of strings containing format
%    specifiers for formatting the numerical output in each column. If a
%    single string, the same specifier is used for every column. If not
%    supplied, the 'g' specifier is used for every column.
% 
%  'RowHeadings - (optional) a cell array of strings for the start of each
%    row. Can be an empty cell array if no row headings are required. If
%    row headings are supplied, the first column will be made wide enough
%    to accomodate all the headings.
% 
%  'ColSep - (optional) A string or character to insert between every 
%    column. The default separation string is ' | ', i.e. a space followed
%    by a vertical bar, followed by a space. A table suitible for inclusion
%    in a LaTeX document can be created using the ' & ' string, for
%    example.
% 
%  'RowEnding' - (optional) An optional string or character to be appended 
%    at the end of every row. Default is an empty string.
%
% Outut
%
%  tablestr - character vector containing the generated text table
%
% Examples
%
% Example 1 - Basic useage
% 
% colheadings = {'number of projects','sales','profit'};
% rowheadings = {'Jimmy Slick', 'Norman Noob'}
% data = [3 rand(1) rand(1); 1 rand(1) rand(1)];
% 
% % To format the first number in each row as a decimal (%d), the second
% % number %16.4f, and the third as %16.5E do the following:
% 
% wid = 16;
% fmt = {'d','.4f','.5E'};
% 
% % In this case 16 will be the field width, and '.5E' is what to use for the
% % options.Format argument
% 
% >> maketablestr(data, 'ColHeadings', colheadings, 'ColWidth', wid, 'Format', fmt, 'RowHeadings', rowheadings);
%             |number of projec |           sales |          profit 
% Jimmy Slick |               3 |          0.4502 |    5.22908E-001 
% Norman Noob |               1 |          0.9972 |    2.78606E-002 
%
% Example 2 - Produce a latex table
% 
% colheadings = {'number of projects','sales','profit'};
% rowheadings = {'Jimmy Slick', 'Norman Noob'}
% data = [3 rand(1) rand(1); 1 rand(1) rand(1)];
% wid = 16;
% options.Format = {'d'};
% 
% colsep = ' & ';
% rowending = ' \\';
%
% >> maketablestr(data, 'ColHeadings', colheadings, 'ColWidth', wid, 'Format', fmt, 'RowHeadings', rowheadings, 'ColSep', colsep, 'RowEnding', rowending);
%
%             & number of projec &            sales &           profit \\
% Jimmy Slick &                3 &    6.948286e-001 &    3.170995e-001 \\
% Norman Noob &                1 &    9.502220e-001 &    3.444608e-002 \\
% 
% Example 3 - Mixed numeric and strings
%
% options.ColHeadings = {'number of projects','sales','profit'};
% options.RowHeadings = {'Jimmy Slick', 'Norman Noob'};
% 
% data = {3, rand(1), rand(1); 
%         1, 'none', 0};
%     
% options.ColWidth = 16;
% options.Format = {'d','.4f','.5E'};
% 
% fileID = 1;
% 
% >> displaytable(data,options.ColHeadings,options.ColWidth,options.Format,options.RowHeadings,fileID);
%             | number of projec |            sales |           profit
% Jimmy Slick |                3 |           0.4854 |     8.00280E-001
% Norman Noob |                1 |             none |     0.00000E+000
%
% See Also: DISP, FPRINTF
%

% Created by Richard Crozier 2018

    options.RowEnding = '';
    options.RowStart = '';
    options.ColSep = ' | ';
    options.RowHeadings = {};
    options.Format = 'g';
    options.ColWidth = [];
    options.ColHeadings = {};
    options.WrapColHeadings = false;
    options.ColHeaderRule = false;
    options.RuleChar = '-';
    options.TopRule = false;
    options.BottomRule = false;
    
    options = parse_pv_pairs (options, varargin);

    tablestr = '';
    
    if ~isempty (options.RowHeadings) ...
            && (~iscellstr(options.RowHeadings) || ~isvector(options.RowHeadings))
        error ('row headings must be vector cell array of strings');
    end
    
    % get the numbers of rows and columns of data
    [nRowD,nColD] = size(data);
    
    % check that sensible format specifiers have been supplied
    if ~iscellstr(options.Format)
        
        if ischar(options.Format)
            options.Format = repmat({options.Format},1,nColD);
        else
            error('options.Format must be a string or cell array of strings.');
        end
        
    elseif isempty(options.Format)
        
        options.Format = repmat({'f'},1,nColD);
        
    elseif numel(options.Format) ~= nColD
        
        if numel(options.Format) == 1
            options.Format = repmat(options.Format,1,nColD);
        else
           error('options.Format does not have the same number of format specifiers as there are columns.');
        end
        
    end
    
    % replace empty format specifiers with 'f'
    for formatstr_ind = 1:numel(options.Format)
        if isempty(options.Format{formatstr_ind})
            options.Format{formatstr_ind} = 'f';
        end
    end
    
    [nRowFms,nColFms] = size(options.Format);
    if(nRowFms>1)
        error ('options.Format can not have more than one row');
    end
    
    rhwid = 0;
    if ~isempty(options.RowHeadings)
        
        if ~iscellstr(options.RowHeadings)
            error('options.RowHeadings must be a cell array of strings');
        end
        
        if numel(options.RowHeadings) ~= size(data, 1)
            error('Rowheadings must be a cell array of strings of the same size as the number of rows in data.')
        end
    
        for r = 1:numel(options.RowHeadings)
            % find the maximum width the first column must be to accomodate
            % the row headings
            rhwid = max(rhwid, length(options.RowHeadings{r}));
        end
        
    end
    
    if isempty (options.ColWidth)

        % no column width supplied, so determine a sensible value for the
        % width of each column

        if isempty (options.ColHeadings)
            
            options.ColWidth = 10;
            
        else
            
            tempoptions.ColWidth = repmat(10, size(options.ColHeadings));
            for colheadind = 1:numel(options.ColHeadings)

                % get a column width which is the minimum to accept the
                % column heading length or the default width specification
                tempoptions.ColWidth(colheadind) = max(length(options.ColHeadings{colheadind}), tempoptions.ColWidth(colheadind));

                if tempoptions.ColWidth < 1
                    error('Column width is less than 1, and the column header is empty.')
                end

            end

            options.ColWidth = tempoptions.ColWidth;
        end
        
    end

    if isscalar(options.ColWidth)
        
        tempoptions.ColWidth = zeros(1, numel(options.Format));

        for formatstr_ind = 1:numel(options.Format)

            % get the number of decimal places specified
            [start_idx, end_idx, extents, matches] = regexp(options.Format{formatstr_ind}, '(?<=\.)\d+');

            if isempty(start_idx)

                % no numbers were found in the format spec, just use the
                % supplied width
                tempoptions.ColWidth(formatstr_ind) = options.ColWidth;

            else

                % some numbers were supplied, use the larger of the width
                % or the number of desired decimal places plus two (to
                % allow for leading number plus decimal point)
                tempoptions.ColWidth(formatstr_ind) = max(options.ColWidth, round(str2double(matches{1})) + 2);

            end

        end

        % replace scalar width with array of width values
        options.ColWidth = tempoptions.ColWidth;
        
    end
        
    if ~isempty(options.ColHeadings)
        
        [nRowCH,nColCH] = size(options.ColHeadings);
        if(nRowCH>1)
            error ('column headings can not have more than one row');
        end
        
        if ~iscellstr(options.ColHeadings)
            error('If not empty, options.ColHeadings must be a cell array of strings');
        end
        
        if(nColCH ~= nColD)
            error ('data must have same number of columns as headings');
        end
    
%         fmt = arrayfun(@(x) ['%',num2str(options.ColWidth(x)),'s |'], 1:nColD, 'UniformOutput', false);

        if options.WrapColHeadings
            
            tempcolheads = options.ColHeadings;
            nheaderrows = 1;
            for colheadind = 1:numel (options.ColHeadings)

                tempcolheads{colheadind} = TextWrapper.wraplines ( options.ColHeadings{colheadind}, ...
                                                                   'Width', options.ColWidth(colheadind) );

                nheaderrows = max ([nheaderrows, numel(tempcolheads{colheadind})]);

            end

            options.ColHeadings = repmat ({''}, nheaderrows, numel (options.ColHeadings));

            for colheadind = 1:size (options.ColHeadings, 2)

                options.ColHeadings(1:numel(tempcolheads{colheadind}), colheadind) = tempcolheads{colheadind};

            end
        
        end
        
        totalwid = numel (options.RowStart) ...
                    + sum (options.ColWidth) ...
                    + numel (options.ColSep)*numel (options.ColWidth) ...
                    + rhwid ...
                    + numel (options.RowEnding);
        horizontal_rule = repmat ( options.RuleChar, 1, totalwid);
        
        
        if options.TopRule

            tablestr = appendstr ( tablestr, '%s\n', horizontal_rule);

        end
        
        % Now loop through the column headings printing them out with the
        % desired column separator        
        for colheadrowind = 1:size (options.ColHeadings, 1)
            
            tablestr = appendstr (tablestr, '%s', options.RowStart);
            
            if ~isempty(options.RowHeadings)
                % TODO allow extra heading for column
                tablestr = appendstr (tablestr, ['%s',options.ColSep], repmat(' ', 1,rhwid));
            end
        
            for colheadcolind = 1:size (options.ColHeadings, 2)

                str = sprintf ( ['%',num2str(options.ColWidth(colheadcolind)),'s'], ...
                                options.ColHeadings{colheadrowind,colheadcolind} );

                if colheadcolind == size (options.ColHeadings, 2)
                    % If we are at the end of a row, don't print the column
                    % separator
                    tablestr = appendstr (tablestr, '%s', str(1:options.ColWidth(colheadcolind)) );
                else
                    % print the column header and column separator
                    tablestr = appendstr (tablestr, ['%s',options.ColSep], str(1:options.ColWidth(colheadcolind)) );
                end

            end

            tablestr = appendstr (tablestr, '%s\n', options.RowEnding);
        
        end

    end
    
    if options.ColHeaderRule
        
        tablestr = appendstr ( tablestr, '%s\n', horizontal_rule);
        
    end
    
    fmt = arrayfun(@(x) ['%',num2str(options.ColWidth(x)),options.Format{x}],1:nColD,'UniformOutput',false);
    
    for data_row_ind = 1:size(data,1)
        
        tablestr = appendstr (tablestr, '%s', options.RowStart);
        
        % first print a row header if necessary
        if ~isempty(options.RowHeadings)
            tablestr = appendstr (tablestr, ['%',num2str(rhwid),'s',options.ColSep], options.RowHeadings{data_row_ind});
        end
            
        % now loop through the data formatting and printing as appropriate
        for data_col_ind = 1:size(data,2)

            if iscell(data)
                % data is a cell array, possibly containing mixed data
                % types

                if ischar(data{data_row_ind,data_col_ind})
                    
                    str = sprintf(['%',num2str(options.ColWidth(data_col_ind)),'s'],data{data_row_ind,data_col_ind});
                    
                elseif isscalar(data{data_row_ind,data_col_ind})
                    
                    % write the number 
                    str = sprintf(fmt{data_col_ind},data{data_row_ind,data_col_ind});
                    
                    if length(str) > options.ColWidth(data_col_ind)

                        % convert to scientific notation as it doesn't fit
                        str =  sprintf(['%',num2str(options.ColWidth(data_col_ind)),'g'],data{data_row_ind,data_col_ind});

                        if length(str) > options.ColWidth(data_col_ind)
                            % fill space with #s as the number stil doesn't
                            % fit in
                            str = repmat('#', 1, options.ColWidth(data_col_ind));
                        end

                    end

                elseif isempty(data{data_row_ind,data_col_ind})
                    
                    % indicate an empty value 
                    str = sprintf(['%',num2str(options.ColWidth(data_col_ind)),'s'],'');
                    
                else
                    % we can only tabulate strings and scalars, so throw an
                    % error
                    error('CROZIER:displaytable:badcellcontents', ...
                          'each cell in cell array data must contain only individual strings or scalar values.')
                
                end
                
                % print the string
                tablestr = appendstr (tablestr, '%s', str(1:options.ColWidth(data_col_ind)));
                
                if data_col_ind < size(data,2)
                    % print column separator
                    tablestr = appendstr (tablestr, options.ColSep);
                end

            else
                
                % data is a numerical matrix

                str = sprintf(fmt{data_col_ind},data(data_row_ind,data_col_ind));

                if length(str) > options.ColWidth(data_col_ind)
                    
                    % convert to scientific notation as it doesn't fit
                    str =  sprintf(['%',num2str(options.ColWidth(data_col_ind)),'g'],data(data_row_ind,data_col_ind));
                    
                    if length(str) > options.ColWidth(data_col_ind)
                        % fill space with #s as the number doesn't fit in
                        str = repmat('#', 1, options.ColWidth(data_col_ind));
                    end
                    
                end
                
                if data_col_ind == size(data,2)
                    % do not print the last column separator at the end of
                    % a row
                    tablestr = appendstr (tablestr, '%s', str(1:options.ColWidth(data_col_ind)));
                else
                    tablestr = appendstr (tablestr, ['%s',options.ColSep], str(1:options.ColWidth(data_col_ind)));
                end
            end
            
        end
        
        % end of line so put in a new row, if there are more rows to be
        % printed
        if data_row_ind < size(data,1)
            tablestr = appendstr (tablestr, '%s\n', options.RowEnding);
        end
        
    end
    
    if options.BottomRule
        
        tablestr = appendstr ( tablestr, '%s', horizontal_rule);
        
    end

end


function tablestr = appendstr (tablestr, formatSpec, varargin)

    tablestr = [tablestr, sprintf(formatSpec, varargin{:})];
    
end
