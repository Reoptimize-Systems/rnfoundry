function help2txtfile (filename, fcn_name, varargin)
% gets fcn help, strips first two chars from each line and writes to a file
%
% Syntax
%
% help2txtfile (filename, fcn_name)
%
% Description
%
% help2txtfile gets the help for a function (or anything which the 'help'
% function can return help for). It strips first two characters from each
% line of the help and writes the result to a file. The first two
% characters are stripped from each line as typically these are two spaces,
% as when printing the help text, the first character (a % character is
% simply replaced with a space, and well-written help files have a space
% after the % character before the start of the text on every line). This
% behaviour can be changed by using the NStripChars option.
%
% Input
%
%  filename - file path to which the help will be written as a text file
%
%  fcn_name - the function or class or class method etc. for which to get
%   the help text. This can be any input accepted by the 'help' function
%
% Addtional arguments may be supplied as parameter-value pairs. The
% available options are:
%
%  'NStripChars' - scalar integer >= 0 giving the number of characters to
%    strip from the start of each line in the help text.
%
%
% See Also: cellstr2txtfile.m, str2txtfile.m
%

    options.NStripChars = 2;
    
    options = parse_pv_pairs (options, varargin);
    
    check.isPositiveScalarInteger (options.NStripChars, true, 'NStripChars', true);
    
    % create the readme files
    readme_txt = str2cellstr (help (fcn_name));
    
    % strip the first two white-space characters from each line
    for ind = 1:numel (readme_txt)

        readme_txt{ind} = readme_txt{ind}((1+options.NStripChars):end);
        
    end

    cellstr2txtfile (filename, readme_txt);

end