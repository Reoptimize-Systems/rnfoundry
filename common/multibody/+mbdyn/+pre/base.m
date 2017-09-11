classdef base < handle
    % base class for the mbdyn preprocessing tools
    
    
    properties (GetAccess = public, SetAccess = protected)
        
        label;
        type;
        
    end
    
    methods (Static)
        
        function out = makeCellIfNot (in)
            if ~iscell (in)
                out = {in};
            else
                out = in;
            end
        end
        
        function ok = checkIsStructuralNode (node, throw)
            
            ok = true;
            if ~isa (node, 'mbdyn.pre.structuralNode')
                ok = false;
                
                if throw
                    error ('Input must be a mbdyn.pre.structuralNode object or subclass')
                end
            end
            
        end
        
        function ok = checkCartesianVector (vec, throw)
            
            ok = true;
            if ischar (vec)
                if ~strcmp (vec, 'null')
                    ok = false;
                
                    if throw
                        error ('position, velocity or angular velocity must be a 3 element numeric column vector, or ''null''')
                    end
                end
            elseif ~((isnumeric (vec) && iscolumn (vec) && size (vec,1) == 3) || isempty (vec))
                
                ok = false;
                
                if throw
                    error ('position, velocity or angular velocity must be a 3 element numeric column vector, or ''null''')
                end
            end
            
        end
        
        function mat = getOrientationMatrix (om)
            
            if isa (om, 'mbdyn.pre.orientmat')
                mat = om.orientationMatrix;
            else
                mat = om;
            end
            
        end
        
        function ok = checkOrientationMatrix (mat, throw)
            
            mat = mbdyn.pre.base.getOrientationMatrix (mat);
            
            ok = mbdyn.pre.base.check3X3Matrix (mat, false);
            
            if ~ok && throw
                error ('orientation matrix must be a 3 x 3 numeric matrix or mbdyn.pre.orientmat object');
            end

        end
        
        function ok = check3X3Matrix (mat, throw)
            
            ok = true;
            if ~((isnumeric (mat) && size (mat,1) == 3 && size (mat,2) == 3) || isempty (mat))
                
                ok = false;
                
                if throw
                    error ('input must be a 3 x 3 numeric matrix')
                end
            end
            
        end
        
        function ok = checkOrientationDescription (odesc, throw)
            % checks if the orientation description choice is valid
            
            if ~ischar (odesc)
                error ('Orientation description must be a char array')
            end
            
            ok = mbdyn.pre.base.checkAllowedStringInputs ( odesc, {'euler123', ...
                                                                   'euler313', ...
                                                                   'euler321', ...
                                                                   'orientation vector', ...
                                                                   'orientation matrix'}, ...
                                                           throw, 'Orientation description');
            
        end
        
        function ok = checkAllowedStringInputs (input, allowedstrs, throw, inputname)
            
            if nargin < 4
                inputname = 'input';
            end
            
            if ~iscellstr (allowedstrs)
                error ('allowed must be a cell array of strings containing a list of allowed values for the input');
            end
            
            ok = true;
            if ~any(strcmp (input, allowedstrs))
                ok = false;
            end
            
            if throw && ~ok
               
                error ('%s must be one of: %s ', inputname, sprintf (['''%s''', repmat(' | ''%s''', 1, numel(allowedstrs)-1)], allowedstrs{:}));
                
            end
            
        end
       
        function str = commaSepList (varargin)
            % generates a comma separated list from input ariables
            
            str = '';
            
            for ind = 1:numel (varargin)
                
                if isnumeric (varargin{ind})
                    
                    % get the transpose of the matrix as it must be written
                    % out row-wise for mbdyn
                    mat = varargin{ind}';
                    
                    if numel (mat) > 1 && size(mat, 1) == size (mat, 2)
                        str = [ str, 'matr, '];
                    end
                    
                    for matind = 1:numel (mat)
                        numstr = mbdyn.pre.base.formatNumber (mat(matind));
                        str = [ str, sprintf('%s, ', numstr) ];
                    end
                    
                elseif ischar (varargin{ind})
                    
                    str = [ str, varargin{ind}, ', '];
                    
                end
                
            end
            
            % strip last comma and space
            str(end-1:end) = [];
            
        end
        
        function numstr = formatNumber (num)
            % fomats a decimal number for pretty output to mbdyn file
            %
            % If it is an integer (to machine precision), it will be output
            % with one decimal place. If a float it will be output to 14
            % decimal places, but with any trailing zeros stripped to
            % shorten it.
            
            if isint2eps (num)
                numstr = sprintf('%d.0', num);
            else
                numstr = sprintf('%.18f', num);
                % strip trailing zeros from decimals
                n = numel (numstr);
                while numstr(n) ~= '.'
                    if numstr(n) == '0'
                        numstr(n) = [];
                    else
                        break;
                    end
                    n = n - 1;
                end
            end
            
        end
        
        function str = addOutputLine (oldstr, str, indentlevel, finalcomma, comment)
            % add a line of output to the output string representing the
            % contents of an MBDyn input file. 
            %
            % Syntax
            %
            %  str = addOutputLine (oldstr, str)
            %  str = addOutputLine (..., indentlevel)
            %  str = addOutputLine (..., indentlevel, finalcomma)
            %  str = addOutputLine (..., indentlevel, finalcomma, comment)
            %
            % Description
            %
            % addOutputLine is a utility function to help create readable,
            % structured MBDyn input files. It is intended to be used in
            % the generateOutputString methods of mbdyn preprocessing
            % objects like elements and nodes. It assists with keeping
            % consistant indentation levels in files. See the examples
            % section for examples of it's use.
            %
            % Input
            %
            %  oldstr - existing string to which new line is to be
            %   appended.  An empty string is acceptable.
            %
            %  str - new string containing text to be appended to oldstr
            %
            %  indentlevel - (optional) scalar integer. The indentation
            %   level of the new line to be added. One indentation level is
            %   four spaces. If the new string contains newline characters,
            %   the same indentation will be inserted before every new
            %   line. The indentation level is absolute, not relative to
            %   oldstr. Default value is zero if not supplied.
            %
            %  finalcomma - (optional) scalar logical. Flag determining
            %    whether to append a comma chachter (',') to the end of the
            %    new line. Default is true if not supplied.
            %
            %  comment - (optional) string which will be inserted as a
            %    comment after newstring.
            %
            % Output
            %
            %  str - string comprised of "oldstr", a newline character, and
            %   then the new string in "str", possibly modified to have the
            %   specified indentation level
            %
            % Examples
            %
            %
            % >> mbdyn.pre.base.addOutputLine ('old string', 'new line of output')
            % 
            % ans =
            % 
            % old string
            % new line of output,
            % 
            % >> mbdyn.pre.base.addOutputLine ('old string', 'new line of output', 1)
            % 
            % ans =
            % 
            % old string
            %     new line of output,
            % 
            % >> mbdyn.pre.base.addOutputLine ('old string', 'new line of output', 2)
            % 
            % ans =
            % 
            % old string
            %         new line of output,
            %
            % % 2 lines of output
            % >> mbdyn.pre.base.addOutputLine ('old string', sprintf('2 new lines\nof output'), 1)
            % 
            % ans =
            % 
            % old string
            %     2 new lines
            %     of output,
            %
            % % 2 lines of output, but don't add a final comma
            % >> mbdyn.pre.base.addOutputLine ('old string', sprintf('2 new lines\nof output'), 1, false)
            % 
            % ans =
            % 
            % old string
            %     2 new lines
            %     of output
            %
            % % using comment
            % >> mbdyn.pre.base.addOutputLine ('old string', sprintf('2 new lines\nof output'), 1, false, 'this is a comment')
            % 
            % ans =
            % 
            % old string
            %     2 new lines
            %     of output # this is a comment
            %
            %
            % See also: mbdyn.pre.base.formatNumber, 
            %           mbdyn.pre.base.commaSepList
            
            if nargin < 4
                finalcomma = true;
            end
            
            if nargin < 3
                indentlevel = 0;
            end
            
            prefix = repmat ('    ', 1, indentlevel);
            
            % add the indentation to any existing newline characters in the
            % new string to be appended so it's indentation level is
            % consistant
            str = strrep (str, sprintf ('\n'), sprintf ('\n%s', prefix));
            
            % joint the string to the existing string with a newline and
            % the appropriate indentation
            str = sprintf ('%s\n%s%s', oldstr, prefix, str);
            
            if finalcomma
                str = [str, ','];
            end
            
            if nargin > 4
                str = [ str, ' # ', comment ];
            end
            
        end
        
        function M = mbdynOrient2Matlab (M)

            % matlabs angles are clockwise 

%             M = [ M(1:3,1:3).', M(1:3,4); ...
%                   0, 0, 0, 1 ];
%             M = [ M(1:3,1:3), M(1:3,4); ...
%                   0, 0, 0, 1 ];
        end
        
        function drawReferences (refs, varargin)
            % draw a collection of mbdyn.pre.reference objects
            %
            % Syntax
            %
            % drawReferences (refs, 'Parameter', value)
            %
            % Description
            %
            % Draws a collection of reference objects in one figure. Each
            % reference is represented as a three-axis coordinate system.
            %
            % Input
            %
            %  refs - cell array of 2 or more mbdyn.pre.reference objects
            %
            % Additional optional arguments can be provided using
            % parameter-value pairs. The available options are:
            %
            %  'PlotAxes' - axes in which to create the plot. If not
            %    supplied a new figure and axes will be created.
            %
            %  'Title' - flag determining whether to add a title to the
            %    plot. Default is true.
            %
            %  'DrawGlobal' - flag determining whether the global axes will
            %    be drawn in the plot for reference. Default is true
            %
            %  'Scale' - scalar value. References are drawn with a default 
            %    size, this option may be used to adjust this by scaling
            %    the length of their axes plots up or down. Default is 1,
            %    no scaling.
            %
            %
            
            options.PlotAxes = [];
            options.Title = true;
            options.DrawGlobal = true;
            options.Scale = 1;
           
            options = parse_pv_pairs (options, varargin);
            
            hax = draw ( refs{1}, ...
                         'PlotAxes', options.PlotAxes, ...
                         'Title', true, ...
                         'DrawGlobal', true, ...
                         'Scale', options.Scale);
           
            text ( refs{1}.pos(1) + 0.1*options.Scale, ...
                   refs{1}.pos(2), ...
                   refs{1}.pos(3) + 0.1*options.Scale, ...
                   refs{1}.name, ...
                   'Interpreter', 'none');
                             
            for ind = 2:numel (refs)
                draw ( refs{ind}, ...
                       'PlotAxes', hax, ...
                       'Title', false, ...
                       'DrawGlobal', false, ...
                       'Scale', options.Scale );
                
                
                if isempty (refs{ind}.name)
                    label = sprintf ('Reference %d', ind);
                else
                    label = refs{ind}.name;
                end
                
                text ( refs{ind}.pos(1) + 0.1*options.Scale, ...
                       refs{ind}.pos(2), ...
                       refs{ind}.pos(3) + 0.1*options.Scale, ...
                       label, ...
                       'Interpreter', 'none');
            end
            
        end
        
    end
    
end