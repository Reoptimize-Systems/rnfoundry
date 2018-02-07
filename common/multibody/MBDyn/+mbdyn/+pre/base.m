classdef base < handle
    % base class for the mbdyn preprocessing tools
    
    
    properties (GetAccess = public, SetAccess = protected)
        
        label;
        type;
        
    end
    
    methods 
        function setLabel (self, label)
            self.label = label;
        end
    end
    
    methods (Static)
        
        function out = makeCellIfNot (in)
            if ~iscell (in)
                out = {in};
            else
                out = in;
            end
        end
        
        function ok = emptyOrCheck (checkfcn, checkvar, checkfcnargs)
            
            ok = true;
            if ~isempty (checkvar)
                ok = feval (checkfcn, checkvar, checkfcnargs{:});
            end
            
        end
        
        function ok = checkIsStructuralNode (node, throw)
            % checks if input is a mbdyn.pre.structuralNode
            %
            % Syntax
            %
            %  ok = checkIsStructuralNode (node, throw)
            %
            % Input
            %
            %  node - value to be tested if it is a mbdyn.pre.structuralNode
            %
            %  throw - logical flag determining whether an error is thrown
            %   by checkIsStructuralNode if node fails check
            %
            % Output
            %
            %  ok - logical flag indicating if check was passed
            %
            
            ok = true;
            if ~isa (node, 'mbdyn.pre.structuralNode')
                ok = false;
                
                if throw
                    error ('Input must be a mbdyn.pre.structuralNode object or subclass')
                end
            end
            
        end
        
        function ok = checkCartesianVector (vec, throw, name)
            % checks if input is a 3 element numeric column vector,
            % suitible for position, angular position, velocity, angular
            % velocity etc. It can also be a string: 'null' which
            % represents [ 0, 0, 0].
            %
            % Syntax
            %
            %  ok = checkCartesianVector (vec, throw)
            %
            % Input
            %
            %  vec - value to be tested if it is a 3 element numeric
            %    column vector, or the keyword 'null'
            %
            %  throw - logical flag determining whether an error is thrown
            %   by checkCartesianVector if vec fails check
            %
            % Output
            %
            %  ok - logical flag indicating if check was passed
            %
            
            if nargin < 3
                name = 'position, velocity or angular velocity';
            end
            
            ok = true;
            if ischar (vec)
                if ~strcmp (vec, 'null')
                    ok = false;
                
                    if throw
                        error ([name ' must be a 3 element numeric column vector, or ''null''']);
                    end
                end
            elseif ~((isnumeric (vec) && iscolumn (vec) && size (vec,1) == 3) || isempty (vec))
                
                ok = false;
                
                if throw
                    error ([name ' must be a 3 element numeric column vector, or ''null''']);
                end
            end
            
        end
        
        function ok = check3ElementNumericVector (vec, throw, name)
            % checks if input is a 3 element numeric vector
            %
            % Syntax
            %
            %  ok = check3ElementNumericVector (vec, throw)
            %
            % Input
            %
            %  vec - value to be tested if it is a 3 element numeric
            %    column vector, or the keyword 'null'
            %
            %  throw - logical flag determining whether an error is thrown
            %   by check3ElementNumericVector if vec fails check
            %
            %  name - optional string used to customise the error message.
            %   The error will be <name> must be a 3 element numeric vector.
            %   Default is 'input' if not supplied.
            %
            % Output
            %
            %  ok - logical flag indicating if check was passed
            %
            
            if nargin < 3
                name = 'input';
            end
            
            ok = true;
            if ~(isnumeric (vec) && isvector (vec) && numel (vec) == 3)
                
                ok = false;
                
                if throw
                    error ('%s must be a 3 element numeric vector', name)
                end
            end
            
        end
        
        function ok = checkNumericScalar (num, throw, name)
            % checks if input is a real scalar value
            %
            % Syntax
            %
            %  ok = checkNumericScalar (vec, throw)
            %  ok = checkNumericScalar (..., name)
            %
            % Input
            %
            %  num - value to be tested if it is a real numeric
            %    scalar
            %
            %  throw - logical flag determining whether an error is thrown
            %   by checkNumericScalar if num fails check
            %
            %  name - optional string used to customise the error message.
            %   The error will be <name> must be a scalar numeric value.
            %   Default is 'value' if not supplied.
            %
            %
            % Output
            %
            %  ok - logical flag indicating if check was passed
            %
            
            if nargin < 3
                name = 'value';
            end
            
            ok = true;
            if ~( isnumeric (num) && isscalar (num) && isreal (num) )
                
                ok = false;
                
                if throw
                    error ('%s must be a scalar numeric value', name);
                end
            end
            
        end
        
        function ok = checkScalarInteger (num, throw, name)
            % checks if input is a scalar integer value to machine precision
            %
            % Syntax
            %
            %  ok = checkScalarInteger (num, throw)
            %  ok = checkScalarInteger (..., name)
            %
            % Input
            %
            %  num - value to be tested if it is a scalar integer
            %
            %  throw - logical flag determining whether an error is thrown
            %   by checkScalarInteger if num fails check
            %
            %  name - optional string used to customise the error message.
            %   The error will be <name> must be a scalar integer (to
            %   machine precision). Default is 'value' if not supplied.
            %
            %
            % Output
            %
            %  ok - logical flag indicating if check was passed
            %
            
            if nargin < 3
                name = 'value';
            end
            
            ok = true;
            if ~( isnumeric (num) && isscalar (num) && isreal (num) && isint2eps (num) )
                
                ok = false;
                
                if throw
                    error ('%s must be a scalar integer (to machine precision)', name);
                end
            end
            
        end
        
        function ok = checkLogicalScalar (tfvalue, throw, name)
            % checks if input is a scalar logical value
            %
            % Syntax
            %
            %  ok = checkLogicalScalar (tfvalue, throw)
            %  ok = checkLogicalScalar (..., name)
            %
            % Input
            %
            %  tfvalue - value to be tested if it is a scalar logical value
            %
            %  throw - logical flag determining whether an error is thrown
            %   by checkLogicalScalar if tfvalue fails check
            %
            %  name - optional string used to customise the error message.
            %   The error will be "<name> must be a scalar logical
            %   (true/false) value". Default is 'value' if not supplied.
            %
            % Output
            %
            %  ok - logical flag indicating if check was passed
            %
            
            if nargin < 3
                name = 'value';
            end
            
            ok = true;
            if ~( islogical (tfvalue) && isscalar (tfvalue) )
                
                ok = false;
                
                if throw
                    error ('%s must be a scalar logical value', name);
                end
            end
            
        end
        
        function mat = getOrientationMatrix (om)
            % gets the raw 3x3 orientation matrix
            %
            % Syntax
            
            if isa (om, 'mbdyn.pre.orientmat')
                mat = om.orientationMatrix;
            else
                mat = om;
            end
            
        end
        
        function ok = checkOrientationMatrix (mat, throw, name)
            % checks if input is a valid orientation matrix
            %
            % Syntax
            %
            %  ok = checkOrientationMatrix (mat, throw)
            %  ok = checkOrientationMatrix (..., name)
            %
            % Input
            %
            %  mat - value to be tested if it is a valid orientation matrix
            %
            %  throw - logical flag determining whether an error is thrown
            %   by checkOrientationMatrix if mat fails check
            %
            %  name - optional string used to customise the error message.
            %   The error will be "<name> must be a 3 x 3 numeric matrix or
            %   mbdyn.pre.orientmat object". Default is 'orientation
            %   matrix' if not supplied.
            %
            % Output
            %
            %  ok - logical flag indicating if check was passed
            %
            
            if nargin < 3
                name = 'orientation matrix';
            end
            
            if isa (mat, 'mbdyn.pre.orientmat')
            
                mat = mbdyn.pre.base.getOrientationMatrix (mat);
                
            end
            
            ok = mbdyn.pre.base.check3X3Matrix (mat, false);
            
            if ~ok && throw
                error ('%s must be a 3 x 3 numeric matrix or mbdyn.pre.orientmat object', name);
            end

        end
        
        function ok = check3X3Matrix (mat, throw, name)
            % checks if input is a 3x3 numeric matrix
            %
            % Syntax
            %
            %  ok = check3X3Matrix (mat, throw)
            %  ok = check3X3Matrix (..., name)
            %
            % Input
            %
            %  mat - value to be tested if it is a 3x3 numeric matrix
            %
            %  throw - logical flag determining whether an error is thrown
            %   by check3X3Matrix if mat fails check
            %
            %  name - optional string used to customise the error message.
            %   The error will be "<name> must be a 3 x 3 numeric matrix".
            %   Default is 'input' if not supplied.
            %
            % Output
            %
            %  ok - logical flag indicating if check was passed
            %
            
            if nargin < 3
                name = 'input';
            end
            
            ok = true;
            if ~( (isnumeric (mat) && size (mat,1) == 3 && size (mat,2) == 3) ...
                    || isempty (mat) )
                
                ok = false;
                
                if throw
                    error ('%s must be a 3 x 3 numeric matrix', name)
                end
            end
            
        end
        
        function ok = checkOrientationDescription (odesc, throw, name)
            % checks if the orientation description choice is valid
            %
            % Syntax
            %
            %  ok = checkOrientationDescription (odesc, throw)
            %
            % Input
            %
            %  odesc - value to be tested if it is a valid orientation
            %   description string. Valid values are: 'euler123',
            %   'euler313', 'euler321', 'orientation vector', 
            %   'orientation matrix
            %
            %  throw - logical flag determining whether an error is thrown
            %   by checkOrientationDescription if odesc fails check
            %
            %  name - optional string used to customise the error message.
            %   The error will be "<name> must be one of: 'euler123' |
            %   'euler313' | 'euler321' | 'orientation vector' |
            %   'orientation matrix'". Default is 'Orientation description'
            %   if not supplied.
            %
            % Output
            %
            %  ok - logical flag indicating if check was passed
            %
            
            if nargin < 3
                name = 'Orientation description';
            end
            
            if ~ischar (odesc)
                error ('Orientation description must be a char array')
            end
            
            ok = mbdyn.pre.base.checkAllowedStringInputs ( odesc, {'euler123', ...
                                                                   'euler313', ...
                                                                   'euler321', ...
                                                                   'orientation vector', ...
                                                                   'orientation matrix'}, ...
                                                           throw, name);
            
        end
        
        function ok = checkAllowedStringInputs (input, allowedstrs, throw, inputname)
            % checks is string is one of a set of allowed values
            %
            % Syntax
            %
            %  ok = checkAllowedStringInputs (input, allowedstrs, throw)
            %  ok = checkAllowedStringInputs (..., inputname)
            %
            % Input
            %
            %  input - value to be tested if it is a valid string. The
            %   string will be compared to the list of valid strings in
            %   allowedstrs.
            %
            %  allowedstrs - cell string array of valid strings for
            %    comparison with input.
            %
            %  throw - logical flag determining whether an error is thrown
            %   by checkOrientationDescription if input fails check
            %
            %  inputname - (optional) string with name to use in error
            %    thrown by checkAllowedStringInputs (if 'throw' is true).
            %    The error will be "<name> must be one of: 'string 1' |
            %   'string 2' | 'string 3'" where allowedstrs is {'string 1',
            %   'string 2', 'string 3'}.
            %
            % Output
            %
            %  ok - logical flag indicating if check was passed
            %
            
            if nargin < 4
                inputname = 'input';
            end
            
            if ~iscellstr (allowedstrs)
                error ('allowedstrs must be a cell array of strings containing a list of allowed values for the input');
            end
            
            assert (ischar (input), '%s must be a char array', inputname);
            
            ok = true;
            if ~any(strcmp (input, allowedstrs))
                ok = false;
            end
            
            if throw && ~ok
               
                error ('%s must be one of: %s ', inputname, sprintf (['''%s''', repmat(' | ''%s''', 1, numel(allowedstrs)-1)], allowedstrs{:}));
                
            end
            
        end
        
        function ok = checkTplDriveCaller (drv, throw, inputname)
            % checks if input is a mbdyn.pre.driveCaller object
            %
            % Syntax
            %
            %  ok = checkTplDriveCaller (drv, throw)
            %  ok = checkTplDriveCaller (drv, throw, inputname)
            %
            % Input
            %
            %  drv - value to be tested if it is a valid
            %   mbdyn.pre.driveCaller object or an object derived from
            %   this class
            %
            %  throw - logical flag determining whether an error is thrown
            %   by checkTplDriveCaller if drv fails check
            %
            %  inputname - (optional) string containing name of variable.
            %   Error message will be customised to include the string. If
            %   not supplied 'input' is used.
            %
            %
            % Output
            %
            %  ok - logical flag indicating if check was passed
            %
            
            if nargin < 3
                inputname = 'input';
            end
            
            ok = true;
            if (~isempty (drv)) && (~isa (drv, 'mbdyn.pre.driveCaller'))
                ok = false;
            end
            
            if ~ok && throw
                error ('%s must be a mbdyn.pre.driveCaller object', inputname);
            end

        end
        
        function ok = checkDrive (drv, throw, inputname)
            % checks if input is a mbdyn.pre.drive object
            %
            % Syntax
            %
            %  ok = checkDrive (drv, throw)
            %  ok = checkDrive (drv, throw, inputname)
            %
            % Input
            %
            %  drv - value to be tested if it is a valid
            %   mbdyn.pre.drive object or an object derived from
            %   this class
            %
            %  throw - logical flag determining whether an error is thrown
            %   by checkDrive if drv fails check
            %
            %  inputname - (optional) string containing name of variable.
            %   Error message will be customised to include the string. If
            %   not supplied 'input' is used.
            %
            %
            % Output
            %
            %  ok - logical flag indicating if check was passed
            %
            
            if nargin < 3
                inputname = 'input';
            end
            
            ok = true;
            if (~isempty (drv)) && (~isa (drv, 'mbdyn.pre.drive'))
                ok = false;
            end
            
            if ~ok && throw
                error ('%s must be a mbdyn.pre.drive object', inputname);
            end

        end
       
        function str = commaSepList (varargin)
            % generates a comma separated list from input variables
            
            str = '';
            
            for ind = 1:numel (varargin)
                
                if isnumeric (varargin{ind})
                    
                    mat = varargin{ind};
                    
                    if numel (mat) > 1 && size(mat, 1) ~= 1 && size (mat, 2) ~= 1
                        % matrix, not vector or scalar
                        if ind > 1
                            % write matrix on it's own lines for clarity
                            str = [ str, sprintf('\n') ];
                        end
                        
                        str = [ str, mbdyn.pre.base.writeMatrix(mat)];
                        
                        % add a comma and newline to end of matrix so rest
                        % of comma separated list continues on next line
                        str = [ str, sprintf(',\n') ];
                        
                    else
                        % handle vectors and scalars
                        for matind = 1:numel (mat)
                            numstr = mbdyn.pre.base.formatNumber (mat(matind));
                            str = [ str, sprintf('%s, ', numstr) ];
                        end
                        
                    end
                    
                elseif ischar (varargin{ind})
                    
                    str = [ str, varargin{ind}, ', '];
                    
                elseif isa (varargin{ind}, 'mbdyn.pre.orientmat')
                    
                    str = [ str, mbdyn.pre.base.writeMatrix(mbdyn.pre.base.getOrientationMatrix (varargin{ind}) )];
                        
                    % add a comma and newline to end of matrix so rest
                    % of comma separated list continues on next line
                    str = [ str, sprintf(',\n') ];
                    
                end
                
            end
            
            % strip last comma and space (or newline if matrix was last)
            str(end-1:end) = [];
            
        end
        
        function str = writeMatrix (mat, usematr)
            % format a matrix to a string for an mbdyn input file
            
            if nargin < 2
                usematr = true;
            end
            
            if usematr
                str = sprintf('matr,\n');
            else
                str = '';
            end

            % write out matrix row-wise for mbdyn
            for rowind = 1:size (mat, 1)
                for colind = 1:size (mat, 2)
                    numstr = mbdyn.pre.base.formatNumber (mat(rowind,colind));
                    str = [ str, sprintf('%s, ', numstr) ];
                end
                % add a newline after each row
                str = sprintf('%s\n', str);
            end
            
            % strip last comma, space and newline
            str(end-2:end) = [];
            
        end
        
        function numstr = formatNumber (num)
            % fomats a decimal number for pretty output to mbdyn file
            %
            % If it is an integer (to machine precision), it will be output
            % with one decimal place. If a float it will be output to 14
            % decimal places, but with any trailing zeros stripped to
            % shorten it.
            
            if isint2eps (num)
                numstr = sprintf ('%d.0', num);
            else
                numstr = sprintf ('%.18f', num);
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
        
        function numstr = formatInteger (num)
            % fomats a non-decimal integer number for output to mbdyn file
            %
            
            if isint2eps (num)
                numstr = sprintf ('%d', num);
            else
                error ('Supplied number is not an integer (to machine precision), cannot format');
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
        
        function nodeType = getNodeType (node)
            
            if isa (node, 'mbdyn.pre.structuralNode6dof') 
                
                nodeType = 'structural';
                
            elseif isa (node, 'mbdyn.pre.structuralNode3dof') 

                nodeType = 'structural';
                
            elseif isa (node, 'mbdyn.pre.abstractNode') 
                
                nodeType = 'abstract';

            elseif isa (node, 'mbdyn.pre.electricNode') 
                
                nodeType = 'electric';
                
            elseif isa (node, 'mbdyn.pre.hydraulicNode') 
                
                nodeType = 'hydraulic';
                
            elseif isa (node, 'mbdyn.pre.parameterNode') 
                
                nodeType = 'parameter';

            elseif isa (node, 'mbdyn.pre.thermalNode') 
                
                nodeType = 'thermal';
                
            else
                error ('Supplied node type is not yet implemented in preprocessor')
            end 
            
        end
        
        function elementType = getElementType (element)
            
            
            if isa (element, 'mbdyn.pre.automaticStructural') 
                
                elementType = 'automatic structural';
                
            elseif isa (element, 'mbdyn.pre.beam')
                
                elementType = 'beam';
                
            elseif isa (element, 'mbdyn.pre.body')
                
                elementType = 'body';
                
            elseif isa (element, 'mbdyn.pre.couple')
                
                elementType = 'couple';
                
            elseif isa (element, 'mbdyn.pre.gravity')
                
                elementType = 'gravity';
                
            elseif isa (element, 'mbdyn.pre.joint')
                
                elementType = 'joint';
                
            elseif isa (element, 'mbdyn.pre.jointRegularization')
                
                elementType = 'joint regularization';
                
            elseif isa (element, 'mbdyn.pre.plate')
                
                elementType = 'plate';
                
            elseif isa (element, 'mbdyn.pre.electric')
                
                elementType = 'electric';
                
            elseif isa (element, 'mbdyn.pre.hydraulic')
                
                elementType = 'hydraulic';
                
            elseif isa (element, 'mbdyn.pre.bind')
                
                elementType = 'bind';
                
            elseif isa (element, 'mbdyn.pre.bulk')
                
                elementType = 'bulk';
                
            elseif isa (element, 'mbdyn.pre.force')
                
                elementType = 'force';
                
            elseif isa (element, 'mbdyn.pre.genel')
                
                elementType = 'genel';
                
            elseif isa (element, 'mbdyn.pre.loadable')
                
                elementType = 'loadable';
                
            elseif isa (element, 'mbdyn.pre.userDefined')
                
                elementType = 'user defined';
                
            else
                error ('Supplied node type is not yet implemented in preprocessor')
            end 
            
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
            %  'ShowLabels' - flag determining whether labels will be added 
            %    for each reference. The labels will be the reference name
            %    (from the name property of the reference), or, if this is
            %    empty, a label based on the index of the reference in the
            %    input cell aray, e.g. Reference 1 for the first label,
            %    Reference 2 for the second etc. Default is true.
            %
            %
            
            options.PlotAxes = [];
            options.Title = true;
            options.DrawGlobal = true;
            options.Scale = 1;
            options.ShowLabels = true;
            options.XLim = [];
            options.YLim = [];
            options.ZLim = [];
           
            options = parse_pv_pairs (options, varargin);
            
            mbdyn.pre.base.checkLogicalScalar (options.Title, true, 'Title' );
            mbdyn.pre.base.checkLogicalScalar (options.DrawGlobal, true, 'DrawGlobal' );
            mbdyn.pre.base.checkLogicalScalar (options.ShowLabels, true, 'ShowLabels' );
            mbdyn.pre.base.checkNumericScalar (options.Scale, true, 'Scale' );
            
            hax = draw ( refs{1}, ...
                         'PlotAxes', options.PlotAxes, ...
                         'Title', options.Title, ...
                         'DrawGlobal', options.DrawGlobal, ...
                         'Scale', options.Scale);
           
            if options.ShowLabels
                if isempty (refs{1}.name)
                    label = sprintf ('Reference %d', 1);
                else
                    label = refs{1}.name;
                end
            else
                label = '';
            end
                
            text ( refs{1}.pos(1) + 0.1*options.Scale, ...
                   refs{1}.pos(2), ...
                   refs{1}.pos(3) + 0.1*options.Scale, ...
                   label, ...
                   'Interpreter', 'none');
                             
            for ind = 2:numel (refs)
                draw ( refs{ind}, ...
                       'PlotAxes', hax, ...
                       'Title', false, ...
                       'DrawGlobal', false, ...
                       'Scale', options.Scale );
                
                if options.ShowLabels
                    if isempty (refs{ind}.name)
                        label = sprintf ('Reference %d', ind);
                    else
                        label = refs{ind}.name;
                    end
                else
                    label = '';
                end
                
                text ( refs{ind}.pos(1) + 0.1*options.Scale, ...
                       refs{ind}.pos(2), ...
                       refs{ind}.pos(3) + 0.1*options.Scale, ...
                       label, ...
                       'Interpreter', 'none');
            end
            
            if ~isempty (options.XLim)
                set (hax, 'XLim', options.XLim);
            end
            
            if ~isempty (options.YLim)
                set (hax, 'YLim', options.YLim);
            end
            
            if ~isempty (options.ZLim)
                set (hax, 'ZLim', options.ZLim);
            end
            
        end
        
    end
    
end


function result = isint2eps(X)
% isint2eps determines if the numbers in a matrix are integers to the limit
% of the machine floating-point relative accuracy for those values

    theMod = mod(X,1);
    
    result = theMod <= eps(theMod);

end

