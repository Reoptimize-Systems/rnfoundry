classdef base < handle
    
    
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
            if ~((isnumeric (vec) && iscolumn (vec) && size (vec,1) == 3) || isempty (vec))
                
                ok = false;
                
                if throw
                    error ('position, velocity or angular velocity must be a 3 element numeric column vector')
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
                        if isint2eps (mat(matind))
                            str = [ str, sprintf('%d.0, ', mat(matind)) ];
                        else
                            numstr = sprintf('%.14f', mat(matind));
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
                            str = [ str, sprintf('%s, ', numstr) ];
                        end
                    end
                    
                elseif ischar (varargin{ind})
                    
                    str = [ str, varargin{ind}, ', '];
                    
                end
                
            end
            
            % strip last comma and space
            str(end-1:end) = [];
            
        end
        
        function str = addOutputLine (oldstr, str, indentlevel, finalcomma, comment)
            
            if nargin < 3
                finalcomma = true;
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

            M = [ M(1:3,1:3).', M(1:3,4); ...
                  0, 0, 0, 1 ];
        end
        
    end
    
    
end