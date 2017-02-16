classdef base < handle
    
    
    properties
        
        label;
        
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
            
            ok = true;
            
            if ~any(strcmp (odesc, {'euler123', 'euler313', 'euler321', 'orientation vector', 'orientation matrix'}))
                
                ok = false;
                
            end
            
            if throw && ~ok
               
                error ('Orientation description must be: euler123 | euler313 | euler321 | orientation vector | orientation matrix');
                
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
                    
                    for matind = 1:numel (mat)
                        str = [ str, sprintf('%g, ', mat(matind)) ];
                    end
                    
                elseif ischar (varargin{ind})
                    
                    str = [ str, varargin{ind}, ', '];
                    
                end
                
            end
            
            % strip last comma and space
            str = sprintf ( [str, '\b\b'] );
            
        end
        
        function str = addOutputLine (oldstr, str, indentlevel, finalcomma, comment)
            
            if nargin < 3
                finalcomma = true;
            end
            
            prefix = repmat ('    ', 1, indentlevel);
            
            str = sprintf ('%s\n%s%s', oldstr, prefix, str);
            
            if finalcomma
                str = [str, ','];
            end
            
            if nargin > 4
                str = [ str, ' # ', comment ];
            end
            
        end
        
    end
    
    
end