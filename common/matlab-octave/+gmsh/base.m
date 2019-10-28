classdef base < handle
    
    properties
        
        tag;
        geoName;
        transformObject;
        netCDFName;
        uid;
        
    end
    
    properties (GetAccess = protected, SetAccess = protected)
       
        shapeData;
        shapeParameters;
        shapeObjects;   % object for the element shape drawing
        needsRedraw; % track whether object needs to be redrawn
        drawAxesH; % handle to figure for plotting
        drawColour;

    end
    
    methods
        
        function str = generateGeoFileStr (self)
            
            str = sprintf ('%s (%s) =', self.geoName, self.formatInteger (self.tag));          
            
        end

        function self = base ()
           
            self.uid = round (rand ()*1e12);
            
        end
        
        function test = eq (a, b)
            
            if isoctave
                test = a.uid == b.uid;
            else
                test = eq@handle (a, b);
            end
            
        end
        
        function setTag (self, tag)
            self.tag = tag;
        end
        
        function checkAxes (self, hax)
            % checks if there is a valid set of axes to plot to, and if
            % not, creates them
            
            % try to figure out if there is a valid set of axes to plot to
            if isempty (hax)
                
                % plot in the existing axes if possible as no new axes
                % supplied
                if self.isAxesHandle (self.drawAxesH)
                    % use existing axes
                    
                    if ~ishghandle (self.drawAxesH)
                        % axes no longer exists, or figure has been closed
                        self.drawAxesH = [];
                        self.deleteAllDrawnObjects ();
                        self.transformObject = [];
                        self.needsRedraw = true;
                        % make a new set of axes to draw to
                        self.makeAxes ();
                    end
                    
                elseif isempty (self.drawAxesH)
                    % make figure and axes
                    self.makeAxes ();
                    self.needsRedraw = true;
                    
                else
                    error ('drawAxesH property is not empty or an axes handle');
                end
            
            elseif self.isAxesHandle (hax)
                % an axes has been supplied, so we plot to this new axes
                
                if isoctave
                    drawnow ();
                end
                
                if ~ishghandle (hax)
                    error ('provided axes object is not valid');
                end
                self.drawAxesH = hax;
                % abandon existing shape objects and transform object
                % TODO: should we delet objects in old axes?
                self.shapeObjects = {};
                self.transformObject = [];
                % we need to redraw, as we're plotting in a different set
                % of axes
                self.needsRedraw = true;
                
            else
                
                error ('hax must be an axes handle or empty');
                
            end
            
            if isempty (self.transformObject) || ~ishghandle (self.transformObject)
                self.transformObject = hgtransform (self.drawAxesH);
            end
            
        end
        
        function makeAxes (self)
            % create axes and transform object
            
            hfig = figure;
            
            if isoctave
                self.drawAxesH = axes;
            else
                self.drawAxesH = axes (hfig);
            end
            
            if ~isempty (self.transformObject) && ishghandle (self.transformObject)
                delete (self.transformObject);
            end
            self.transformObject = [];
            self.needsRedraw = true;
            
        end
        
        function deleteAllDrawnObjects (self)
            
            for ind = 1:numel (self.shapeObjects)
                if ~isempty (self.shapeObjects{ind}) ...
                        && ishghandle (self.shapeObjects{ind})

                    delete (self.shapeObjects{ind});

                end
            end
            self.shapeObjects = {};
            
        end
        
        
        
    end
    
    methods (Static)
        
        function numstr = formatNumber (num)
            % fomats a decimal number for pretty output to gmsh geo file
            %
            % If it is an integer (to machine precision), it will be output
            % with one decimal place. If a float it will be output to 14
            % decimal places, but with any trailing zeros stripped to
            % shorten it.
            
            if check.isint2eps (num)
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
        
        function numstr = formatInteger (varargin)
            % fomats a non-decimal integer number for output to gmsh geo file
            %
            
            numstr = cell(1, numel(varargin));
            
            for ind = 1:numel (varargin)
                num = varargin{ind};
                
                if check.isint2eps (num)
                    numstr{ind} = sprintf ('%d', num);
                else
                    error ('Supplied number is not an integer (to machine precision), cannot format');
                end
            end
            
            if numel (varargin) == 1
                numstr = numstr{1};
            end
            
        end
        
        
        function str = commaSepList (varargin)
            % generates a comma separated list from input variables
            
            str = '';
            
            for ind = 1:numel (varargin)
                
                if isnumeric (varargin{ind})
                    
                    mat = varargin{ind};
                    
                    % handle vectors and scalars
                    for matind = 1:numel (mat)
                        numstr = mbdyn.pre.base.formatNumber (mat(matind));
                        str = [ str, sprintf('%s, ', numstr) ];
                    end
                    
                elseif ischar (varargin{ind})
                    
                    str = [ str, varargin{ind}, ', '];
                    
                elseif isa (varargin{ind}, 'gmsh.base')
                    
                    for obj_ind = 1:numel (varargin{ind})
                    
                        str = [ str, varargin{ind}.generateMBDynInputString() ];

                        % add a comma and newline to end of matrix so rest
                        % of comma separated list continues on next line
                        str = [ str, sprintf(',\n') ];
                    
                    end
                    
                end
                
            end
            
            % strip last comma and space (or newline if matrix was last)
            str(end-1:end) = [];
            
        end
        
    end
    
end