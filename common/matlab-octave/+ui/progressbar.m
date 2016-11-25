classdef progressbar < handle
    
    properties
        % Vizualization parameters
        strPercentageLength;   %   Length of percentage string (must be >5)
        strDotsMaximum;%   The total number of dots in a progress bar
        textLeader;
        firstRun;
        strCR;
        lastProgress;
        
    end
    
    methods
        
%         function self = progressbar (varargin)
%             
%             init (self, varargin)
%             
%         end
        
        function init (self, varargin)
            
            options.StrPercentageLength = 10;
            options.StrDotsMaximum = 50;
            options.TextLeader = '';
            
            options = parse_pv_pairs (options, varargin);
            
            self.strPercentageLength = options.StrPercentageLength;
            self.strDotsMaximum = options.StrDotsMaximum;
            self.textLeader = options.TextLeader;
            
            self.firstRun = true;
            self.strCR = '';
            self.lastProgress = 0;
            
        end
        
        function done (self, str)
            
            if nargin < 2
                str = '';
            end
            
            fprintf([str, '\n']);
            
        end
        
        function dispProgress (self, c)
            
            % Progress bar - normal progress
            c = floor (c);
            
            percentageOut = [num2str(c) '%%'];
            
            percentageOut = [ percentageOut, ...
                              repmat(' ', 1, self.strPercentageLength-length(percentageOut)-1)];
            
            nDots = floor (c/100*self.strDotsMaximum);
            
            dotOut = ['[' repmat('.',1,nDots) repmat(' ',1,self.strDotsMaximum-nDots) ']'];
            
            strOut = [self.textLeader, percentageOut, dotOut];

            % Print it on the screen
            if self.firstRun
                % Don't do carriage return during first run
                fprintf (strOut);
            else
                % Do it during all the other runs
                fprintf ([self.strCR, strOut]);
            end

            % Update carriage return
            self.strCR = repmat ('\b', 1, length(strOut)-1);
            
            self.firstRun = false;
            
            self.lastProgress = c;
            
        end
        
    end
    
end