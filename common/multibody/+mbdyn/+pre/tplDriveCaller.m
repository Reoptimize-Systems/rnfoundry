classdef tplDriveCaller < mbdyn.pre.driveCaller
    
    properties
        
        type;
        driveCaller;
    end
    
    methods
        
        function self = tplDriveCaller (type, drivecallers, varargin)
            
            switch type
                
                case 'null'
                    
                case 'single'
                    options.Entity = [];
                    
                case 'component'
                    
                otherwise
                    error ('tplDriveCaller type must be null | single | component');
                    
            end
            
            options = parse_pv_pairs (options, varargin);
            
            
            switch type
                
                case 'null'
                    
                case 'single'
                    
                    
                case 'component'
                    
            end 
            
            self.type = type;
            self.driveCallers = drivecallers;
            
        end
        
        
        function str = generateOutputString (self)

            % delete newline character and space from start
            str = self.addOutputLine (str , 'revolute pin', 0, true);
            
            str = self.addOutputLine (str, sprintf('%d', self.node.label), 2, true, 'node label');
            
            out = self.makeCellIfNot (self.relativeOffset);
            str = self.addOutputLine (str, self.commaSepList ('position', out{:}), 3, true, 'node relative position' );
            
            if ~isempty (self.nodeRelativeOrientation)
                out = self.makeCellIfNot (self.nodeRelativeOrientation);
                str = self.addOutputLine (str, self.commaSepList ('orientation', out{:}), 3, true, 'node relative orientation');
            end
            
            out = self.makeCellIfNot (self.pinPosition);
            addcomma = ~(isempty (self.absolutePinOrientation) && isempty (self.initialTheta));
            str = self.addOutputLine (str, self.commaSepList ('position', out{:}), 2, addcomma, 'pin absolute position');
            
            if ~isempty (self.absolutePinOrientation)
                addcomma = ~isempty (self.initialTheta);
                out = self.makeCellIfNot (self.absolutePinOrientation);
                str = self.addOutputLine (str, self.commaSepList ('orientation', out{:}), 2, addcomma, 'pin absolute orientation');
            end
            
            if ~isempty (self.initialTheta)
                out = self.makeCellIfNot (self.initialTheta);
                str = self.addOutputLine (str, self.commaSepList ('initial theta', out{:}), 2, false, 'initial theta');
            end
            
            str = self.addOutputLine (str, ';', 1, false, 'end revolute pin');
        end
        
    end
    
end