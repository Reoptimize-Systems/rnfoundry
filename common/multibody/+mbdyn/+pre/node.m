classdef node < mbdyn.pre.base
    
    properties (GetAccess = public, SetAccess = protected)
        
        humanReadableLabel;
        
%         uniqueTextLabel;
        
        scale;

        output;
        
    end

    methods
        
        function self = node (varargin)
            
            options.HumanReadableLabel = '';
            options.Output = true;
            options.Scale = [];
%             options.UniqueTextLabelPrefix = 'node';
            
            options = parse_pv_pairs (options, varargin);
            
            if ~ischar (options.HumanReadableLabel)
                error ('''HumanReadableLabel'' must be a char array');
            end
            
            if ~islogical (options.Output)
               error ('''Output'' must be a boolean true/false value');
            end
            
            self.humanReadableLabel = options.HumanReadableLabel;
            self.output = options.Output;
            self.scale = options.Scale;
            
        end
        
        
    end

end
