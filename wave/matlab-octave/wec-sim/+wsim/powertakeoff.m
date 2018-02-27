classdef powertakeoff < handle
   
    properties 
        id;
        
    end
    
    properties (GetAccess=public, SetAccess=private)
        uniqueLoggingNames;
        referenceNode;
        otherNode;
    end
    
    properties (GetAccess=protected, SetAccess=private)
        logger;
        loggerReady;
    end
    
    methods (Abstract)
        
        forceAndMoment (self);
%         forceSize (self);
%         initDataStructure (self);
        loggingSetup(self, logger);
        logData (self)
        
    end
    
    methods
        
        function self = powertakeoff (reference_node, other_node)
            
            if ~isa (reference_node, 'mbdyn.pre.structuralNode')
                error ('reference_node must be an mbdyn.pre.structuralNode')
            end
            
            if ~isa (other_node, 'mbdyn.pre.structuralNode') 
                error ('other_node must be an mbdyn.pre.structuralNode')
            end
            
            self.referenceNode = reference_node;
            self.otherNode = other_node;
            
            self.loggerReady = false;
%             self.id = 0;
            self.logger = [];
            
        end
        
        function set.id (self, newid)
            % set new integer id for power take off
            
            check.isPositiveScalarInteger (newid, true, 'id');
            
            self.id = newid;
        end
        
    end
    
    methods (Access=protected)
        
        function initLogging (self, info, logger)
            
            for ind = 1:numel (info.AvailableNames)
                info.AvailableNames{ind} = sprintf ('PTO_%d_%s', self.id, info.AvailableNames{ind});
            end
            
            self.uniqueLoggingNames = info.AvailableNames;
            
            if ~isempty (logger)
                
                assert (isa (logger, 'wsim.logger'), ...
                    'logger must be a wsim.logger object (or empty)');
            
                initLoggerObj (self, info, logger)
                
            end

            
        end
        
        function initLoggerObj (self, info, logger)
            
            for ind = 1:info.NAvailable
                
                % putting IndepVar in will automatically allocate enough
                % space for variable, same length as indepvar, if there is
                % no Indep var, preallocated length will be 1
                logger.addVariable ( info.AvailableNames{ind}, ...
                                     info.Sizes{ind}, ...
                                     'Desc', info.Descriptions{ind}, ...
                                     'Indep', info.IndepVars{ind} );
                
            end
            
            self.logger = logger;
            
            self.loggerReady = true;
            
        end
        
    end
    
end