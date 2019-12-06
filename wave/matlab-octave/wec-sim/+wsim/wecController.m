classdef wecController < handle
% base class for WEC controllers
%
    
    properties
        
        wecSimObj; % the wsim.wecSim object running the simulation
        
    end
    
    properties (GetAccess=public, SetAccess=protected)
        
        usingSampleDelays;
        sampleDelayInfo;
        
    end
    
    
    methods
        
        function self = wecController (varargin)
            % wsim.wecController constructor
            %
            % Syntax
            %
            % ctrlobj = wsim.wecController ()
            % ctrlobj = wsim.wecController ('Parameter', Value)
            %
            % Description
            %
            % The wsim.wecController class is a base class for all
            % WEC controllers
            %
            % Input
            %
            % Addtional arguments may be supplied as parameter-value pairs.
            % The available options are:
            %
            %  'SampleDelays' - empty matrix or structure containing 
            %    specification of sample delayed logged variables to be
            %    used by the controller to generate the control outputs. If
            %    a structure it must contain at least the following fields:
            %
            %    SampleDelays : 
            %
            %    InitialValues : 
            %
            % Output
            %
            %  ctrlobj - wsim.wecController object
            %
            %
            %
            % See Also: 
            %
            
            options.SampleDelays = [];
            
            options = parse_pv_pairs (options, varargin);
            
            if isempty (options.SampleDelays)
                
                self.usingSampleDelays = false;
                self.sampleDelayInfo = [];
                
            else
                
                check.structHasAllFields (options.SampleDelays, {'SampleDelays', 'InitialValues'}, true, 'SampleDelays');
                
                self.sampleDelayInfo = options.SampleDelays;
                
                if all (self.sampleDelayInfo.SampleDelays == 0)
                    self.usingSampleDelays = false;
                else
                    self.usingSampleDelays = true;
                end
            end 
            
        end
        
        function info = loggingInfo (self)
            
            info = [];
            
        end
        
        function [value, info] = ptoControlOutput (self, id, pto_internal_variables)
            
            value = 0;
            info = struct ();
            
        end
        
        function start (self)
            
            
        end
        
        function advanceStep (self)
            
            
        end
        
        function finish (self)
            
            
        end
        
    end
    
    methods (Access = protected)
        
%         function out = timeDelayedLogVars (self, varnames, timevarnames, delays, initialvals)
%             
%             % copy intial vals for output as this will be correct size and
%             % hold the right default values which we will only overwrite
%             % once the current sim time 
%             out = initialvals;
%             
%             for varind = 1:numel (varnames)
%                 
%                 if self.logger.data.(timevarnames{ind})(self.logger.info.(timevarnames{ind}).LastLogIndex) - delays(varind) > self.wecSimObj.simInfo.TStart
%                     
%                 end
%                 
%             end
%             
%         end
        
        function out = sampleDelayedLogVars (self, varnames, sampledelays, initialvals)
            % returns logged variable delayed by a number of time steps
            %
            % Syntax
            %
            % out = sampleDelayedLogVars (ctrlobj, varnames, sampledelays, initialvals)
            %
            % Description
            %
            % sampleDelayedLogVars returns logged variables, but delayed by
            % a specified number of time steps. Until that number of steps
            % has passed an initial value for the variable is returned.
            %
            % Input
            %
            %  ctrlobj - wsim.wecController object
            %
            %  varnames - cell array of character vectors containing the
            %   names of the variables to be returned. These must be names
            %   of variables logged by the logger in the wsim.wecSim
            %    object this controller is added to.
            %
            %  sampledelays - vector of integer of the same length as
            %   varnames. Each element is the number of samples to delay
            %   the output for, for the corresponding variable in
            %   'variables' described above.
            %
            %  initialvals - vector of values of the same length as
            %   varnames. Each element is the value which will be reported
            %   for the corresponding variable in 'variables' described
            %   above until the specified number of samples has passed.
            %
            % Output
            %
            %  out - cell array of values for each of the variables
            %   requested in 'varnames'. Either the specified initial
            %   value, or the value logged the specified number of time
            %   steps ago.
            %
            %
            % See also:
            %
            
            % copy initial vals for output as this will be correct size and
            % hold the right default values which we will only overwrite
            % once the current sim time is greater than the sample delay
            out = initialvals;
            
            for varind = 1:numel (varnames)
                
                if self.wecSimObj.logger.info.(varnames{varind}).LastLogIndex >= sampledelays(varind)
                    
                    % get value at sampledelays - 1, as a value of 0 will
                    % give the last logged value, which is from the
                    % previous time step, i.e. a delay of one time step
                    % from the current time
                    out{varind} = getPrevLoggedVal(self.wecSimObj.logger, varnames{varind}, sampledelays(varind)-1);
                    
                end
                
            end
            
        end
        
    end
    
end