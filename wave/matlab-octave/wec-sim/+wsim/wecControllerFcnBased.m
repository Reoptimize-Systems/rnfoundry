classdef wecControllerFcnBased < wsim.wecController
% base class for CEORl WEC controllers
    
    properties
        
        numControlParams;
        paramsMinMax;
        outputFcns;
        ptoVariableNames;
        ptoVarSampleDelays
        wecSimVariableNames;
        wecSimVarSampleDelays;
        
    end
    
    methods
        
       function self = wecControllerFcnBased (fcn, ptovarnames, wecsimvarnames, varargin)
            % ceorl.wecController constructor
            %
            % Syntax
            %
            % ctrlobj = ceorl.wecController ()
            % ctrlobj = ceorl.wecController ('Parameter', Value)
            %
            % Description
            %
            % The ceorl.wecController class is a base class for all ceorl
            % WEC controllers
            %
            % Input
            %
            % Addtional arguments may be supplied as parameter-value pairs.
            % The available options are:
            %
            %  'SampleDelays' - empty matrix or structure containing 
            %    specification of sample delayed logged variables to be
            %    used by the controller to generate the control outputs
            %
            % Output
            %
            %  ctrlobj - ceorl.wecController object
            %
            %
            %
            % See Also: 
            %
            
            options.PTOVarSampleDelays = struct ('SampleDelays', zeros(1, numel(ptovarnames)), ...
                                                 'InitialValues', { cell(1, numel(ptovarnames)) } );
                                             
            options.WECSimVarSampleDelays = struct ('SampleDelays', zeros(1, numel(wecsimvarnames)), ...
                                                    'InitialValues', { cell(1, numel(wecsimvarnames)) } );
            
            options = parse_pv_pairs (options, varargin);
            
            if iscell (fcn)
                for ind = 1:numel(fcn)
                    assert (isa ('function_handle', fcn{ind}),...
                            'fcn{%d} is not a function handle.', ...
                            ind);
                end
            elseif isa (fcn, 'function_handle')
                fcn = {fcn};
            else
                error ('fcn must be a function returning a force or torque, or a cell array of such functions');
            end
            
            sample_delays = struct ();
            
            sample_delays.SampleDelays = [ options.PTOVarSampleDelays.SampleDelays, ...
                                           options.WECSimVarSampleDelays.SampleDelays ];
                                       
            sample_delays.InitialValues = [ options.PTOVarSampleDelays.InitialValues, ...
                                            options.WECSimVarSampleDelays.InitialValues ];
                
            self = self@wsim.wecController ('SampleDelays', sample_delays);
            
            self.ptoVariableNames = ptovarnames;
            self.ptoVarSampleDelays = options.PTOVarSampleDelays;
            
            self.wecSimVariableNames = wecsimvarnames;
            self.wecSimVarSampleDelays = options.WECSimVarSampleDelays;
            
            self.outputFcns = fcn;
            
       end
        
       function value = ptoControlOutput (self, pto_id, pto_vars)
            % returns the desired force based on the current damping value and pto velocity
            %
            % Syntax
            %
            % value = ptoControlOutput (ctrlobj, pto_id, pto_vars)
            %
            % Input
            %
            %  ctrlobj -
            %   ceorl.point_absorber.wecControllerSinglePTOActiveDamping
            %   object
            %
            %  pto_id - unused (but provided to match the standard syntax
            %    of the ptoControlOutput method)
            %
            %  pto_vars - structure containing the field 'RelativeVelocity'
            %    which contains the relative velocity of the two parts of
            %    the linear PTO.
            %
            %
            % Output
            %
            %  value - the force to be applued based on the damping
            %   coefficient and the relative velocity of the PTO components
            %
            %
            % See also:
            %
            
            fcn_input_vars = {};
            
            for ind = 1:numel(self.ptoVariableNames)
                
                if self.usingSampleDelays && self.ptoVarSampleDelays.SampleDelays(ind) ~= 0
                    fcn_input_vars(end+1) = sampleDelayedLogVars ( self, ...
                                                          {['PTO_', int2str(pto_id), '_', self.ptoVariableNames{ind}]}, ...
                                                          self.ptoVarSampleDelays.SampleDelays(ind), ...
                                                          self.ptoVarSampleDelays.InitialValues(ind) );
                else
                    fcn_input_vars(end+1) = pto_vars.(self.ptoVariableNames{ind});
                end
            
            end
            
            for ind = 1:numel(self.wecSimVariableNames)
                
                if self.usingSampleDelays && self.wecSimVarSampleDelays.SampleDelays(ind) ~= 0
                    fcn_input_vars(end+1) = sampleDelayedLogVars ( self, ...
                                                          self.wecSimVariableNames(ind), ...
                                                          self.wecSimVarSampleDelays.SampleDelays(ind), ...
                                                          self.wecSimVarSampleDelays.InitialValues(ind) );
                else
                    fcn_input_vars(end+1) = self.wecSimObj.(['last', self.wecSimVariableNames{ind}]);
                end
            
            end
            
            value = feval (self.outputFcns{pto_id}, fcn_input_vars{:});
            
        end
        
        function setControlParameters (self, new_params)
            
            
            
        end
        
    end
    
end