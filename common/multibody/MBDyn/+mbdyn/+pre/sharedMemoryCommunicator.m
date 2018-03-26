classdef sharedMemoryCommunicator < mbdyn.pre.externalFileCommunicator
    
    properties (GetAccess = public, SetAccess = private)
        create;
        sharedMemoryName;
    end
    
    methods
        
        function self = sharedMemoryCommunicator (name, create, varargin)
            
            options.SleepTime = [];
            options.Coupling = [];
            options.SendAfterPredict = 'yes';
            
            options = parse_pv_pairs (options, varargin);
            
            
            self = self@mbdyn.pre.externalFileCommunicator ( ...
                        'SleepTime', options.SleepTime, ...
                        'Coupling', options.Coupling, ...
                        'SendAfterPredict', options.SendAfterPredict );
                    
            self.type = 'shared memory';
            
            if ischar (name)
                self.sharedMemoryName = name;
            else
                error ('Shared memory region name in ''name'' must be a string');
            end

            self.checkAllowedStringInputs (create, {'yes', 'no'}, true, 'Create');
            
            self.create = create;
            
            self.commMethod = 'shared memory';
            
        end
        
        function str = generateMBDynInputString (self)
            
            str = generateMBDynInputString@mbdyn.pre.externalFileCommunicator(self);
            
%             if ~isempty (self.create)
                
                str = self.addOutputLine (str, self.commaSepList ('create', self.create), 1, true, 'will mbdyn create shared memory?');
                
%             end

            addcomma = any ( ~[ isempty(self.sleepTime);
                                isempty(self.coupling);
                                isempty(self.sendAfterPredict) ] );
                            
            str = self.addOutputLine (str, self.commaSepList ('name', ['"', self.sharedMemoryName, '"']), 1, addcomma, 'shared memory region name');

            if ~isempty (self.sleepTime)
                addcomma = any ( ~[ isempty(self.coupling);
                                    isempty(self.sendAfterPredict) ] );
                                
                str = self.addOutputLine (str, self.commaSepList ('sleep time', self.sleepTime), 1, addcomma);
            end
            
            if ~isempty (self.coupling)
                addcomma = ~isempty(self.sendAfterPredict);
                                
                str = self.addOutputLine (str, self.commaSepList ('coupling', self.coupling), 1, addcomma);
            end
            
            if ~isempty (self.sendAfterPredict)
                str = self.addOutputLine (str, self.commaSepList ('send after predict', self.sendAfterPredict), 1, false);
            end
            
        end
        
        function comminfo = commInfo (self)
            % gets communication info for the shared memory communicator.
            %
            %
            
            comminfo.commMethod = self.commMethod;
            
            if strcmp (self.commMethod, 'shared memory')
                comminfo.shared_mem_name = self.sharedMemoryName;
            else
                error ('unrecognised communication type');
            end

        end
        
    end
    
end