classdef sharedMemoryCommunicator < mbdyn.pre.externalFileCommunicator
    
    properties (GetAccess = public, SetAccess = private)
        create;
        sharedMemoryName;
    end
    
    methods
        
        function self = sharedMemoryCommunicator (varargin)
            
            options.SleepTime = [];
            options.Coupling = [];
            options.SendAfterPredict = 'yes';
            options.Create = [];
            options.SharedMemoryName = 'auto';
            
            options = parse_pv_pairs (options, varargin);
            
            if ~isempty (options.Create)
                mbdyn.pre.base.checkAllowedStringInputs (options.Create, {'yes', 'no'}, true, 'Create');
            end
            
            assert (ischar (options.SharedMemoryName), 'SharedMemoryName must be a character vector');
            
            if strcmpi (options.SharedMemoryName, 'auto')
                
                if isempty (options.Create) || strcmpi (options.Create, 'no'), ...
                    error ('If Create is ''no'' or empty, SharedMemoryName cannot be ''auto''');
                end
                
                if ~mbdyn.pre.base.isOctave ()
                    % make sure the random number seed is differrent in
                    % different matlab instances to avoid name clashes
                    rng('shuffle');
                end
                
                % make the path with the random name
                options.SharedMemoryName = sprintf ('mbdyn_shared_memory_%d', randi (100000) );
                
            end

            self = self@mbdyn.pre.externalFileCommunicator ( ...
                        'SleepTime', options.SleepTime, ...
                        'Coupling', options.Coupling, ...
                        'SendAfterPredict', options.SendAfterPredict );
                    
            self.type = 'shared memory';
            
            self.sharedMemoryName = options.SharedMemoryName;
            self.create = options.Create;
            
            self.commMethod = 'shared memory';
            
        end
        
        function str = generateMBDynInputString (self)
            % generates MBDyn input string for shared memory communicators
            % 
            % Syntax
            %  
            % str = generateMBDynInputString (shc)
            %  
            % Description
            %  
            % generateMBDynInputString is a method shared by all MBDyn
            % components and is called to generate a character vector used
            % to construct an MBDyn input file.
            %  
            % Input
            %  
            %  shc - mbdyn.pre.sharedMemoryCommunicator object
            %  
            % Output
            %  
            %  str - character vector for insertion into an MBDyn input
            %   file.
            %
            
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