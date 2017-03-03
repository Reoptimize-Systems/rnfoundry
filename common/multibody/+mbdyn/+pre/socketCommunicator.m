classdef socketCommunicator < mbdyn.pre.externalFileCommunicator
    
    properties
        create;
        port;
        path;
        host;
    end
    
    methods
        
        function self = socketCommunicator (varargin)
            
            options.SleepTime = [];
            options.Coupling = [];
            options.SendAfterPredict = 'yes';
            options.Create = [];
            options.Host = [];
            options.Port = [];
            options.Path = [];
            
            options = parse_pv_pairs (options, varargin);
            
            
            self = self@mbdyn.pre.externalFileCommunicator ( ...
                        'SleepTime', options.SleepTime, ...
                        'Coupling', options.Coupling, ...
                        'SendAfterPredict', options.SendAfterPredict );
                    
            self.type = 'socket';

            if ~isempty (options.Create)
                self.checkAllowedStringInputs (options.Create, {'yes', 'no'}, true, 'Create');
            end
            
            if ~isempty (options.Path) && ~isempty (options.Port)
                error ('You cannot specify both path and port option for the socket');
            end
            
            if isempty (options.Path) && isempty (options.Port)
                error ('You must specify either Path or Port option for the socket');
            end
            
            self.create = options.Create;
            self.path = options.Path;
            self.port = options.Port;
            self.host = options.Host;
            
        end
        
        function str = generateOutputString (self)
            
            str = generateOutputString@mbdyn.pre.externalFileCommunicator(self);
            
            if ~isempty (self.create)
                
                str = self.addOutputLine (str, self.commaSepList ('create', self.create), 1, true, 'reference node');
                
%                 if strcmp (self.Create, 'no')
%                     
%                 end
                
            end
            
            addcomma = any ( ~[ isempty(self.sleepTime);
                                isempty(self.coupling);
                                isempty(self.sendAfterPredict) ] );
                            
            if isempty (self.path)
                % use port
                strline = self.commaSepList ('port', self.port);
                if ~isempty (self.host)
                    strline = [strline, self.commaSepList('host', self.host)];
                end
                str = self.addOutputLine (str, strline, 1, addcomma);
            else
                % use path
                str = self.addOutputLine (str, self.commaSepList ('path', ['"', self.path, '"']), 1, addcomma);
            end
            
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
        
    end
    
end