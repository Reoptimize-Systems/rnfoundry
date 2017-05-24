classdef stateSpace < handle
    
    properties (GetAccess = public, SetAccess = private)
        
        A;
        B;
        C;
        D;
        x0;
        x;
        
    end
    
    methods
        
        function SS = stateSpace (A, B, C, D, x0)
            
            SS.A = A;
            SS.B = B;
            SS.C = C;
            SS.D = D;
            
            SS.x0 = x0(:);
            
            % initialise the state-space system
            reset (SS);
            
        end
        
        function reset (SS)
            % reset the state-space system to its initial conditions
            
            SS.x = SS.x0;
        end
        
        function xdot = derivatives (SS, u)
            % get the state-space system derivatives
            
            xdot = SS.A * SS.x + SS.B * u;
            
        end
        
        function y = outputs (SS, u)
            % get the state-space system outputs
            
            y = SS.C * SS.x + SS.D * u;
            
        end
        
        function update (SS, x)
           % update the state space vector to a new value
           
           SS.x = x;
            
        end
        
        function status = outputfcn (SS, t, x, flag, varargin)
            % OutputFcn for calling in matlab's ode solvers after a time
            % step is completed
            status = 0;
            
            if isempty (flag)
                
                SS.update (x(:,end));
                
            elseif strcmp (flag, 'init')
                
                % input is tspan and y0
                SS.update (x(:));
                
            elseif strcmp (flag, 'done')
                
                
            end
            
        end
        
    end
    
end