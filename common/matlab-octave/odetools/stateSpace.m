classdef stateSpace < handle
    
    properties (GetAccess = public, SetAccess = private)
        
        A;
        B;
        C;
        D;
        x0;
        x;
        
        integrationReady;
        t0;
        t;
        
    end
    
    properties (GetAccess = private, SetAccess = private)
        
        % integration variables
        a;
        b;
        c;
        nx;
        xk; % Current state
        tk;
        s;  % Length of weights
        d;  % Matrix of derivatives
        updateIntegration;
        ufcn;
        
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
            SS.integrationReady = false;
            SS.updateIntegration = false;
            
        end
        
        function xdot = derivatives (SS, u)
            % get the state-space system derivatives
            
            xdot = SS.A * SS.x + SS.B * u;
            
        end
        
        function y = outputs (SS, u)
            % get the state-space system outputs
            
            y = SS.C * SS.x + SS.D * u;
            
        end
        
        function initIntegration (SS, t0, ufcn, a, b, c)
            
            SS.t0 = t0;
            
            if nargin < 5
                SS.a = [ 0, 0; 0.5, 0 ];
                SS.b = [ 0, 1 ];
                SS.c = [ 0, 0.5];
            else
                SS.a = a;
                SS.b = b;
                SS.c = c;
            end
            
            if nargin < 3
                SS.ufcn = @(t, x) 0;
            else
                SS.ufcn = ufcn;
            end
            
            SS.nx = numel(SS.x0);   % Number of states
            SS.xk = SS.x0(:);       % Current state

            SS.s = length(SS.b);       % Length of weights
            SS.d = zeros(SS.nx, SS.s);    % Matrix of derivatives
            
            SS.updateIntegration = true;
            SS.integrationReady = true;
            
        end
        
        function stepIntegrate (SS, dt)
            
            if SS.integrationReady
            
                % Current time
                SS.tk = SS.t + dt;

                % Calculate derivatives.
                SS.d(:, 1) = dt * SS.derivatives (SS.ufcn(SS.tk, SS.x));
                
                for z = 2:SS.s
                    
                    dxk = sum(bsxfun(@times, SS.a(z, 1:z-1), SS.d(:, 1:z-1)), 2);
                    
                    SS.d(:, z) = dt * SS.derivatives (SS.ufcn(SS.tk + SS.c(z) * dt, SS.xk + dxk));
                    
                end
                
            else
                error ('State-space integrator has not been initialised')
            end
            
        end
        
        function update (SS, x)
           % update the state space vector to a new value
           
           if SS.updateIntegration
               % update the current time
               SS.t = SS.tk;
               
               % Update the state.
               for z = 1:SS.s
                   SS.xk = SS.xk + SS.b(z) * SS.d(:, z);
               end
               
               SS.x = SS.xk;
           else
               SS.x = x;
           end
            
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