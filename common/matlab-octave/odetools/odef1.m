function [varargout] = odef1(odefun,tspan,y0,varargin)
%ODEF1  Solve differential equations with a non-adaptive method of order 1.
%   Y = ODE1(ODEFUN,TSPAN,Y0) with TSPAN = [T1, T2, T3, ... TN] integrates 
%   the system of differential equations y' = f(t,y) by stepping from T0 to 
%   T1 to TN. Function ODEFUN(T,Y) must return f(t,y) in a column vector.
%   The vector Y0 is the initial conditions at T0. Each row in the solution 
%   array Y corresponds to a time specified in TSPAN.
%
%   Y = ODE1(ODEFUN,TSPAN,Y0,P1,P2...) passes the additional parameters 
%   P1,P2... to the derivative function as ODEFUN(T,Y,P1,P2...). 
%
%   This is a non-adaptive solver. The step sequence is determined by TSPAN.
%   The solver implements the forward Euler method of order 1.   
%
%   Example 
%         tspan = 0:0.1:20;
%         y = ode1(@vdp1,tspan,[2 0]);  
%         plot(tspan,y(:,1));
%     solves the system y' = vdp1(t,y) with a constant step size of 0.1, 
%     and plots the first component of the solution.   
%

    if ~isvector(tspan)
      error('tspan should be a vector of integration time steps.');
    end

    if ~isvector(y0)
      error('y0 should be a vector of initial conditions.');
    end

    % get the time step sizes
    h = diff(tspan);
    
    if any(sign(h(1))*h <= 0)
        error('Entries of tspan are not in order.') 
    end  

    if (nargin >= 4)
    
        if (~isstruct (varargin{1}))
            % there is no options structure
            odeoptions = odeset ();
            odefcnargs = varargin;
        elseif (numel (varargin) > 1)
            % varargin{1} should be the options structure (should really do some 
            % checking here)
            odeoptions = varargin{1};
            odefcnargs = varargin(2:end);
        else
            odeoptions = varargin{1};
            odefcnargs = {};
        end
        
    else
        vodeoptions = odeset (); 
        odefcnargs = {};
        
    end
    
    if ~isempty (odeoptions.OutputFcn)
        % check the output function is actually a function
        if ischar (odeoptions.OutputFcn)
            odeoptions.OutputFcn = str2func (odeoptions.OutputFcn);
        end
        
        if ~is_function_handle (odeoptions.OutputFcn)
            haveoutputfcn = true;
        else
            error ('''OutputFcn'' is not a function.');
        end
    else
        haveoutputfcn = false;
    end
    
    % TODO: add event functions
    haveeventfcn = false;
    foundevents = {[], [], [], []};
    
    y0 = y0(:);   % Make y0 a column vector.
    
    try
      f0 = feval(odefun,tspan(1),y0,odefcnargs{:});
    catch
      msg = ['Unable to evaluate the ODEFUN at t0,y0. ',lasterr];
      error(msg);  
    end  

    if ~isequal(size(y0),size(f0))
      error('Inconsistent sizes of Y0 and f(t0,y0).');
    end  
    
    if haveoutputfcn
        feval (odeoptions.OutputFcn, tspan, y0, 'init', vfunarguments{:});  
    end

    neq = numel(y0);
    N = numel(tspan);
    T = tspan(:);
    Y = zeros(neq,N);

    Y(:,1) = y0;
    for ind = 1:N-1 
    
        Y(:,ind+1) = Y(:,ind) + h(ind)*feval(odefun,T(ind),Y(:,ind),odefcnargs{:});
        
        % run output function if supplied
        if haveoutputfcn
            feval (odeoptions.OutputFcn, T(ind+1), Y(:,ind+1), [], odefcnargs{:});  
        end
        
        % TODO: check event function
        
    end
    
    if haveoutputfcn
        feval (odeoptions.OutputFcn, [], [], 'done', vfunarguments{:});  
    end
    
    % now supply the outputs
    if (nargout == 1)                 % Sort output variables, depends on nargout
      varargout{1}.x = T.';           % Time stamps are saved in field x
      varargout{1}.y = Y;             % Results are saved in field y
      varargout{1}.solver = 'odef1';  % Solver name is saved in field solver
      if (haveeventfcn) 
        varargout{1}.ie = foundevents{2};  % Index info which event occured
        varargout{1}.xe = foundevents{3};  % Time info when an event occured
        varargout{1}.ye = foundevents{4};  % Results when an event occured
      end
      if (strcmpi (odeoptions.Stats, 'on'))
        varargout{1}.stats = struct;
        varargout{1}.stats.nsteps   = N;
        varargout{1}.stats.nfailed  = 0;
        varargout{1}.stats.nfevals  = N;
      end
    elseif (nargout == 2)
      varargout{1} = T;     % Time stamps are first output argument
      varargout{2} = Y.';   % Results are second output argument
    elseif (nargout == 5)
      varargout{1} = T;     % Same as (nargout == 2)
      varargout{2} = Y.';   % Same as (nargout == 2)
      varargout{3} = [];
      varargout{4} = [];
      varargout{5} = [];
      if (haveeventfcn) 
        varargout{3} = foundevents{3};  % Event occurance times
        varargout{4} = foundevents{4};  % Values when events occured
        varargout{5} = foundevents{2};  % Event index
      end
    end
  
end
