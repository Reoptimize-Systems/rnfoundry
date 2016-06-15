classdef odederiv < handle
    % calculate numerical derivative of quantity for use in ode solver
    % routine
    
   properties (GetAccess = public, SetAccess = private)
       ulast;
       tlast;
       strict;
   end
   
   methods
       
       function self = odederiv (t0, u0, varargin)
           % constructor for odederiv class
           
           options.strict = false;
           
           options = parse_pv_pairs (options, varargin);
           
           self.strict = options.strict;
           
           self.ulast = u0;
           
           self.tlast = t0;
           
       end
       
       function dudt = derivative (self, t, u)
           % get the current value of the derivative
           %
           % call in ode function
           %
           % Syntax
           %
           % odederiv.derivative (t, u)
           %
           % Description
           %
           % derivative returns the numerical derivative of the supplied
           % quantity with respect to it's previous value such that
           %
           % du/dt = (u - prev u) / (t - prev t)
           %
           % derivative is vectorised so u may be a column of several
           % quantities for which the derivative is to be calculated, in
           % this case t is a scalar and u is a column vector of separate
           % quantities all at the point in time.
           %
           % In addition, derivative can take a row vector of n time values
           % and a matrix of values of u, e.g.
           %
           % t = [ t1, t2, t3 ]
           % u = [ u11, u12, u13;
           %       u21, u22, u23;
           %       u31, u32, u33 ]
           %
           % In this case the columns of u are the values of the quantities
           % at each time point in t. derivative will return a matrix
           %
           % dudt = [ d u11 / d t1, d u12 / d t2, d u13 / d t3;
           %          d u21 / d t1, d u22 / d t2, d u23 / d t3;
           %          d u31 / d t1, d u32 / d t2, d u33 / d t3 ]
           %
           % This facilitates vectorisation of ode solution functions
           %
           
           if any(t == self.tlast)
               dudt = zeros (size (u));
           else
               % rows of u are the separate components, the columns of u
               % are the new values of the components at each new value of
               % t. 
               dudt = bsxfun (@rdivide, bsxfun (@minus, u, self.ulast ), (t - self.tlast));
           end
           
       end
       
       function update (self, t, u)
           % update the previous values of t and u
           %
           % call in Outputfcn
           
           if self.strict
               if numel(t) > 1
                   error ('When updating, t must be a scalar');
               end
               if ~samesize(u, self.ulast)
                   error ('When updating, new value of u must be the same size as the previous value');
               end
           end
           
           self.ulast = u;
           self.tlast = t;
           
       end
       
       function status = outputfcn (self, t, u, flag, varargin)
           
           status = 0;
           
           if isempty (flag)
               
               self.update (t(end), u(:,end));
               
           elseif strcmp (flag, 'init')
               
           elseif strcmp (flag, 'done')
               
           end
           
       end
       
   end
    
end

% Test scraps
%
% od = odederiv (0, [1;2;3])
% 
% od.derivative (1, [4;4;4])
% 
% ans =
% 
%      3.0000e+000
%      2.0000e+000
%      1.0000e+000
% 
% od.derivative ([1, 2], [4, 5; 4, 5; 4, 5])
% 
% ans =
% 
%      3.0000e+000     2.0000e+000
%      2.0000e+000     1.5000e+000
%      1.0000e+000     1.0000e+000
% 
% od.derivative ([1, 2], [4, 4; 4, 4; 4, 4])
% 
% ans =
% 
%      3.0000e+000     1.5000e+000
%      2.0000e+000     1.0000e+000
%      1.0000e+000   500.0000e-003

