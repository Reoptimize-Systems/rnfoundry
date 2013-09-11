classdef odeinternals < handle
    % class for extracting internally calculated values in an ode solution
    % function
    %
    % 
    
   properties (SetAccess = private)
       t_intermediate = [];
       vals_intermediate = [];
       y_intermediate = [];
   end
   
   properties (SetAccess = private, GetAccess = public)
      t_output = [];
      vals_output = [];
   end
   
   methods
       
       % constructor
       function obj = oderesults()
           this.reset();
       end
       
       function reset(this)
           % reset the odeinternals object by clearing all the history
           % values
           
           this.clearoutput();
           this.clearintermediate()
       end
       
       function clearoutput(this)
           
           this.t_output = [];
           this.vals_output = [];
           
       end
       
       function clearintermediate(this)
           
           this.t_intermediate = [];
           this.y_intermediate = [];
           this.vals_intermediate = [];
           
       end
       
       function addintermediatestep(this, t, y, vals)
           
           this.t_intermediate = [this.t_intermediate, t];
           this.y_intermediate = [this.y_intermediate, y(:)];
           this.vals_intermediate = [this.vals_intermediate, vals(:)];
           
       end
       
       function addoutputstep(this, t, y)
           % adds an output step to the history
           
           ind = find(this.t_intermediate == t & all(bsxfun(@eq, this.y_intermediate, y(:))));
           
           if isempty(ind)
               error('Output time step could not be located in the list of ');
           end
           
           this.t_output = [this.t_output, t];
           this.vals_output = [this.vals_output, this.vals_intermediate(:,ind)];
           
           % now reset the temporary history from the previous step
           this.clearintermediate();
           
       end
       
   end
   
end