function [sys, x0, str, ts] = tubular_sfcn(t, x, u, flag, Mdesign, simoptions)

    switch flag
        
        case 0 % initialize
            
            % str needs to be set to null. This is reserved for use in
            % future versions of Simulink.
            str = [] ;
            
            % ts is an array of two columns to specify sampling time and
            % time offsets for discrete systems. As this function is a
            % continuous system, we set this to [0 0] during initiation
            % phase.
            ts = [0 0];
            
            % The initial conditions of the variables
            x0 = [];
            
            s = simsizes;
            
            % Number of continuous variables to be solved
            s.NumContStates = 0;
            
            % number of discrete states to be solved, zero as all variables
            % are continuous
            s.NumDiscStates = 0;
            
            % Set the number of outputs
            s.NumOutputs = 4;
            
            % We will require the current as an input to determine the
            % correct forces
            s.NumInputs = 3;
            
            % current output will depend on voltage at each time step so we
            % require direct algebraic feed through of input to output.
            s.DirFeedthrough = 1 ;
            
            % number of sample times. ( for continuous process, this is set
            % to 1)
            s.NumSampleTimes = 1 ;
            
            % Using the command simsizes again with the structure variable
            % as the argument actually translates the values in the
            % structure, s, into a row vector which gets sent to Simulink
            % via sys
            sys = simsizes(s);
            
        case 1 % derivatives
            
            % No derivatives to be solved at the minute
            sys = [];
            
        case 3 % output

            % Set the coil current
            Icoil = u;
            
            % extract the outputs
            [EMF, Fy] = slink_TM(y, v, Icoil, design, ...
                            simoptions.dpsidxfun, simoptions.Fyfun);

            % pass them back to simulink
            sys = [EMF(1), EMF(2), EMF(3), Fy];
            
        case {2 4 9} 
            
            % 2:discrete
            % 4:calcTimeHit
            % 9:termination
            sys =[];
            
        otherwise
            
            error(['unhandled flag =',num2str(flag)]);
            
    end

end