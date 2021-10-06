function [sys, x0, str, ts] = machineidealpto_sfcn(t, x, u, flag, design, simoptions)
% machineidealpto_sfcn: level-1 matlab s-function for the simulation of a machine
% producing the ideal power take-off force
%
% Takes the standard sfunction inputs, t, x, u, and flag, and also required
% the presence in the workspace of some additional data passed in in
% design, and simoptions descibing the machine.  
%
%
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
            
            % Number of continuous variables to be solved, currently none 
            s.NumContStates = 0;
            
            % number of discrete states to be solved, zero as all variables
            % are continuous, and there are currently no derivatives to be 
            % solved in this case anyway
            s.NumDiscStates = 0;
            
            % Set the number of outputs, the EMF and the coil currents make
            % 6 outputs, plus another for the force making 7 in total
            s.NumOutputs = 12;
            
            % We take the following inputs
            %
            % Fpto - the desired power take-off force
            %
            % xBh - position of the buoy in heave 
            %
            % xBs - position in the buoy in surge 
            %
            % vBh - velocity of the buoy in heave
            %
            % vBs - the velocity of the buoy in surge
            %
            s.NumInputs = 5;
            
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
            
            % there are currently no derivatives to calculate
            sys = [ ];
            
        case 3 % output

            % Get the inputs
            
            Fpto = u(1);
            
            xBh = u(2);
            
            xBs = u(3);
            
            vBh = u(4);
            
            vBs = u(5);
            
            % Calculate the machine outputs
            [EMF, Icoils, limitIcoils, Force, Ploss, limitPloss] = machineidealfpto_slotless(design, simoptions, Fpto, xBh, xBs, vBh, vBs);
            
            % pass them back to simulink
            sys = [ EMF; Icoils; limitIcoils; Force; Ploss; limitPloss ]';
            
        case {2 4 9} 
            
            % 2:discrete
            % 4:calcTimeHit
            % 9:termination
            sys =[];
            
        otherwise
            
            error(['unhandled flag =',num2str(flag)]);
            
    end

end