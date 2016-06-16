function varargout = prescribedmotodeforcefcn_linear_mvgarm(t, x, design, simoptions)
% solves the right had side of the differential equations for a linear
% generator system with a moving armature and translator
%
% Syntax:
%
% Input
%
% x is a vector of values of the state of the sytem at time t. The members
% of x are the following:
%
%
% design is a standard machine design structure popluated with all the
% necessary information to perform the simulation. 
%
% 'simoptions'  is a structure determining the buoy and sea parameters used
% in the simulation, and also any penalties to be used in the scoring of
% the design. Essential fields are:

%

    if nargin == 0
        
        % get the output names
        varargout{1} = {'ydot', 'dpsidxR', 'EMF', 'xT', 'vT', 'Fpto', ...
                        'FfA', 'FaddA', 'FaddT', 'FaddEBD', 'RPhase'};
        return;
    end

    % Change the x members into more useful variables names, MATLAB will
    % optimise away any memory penalty associated with this I think
    xA  = x(design.Phases+1);
    vA  = x(design.Phases+2);
    
    Iphases = x(1:design.Phases);
    
    Icoils = Iphases ./ design.Branches;

    % Initialize dx with zeros
    dx = zeros(size(x));
    
    % get the velocity and position at the current time
    [xTtemp, vTtemp] = prescribedmotvelpos(t, simoptions);
    
%     % Interpolate the data set (times,xT) at current time
%     xTtemp = interp1(simoptions.drivetimes, simoptions.xT, t); 
%     
%     % Interpolate the data set (times,vT) at current time
%     vTtemp = interp1(simoptions.drivetimes, simoptions.vT, t); 

    simoptions.BuoySim.tether_length = 1000;

    % determine the machine outputs
    [dpsidxR, EMF, FA, FfeaVec, xT, vT, design] = machineodesim_linear_mvgarm(design, simoptions, Icoils, xA, vA, xTtemp, vTtemp, 0, 0);

    % find the derivative of the coil current (solving the differential
    % equation describing the simple output circuit)
    dx(1:design.Phases,1) = circuitode_linear(Iphases, EMF, design);
    
    % determine the forces due to the magnets and electrical forces at
    % the relative position xR with the current values of J. Forces are
    % fitted to a 1m stack length, so we adjust for this by multiplying by
    % ls, the actual stack length in m
    %Ffea = sum(intbpolyshearforce_AC(design, J, pos)) .* design.Poles(1);

    % Calculate the drag forces on the translator
    % Fdrag = sign(vT) .* 0.5 .* realpow(vT,2) .* simoptions.BuoySim.BuoyParameters.rho .* design.Cd .* design.DragArea;
     
    [FaddT, FaddA, ForceBD] = feval(simoptions.ODESim.ForceFcn, design, simoptions, xT, vT, xA, vA, xTtemp, 0, vTtemp, 0, simoptions.ODESim.ForceFcnArgs{:});
    
    % Find the frictional force that would cause the armature to move at  
    % the same speed as the translator 
    
%     if abs(vA) > 0.00001 %2 * max(simoptions.abstols(5))
%         
%         % if the armature is moving, the frictional force is simply
%         % whatever the value of mu N is.
%         
%         FfA = Fa(2) * -sign(vA);
%         dx(design.Phases+1,1) = vA;
%         aA = (Fpto + (sum(Fa) - Fa(2)) + FfA - design.weightA) ./ design.massA;
%         
%     elseif abs(Fpto + sum(Fa) - Fa(2)) < Fa(2)
%         
%         % FfA = Ffriction * (2/pi) * atan(epsilon * -vR) / (1 + dconst * abs(vR));
%         
%         % if the armature is not moving, and the net forces are less than
%         % the frictional force, the acceleration will be zero, and the
%         % frictional force equal to the net of the other forces
%         FfA = -(sum(Fa) - Fa(2));
%         dx(design.Phases+1,1) = vA;
%         aA = 0;
% 
%     else
%         
%         % if the armature is not moving, but the frictional force is less
%         % than the net forces, there will be an acceleration
%         FfA = Fa(2) * -sign((sum(Fa) - Fa(2)));
%         dx(design.Phases+1,1) = vA;
%         aA = (Fpto + (sum(Fa) - Fa(2)) + FfA - design.weightA) ./ design.massA;
%         
%     end    
%     
%     Fa(2) = FfA;
    
    % determine the acceleration of the armature
%     aA = (Fpto + sum(Fa) - Fa(2) - design.weightA) ./ design.massA;
%     dx(design.Phases+1,1) = vA;

    [FfA, aA] = odefriction(vA, design.massA, design.AngleFromHorizontal, design.mu_fA, [FaddA, FA]);
    
    % determine the acceleration of the armature
%     aA = (Fpto + sum(Fa) - design.weightA) ./ design.massA;
    
    dx(design.Phases+1,1) = vA;
    dx(design.Phases+2,1) = aA;
    
    % ************************************************************************

    % Now assign the outputs
    varargout{1} = dx;

    % To record the forces
    if nargout > 1
        
        % rate of change of flux linkage with displacement 
        varargout{2} = dpsidxR;
        % per-coil induced voltage 
        varargout{3} = EMF;
        % vertical displacement of the translator  
        varargout{4} = xT;
        % vertical velocity of the translator
        varargout{5} = vT;
        % total electromagnetic forces acting on effector
        varargout{6} = -FA;
        % frictional forces on the armature
        varargout{7} = FfA;
        % ouput the additional armature forces
        varargout{8} = FaddA;
        % output the additional forces transmitted to the translator
        varargout{9} = FaddT;
        % get the breakdown of forces
        varargout{10} = ForceBD;
        % output the Phases resistance
        varargout{11} = diag(design.RPhase);
        
    end

end

