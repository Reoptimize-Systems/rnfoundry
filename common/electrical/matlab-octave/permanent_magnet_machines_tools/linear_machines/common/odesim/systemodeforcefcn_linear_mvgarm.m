function varargout = systemodeforcefcn_linear_mvgarm(t, x, design, simoptions)
% systemodeforcefcn_linear_mvgarm: solves the right had side of the
% differential equations for the combined snapper and buoy system
%
% Syntax:
%
% [dx] = systemode_linear(t, x, design, simoptions)
% 
% [..., FBDh] = ] = systemode_linear(t, x, design, simoptions)
%     
% [..., bouyancy_force] = systemode_linear(t, x, design, simoptions)
% 
% [..., excitation_force_heave] = systemode_linear(t, x, design, simoptions)
% 
% [..., excitation_force_surge] = systemode_linear(t, x, design, simoptions)
% 
% [..., Hrad_force] = systemode_linear(t, x, design, simoptions)
% 
% [..., Srad_force] = systemode_linear(t, x, design, simoptions)
% 
% [..., dpsidxR] = systemode_linear(t, x, design, simoptions)
% 
% [..., EMF] = systemode_linear(t, x, design, simoptions)
% 
% [..., Ffea] = systemode_linear(t, x, design, simoptions)
% 
% [..., Fs] = systemode_linear(t, x, design, simoptions)
% 
% [..., Ffea_surge] = systemode_linear(t, x, design, simoptions)
%
% Input
%
% x is a vector of values of the state of the sytem at time t. The members
% of x are the following:
%
% x(1) = Buoy position - heave xBh
% x(2) = Bouy velocity - heave vBh
% x(3) = Bouy position - surge xBs
% x(4) = Buoy velocity - surge vBs
% x(5) = I, coil current in phase 1
% x(6) = I, coil current in phase 2
% x(7) = I, coil current in phase 3
% x(8) = I(1) (heave)
% x(9) -- x(N) = I(2) - I(N-6) (heave)
% x(N+1) -- x(N+M) = I(1) - I(M) (surge)
%
% design is a standard machine design structure popluated with all the
% necessary information to perform the simulation. 
%
% 'simoptions'  is a structure determining the buoy and sea parameters used
% in the simulation, and also any penalties to be used in the scoring of
% the design. Essential fields are:
%
%	- BuoyParameters - 
%
%   a structure containing relevent detailsabout the characteristics of the
%   buoy such as radius, draft etc. See defaultbuoyparameters for details. 
%   It should contain the following members:
%
%   Halpha - alpha coeficients in heave
%   Hbeta - beta coefficients in heave
%   Salpha - alpha coeficients in surge
%   Sbeta - beta coeficients in surge
%   drag_coefficient - coefficient of drag on the buoy
%   a - radius of the buoy
%   g - local gravitational constant (in case we deploy on Mars or Moon in 
%       future)
%   rho - density of seawater (see above)
%
%	- SeaParameters - 
%
%   optional structure containing relevent details of the sea to be used in
%   the simulation. If not supplied the structure returned by the function
%   defaultseaparameters will be used. See defaultseaparameters for further
%   details. It should contain the following members:
%
%   amp - vector of wave amplitudes for each wave sinusoid being applied
%   phase - vector of phase values for each sinusoid being applied
%   sigma - vector of angular frequencies for each sinusoid
%   wave_number -  vector of wave numbers for the frequencies
%   tether_length - the length of the tether from the buoy to the hawse
%                   hole
%
% Output:
%
%   dx - solution to the rhs of differential equations describing the
%   system at the current time step
%
%   EMF - per-coil induced voltage at the current time step
%
%   Ffea - total electromagnetic forces acting on armature at the current
%   time step
%
%   xT - vertical displacement of the translator at the current time step
%
%   vT - vertical velocity of the translator at the current time step
%
%   Ffeah - heave component of generator forces on buoy at the current time
%   step
%
%   Ffeas - surge component of generator forces on buoy at the current time
%   step
%
%   FBDh - buoy drag forces in heave at the current time step
%
%   FBDs - buoy drag forces in surge at the current time step
%
%   bouyancy_force - simple buoyancy force at the current time step
%
%   excitation_force_heave - excitation force in heave at the current time
%   step
%
%   excitation_force_surge - excitation force in surge at the current time
%   step
%
%   radiation_force_heave - radiation force in heave at the current time
%   step
%
%   radiation_force_surge - radiation force in surge at the current time
%   step
%
%   dlambdaVdxR -rate of change of flux linkage with relative displacement
%   at the current time step
%

    if nargin == 0

        % get the output names
        varargout{1} = {'ydot', 'EMF', 'Fpto', 'xT', 'vT', ...
                        'Ffea_heave', 'Ffea_surge', 'FBDh', 'FBDs', ...
                        'buoyancy_force', 'excitation_force_heave', ...
                        'excitation_force_surge', 'radiation_force_heave', ...
                        'radiation_force_surge', 'dpsidxR', 'FfA', ...
                        'FaddA', 'FaddT', 'FaddEBD', 'wave_height', 'RPhase'};
        return;
    end
    
%     if (t > 17 && any(abs(x(1:6)) > 10)) || t >  17.496919605
%         keyboard
%     end
    
    % Change the x members into more useful variables names, MATLAB will
    % optimise away any memory penalty associated with this I think
    xBh = x(1);
    vBh = x(2);
    xBs = x(3);
    vBs = x(4);
    xA  = x(5);
    vA  = x(6);
    Iphases = x(7:6+design.phases);
    
    Icoils = Iphases ./ design.Branches;

    % Initialize dx with zeros
    dx = zeros(size(x));

    % The differentials of the positions are the velocities
    dx(1,1) = vBh;
    dx(3,1) = vBs;
    dx(5,1) = vA;

    % determine the machine outputs
    [dpsidxR, EMF, FA, FfeaVec, xT, vT, design] = ...
        machineodesim_linear_mvgarm(design, simoptions, Icoils, xA, vA, xBh, vBh, xBs, vBs);

    % find the derivative of the coil current (solving the differential
    % equation describing the simple output circuit)
    dx(7:6+design.phases,1) = circuitode_linear(Iphases, EMF, design);

%     if any([xT, vT, xA, vA, xBh, xBs, vBh, vBs] > 30)
%         keyboard
%     end
    
    % Calculate the drag forces on the translator
    % Fdrag = sign(vT) .* 0.5 .* realpow(vT,2) .* simoptions.BuoyParameters.rho .* design.Cd .* design.DragArea;
    [FaddT, FaddA, FaddEBD] = feval(simoptions.forcefcn, design, simoptions, xT, vT, xA, vA, xBh, vBh, xBs, vBs, simoptions.forcefcnargs{:});
    
%     if abs(vA) > 0.00001 %2 * max(simoptions.abstols(5))
%         
%         % if the armature is moving, the frictional force is simply
%         % whatever the value of mu N is.
%         
%         FfA = Fa(2) * -sign(vA);
%         dx(design.phases+1,1) = vA;
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
%         dx(design.phases+1,1) = vA;
%         aA = 0;
% 
%     else
%         
%         % if the armature is not moving, but the frictional force is less
%         % than the net forces, there will be an acceleration
%         FfA = Fa(2) * -sign((sum(Fa) - Fa(2)));
%         dx(design.phases+1,1) = vA;
%         aA = (Fpto + (sum(Fa) - Fa(2)) + FfA - design.weightA) ./ design.massA;
%         
%     end    
    
    [FfA, aA] = odefriction(vA, design.massA, design.AngleFromHorizontal, design.mu_fA, [FaddA, FA]);
    
%     Fa(2) = FfA;
    % determine the acceleration of the armature
%     aA = (Fpto + sum(Fa) - design.weightA) ./ design.massA;
    
    dx(5,1) = vA;
    dx(6,1) = aA;
    
    % ********************************************************************
    % Buoy Acceleration, these calcs assume all of the mass is in the buoy,
    % this should be a reasonable assumption provided the mass of the
    % translator is not too significant in comparison

    Fexternal = FfeaVec + FaddT;
    
    % Determine the buoy/sea interaction forces and derivatives
    [dx, buoyancy_force, ...
        excitation_force_heave, ...
        excitation_force_surge, ...
        radiation_force_heave, ...
        radiation_force_surge, ...
        FBDh, FBDs, wave_height] = buoyodesim(t, x, dx, xBh, vBh, vBs, simoptions, Fexternal);

    % ************************************************************************
    
%     if any(isnan([x; dx; excitation_force_heave; excitation_force_surge; ...
%                   radiation_force_heave; radiation_force_surge; ...
%                   FBDh; FBDs; wave_height; Fexternal'; FfA; aA; ...
%                   dpsidxR; EMF; FA; FfeaVec'; xT; vT; ...
%                   FaddT'; FaddA'; Icoils'; xA; vA; xBh; xBs; vBh; vBs])) || ...
%         any(isinf([x; dx; excitation_force_heave; excitation_force_surge; ...
%                   radiation_force_heave; radiation_force_surge; ...
%                   FBDh; FBDs; wave_height; Fexternal'; FfA; aA; ...
%                   dpsidxR; EMF; FA; FfeaVec'; xT; vT; ...
%                   FaddT'; FaddA'; Icoils'; xA; vA; xBh; xBs; vBh; vBs]));
%               
%         keyboard;
%     
%     end
    
    % Now assign the outputs
    varargout{1} = dx;    

    % To record the forces
    if nargout > 1
        
        % per-coil induced voltage 
        varargout{2} = EMF;
        % total electromagnetic forces acting on armature
        varargout{3} = -FA;
        % vertical displacement of the translator  
        varargout{4} = xT;
        % vertical velocity of the translator
        varargout{5} = vT;
        % heave component of generator forces
        varargout{6} = Fexternal(1);
        % surge component of generator forces
        varargout{7} = Fexternal(2);
        % buoy drag forces in heave
        varargout{8} = FBDh;
        % buoy drag forces in surge
        varargout{9} = FBDs;
        % simple buoyancy force
        varargout{10} = buoyancy_force;
        % excitation force in heave
        varargout{11} = excitation_force_heave;
        % excitation force in surge
        varargout{12} = excitation_force_surge;
        % radiation force in heave
        varargout{13} = radiation_force_heave;
        % radiation force in surge
        varargout{14} = radiation_force_surge;
        % rate of change of flux linkage with displacement 
        varargout{15} = dpsidxR;
        % output armature friction
        varargout{16} = FfA;
        % ouput the armature forces
        varargout{17} = FaddA;
        % output the additional forces transmitted to the buoy from the
        % machine
        varargout{18} = FaddT;
        % output a breakdown of all the forces in the machine
        varargout{19} = FaddEBD;
        % wave height
        varargout{20} = wave_height;
        % phase resistance
        varargout{21} = diag(design.RPhase);
        
    end

end

