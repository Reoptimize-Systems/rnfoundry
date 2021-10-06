function varargout = prescribedmotodetorquefcn_dqactiverect_ROTARY (t, x, design, simoptions)
% solves the right had side of the differential
% equations for a machine moving with a prescribed motion
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

    if nargin == 0
        % get the output names
        varargout{1} = { 'ydot', 'EMF', 'thetaT', 'omegaT', ...
                         'Tqpto', 'D', 'Vsref', 'PI_vals' };
        return;
    end

    % Change the x members into more useful variables names, MATLAB will
    % optimise away any memory penalty associated with this I think    
    Idq = x;
    
    % get the velocity and position at the current time
    [thetaE, omegaE] = prescribedmotomegatheta (t, simoptions);

    % determine the machine outputs
    % Get the EMF and forces using the core machine simulation function
    Tqpto = 1.5 * design.Poles/2 * design.FluxLinkagePhasePeak * Idq(2);

    omega_flux = omegaE .* design.Poles / 2;

    Vo = design.MachineSidePowerConverter.Vdc;

    if t <= simoptions.FOCTApp
        Vsdqref = [ 0; 0 ];
        D = zeros (2, 1);
    else

        dt = calcDt (design.FOControl.PI_d, t);
        Vsdqref(1,1) = calculate ( design.FOControl.PI_d, Idq(1), simoptions.Isdref, dt );
        Vsdqref(2,1) = calculate ( design.FOControl.PI_q, Idq(2), simoptions.Isqref, dt );
          
        if any (abs(Vsdqref) > Vo)
            
            scalefac = Vo / max(abs(Vsdqref));
        
            Vsdqref = Vsdqref * scalefac;
            
        end
        
        D = Vsdqref ./ Vo;
        
    end
    
    EMF = [ 0; design.FluxLinkagePhasePeak * omega_flux ];
    
    dx = dqactiverectifiercurrentderiv ( Idq, ...
                                         EMF, ...
                                         design.MachineSidePowerConverter.REquivalent(1), ...
                                         design.L(1), ...
                                         D, ...
                                         Vo, ...
                                         omega_flux );
     
    % ************************************************************************

    % Now assign the outputs
    varargout{1} = dx;

    % To record the forces
    if nargout > 1
        
        % per-coil induced voltage 
        varargout{2} = EMF;
        % vertical displacement of the translator  
        varargout{3} = thetaE;
        % vertical velocity of the translator
        varargout{4} = omegaE;
        % output the additional forces transmitted to the primemover
        varargout{5} = Tqpto;
        varargout{6} = D;
        varargout{7} = Vsdqref;
        varargout{8} = [design.FOControl.PI_d.error, design.FOControl.PI_q.error, design.FOControl.PI_d.integral, design.FOControl.PI_q.integral ];
        
    end

end

