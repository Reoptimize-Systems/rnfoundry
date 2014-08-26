function varargout = prescribedmotodetorquefcn_ROTARY(t, x, design, simoptions)
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
        varargout{1} = {'ydot', 'dpsidthetaR', 'EMF', 'thetaT', 'omegaT', 'Tqpto', 'TqaddE', 'TqaddEBD', 'RPhase'};
        return;
    end

    % Change the x members into more useful variables names, MATLAB will
    % optimise away any memory penalty associated with this I think    
    Iphases = x(1:design.Phases);
    
    Icoils = Iphases ./ design.Branches;

    % Initialize dx with zeros
    dx = zeros(size(x));
    
    % get the velocity and position at the current time
    [thetaE, omegaE] = prescribedmotomegatheta(t, simoptions);

    % determine the machine outputs
    % Get the EMF and forces using the core machine simulation function
    [Tqeff, Tqreac, EMF, dpsidthetaR, design] = ...
        machineodesim_AM(design, simoptions, thetaE, 0, omegaE, 0, Icoils);
    
    % find the derivative of the coil current (solving the differential
    % equation describing the simple output circuit)
    dx(1:design.Phases,1) = circuitode_linear(Iphases, EMF, design);
     
    % call the supplied additional force function
    [TqaddE, TorqueBD] = ...
        feval(simoptions.torquefcn, design, simoptions, thetaE, omegaE, EMF, Iphases, simoptions.forcefcnargs{:});
    
    % ************************************************************************

    % Now assign the outputs
    varargout{1} = dx;

    % To record the forces
    if nargout > 1
        
        % rate of change of flux linkage with displacement 
        varargout{2} = dpsidthetaR;
        % per-coil induced voltage 
        varargout{3} = EMF;
        % vertical displacement of the translator  
        varargout{4} = thetaE;
        % vertical velocity of the translator
        varargout{5} = omegaE;
        % total electromagnetic forces acting on effector
        varargout{6} = Tqeff;
        % output the additional forces transmitted to the primemover
        varargout{7} = TqaddE;
        % output a breakdown of the additional forces transmitted to the
        % prime mover
        varargout{8} = TorqueBD;
        % output the phase resistance at each time step
        varargout{9} = diag(design.RPhase)';
        
        
    end

end

