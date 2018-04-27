function varargout = test_activerectifiercurrentderiv_evalfcn (t, x, design, simoptions)
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
        varargout{1} = { 'ydot', 'thetaT', 'omegaT', ...
                         'Vabc', 'D', 'Idq', 'Tqpto', 'EMF' };
        return;
    end

    % Change the x members into more useful variables names, MATLAB will
    % optimise away any memory penalty associated with this I think    
    Iphases = x(1:design.Phases);
    
    Icoils = Iphases ./ design.Branches;

    % Initialize dx with zeros
    dx = zeros (size (x));
    
    % get the velocity and position at the current time
    [thetaE, omegaE] = prescribedmotomegatheta (t, simoptions);

    % determine the machine outputs
    % Get the EMF and forces using the core machine simulation function
    [Tqeff, ~, EMF,  ~, design] = ...
        machineodesim_AM (design, simoptions, thetaE, 0, omegaE, 0, Icoils);

    theta_flux = thetaE .* design.Poles / 2 + pi();

    Vo = design.MachineSidePowerConverter.Vdc;
    
    if t <= simoptions.FOCTApp
        Vabc = zeros (size (Iphases));
        D = zeros (numel (Iphases), 1);
    else

        Vabc = dq02abc ( [simoptions.Vsdref; simoptions.Vsqref], theta_flux, true );
        
        D =  Vabc ./ Vo;
        
    end
    
    VnN = 0;
    
    dx(1:design.Phases,1) = activerectifiercurrentderiv ( Iphases(:), ...
                                                          EMF(:), ...
                                                          design.MachineSidePowerConverter.REquivalent, ...
                                                          design.L, ...
                                                          D, ...
                                                          Vo, ...
                                                          VnN );
     
    
    % ************************************************************************

    % Now assign the outputs
    varargout{1} = dx;

    % To record the forces
    if nargout > 1
        
        % vertical displacement of the translator  
        varargout{2} = thetaE;
        % vertical velocity of the translator
        varargout{3} = omegaE;
        
        varargout{4} = Vabc;
        
        varargout{5} = D;
        
        varargout{6} = abc2dq0 (Iphases, theta_flux, true);
        
        varargout{7} = Tqeff;
        
        varargout{8} = EMF;
        
    end

end

