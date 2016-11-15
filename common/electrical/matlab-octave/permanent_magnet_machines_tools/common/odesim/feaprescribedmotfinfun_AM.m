function [design, simoptions] = feaprescribedmotfinfun_AM (design, simoptions)
% common function for setting up an ode simulation of an electrical machine
% based on full finite element analysis at every time step
%
% Syntax
%
% [design, simoptions] = feaprescribedmotfinfun_AM (design, simoptions)
%
%

    % set up the numerical derivative calculation of the flux linkage using
	% an odederiv object
    flodederiv = odederiv (0, simoptions.ODESim.SolutionComponents.PhaseFluxLinkages.InitialConditions', ...
        'DoPlot', true);
    
	simoptions.ODESim.SolutionComponents.PhaseFluxLinkages.NumericalDerivatives = ...
        flodederiv;
    
    % we also create an output function for these components to update the
    % derivatives after each time step
    simoptions.ODESim.SolutionComponents.PhaseFluxLinkages.OutputFcn = ...
        @flodederiv.outputfcn;
    
     simoptions.ODESim.SolutionComponents.PhaseFluxLinkages.ResetFcn = ...
        @flodederiv.reset;
    
    simoptions.ODESim = setfieldifabsent (simoptions.ODESim, 'OutputFcn', 'odesimoutputfcns_AM');
    
%     simoptions.ODESim.Solver = @odef1maxstep;
%     simoptions.ODESim.Solver = @rkfixed;
%     simoptions.ODESim.Solver = @ode15s;

end