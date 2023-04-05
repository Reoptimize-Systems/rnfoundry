function [design, simoptions, T, Y, results] = evaluatedesign_ROTARY(design, simoptions)
% simulates and evaluates the design of a rotary pm machine
%
% Syntax
%
% [design, simoptions, T, Y, results] = evaluatedesign_ROTARY (design, simoptions)
%
% Description
%
% evaluatedesign_ROTARY simulates a rotary permanent magnet machine
% according to the specified simulation parameters, and evaluates the
% design according to supplied scoring criteria and data.
% evaluatedesign_ROTARY is generally intended to be called by a higher
% level function specific to a particular type of radial flux machine,
% e.g. evaluatedesign_RADIAL.
%
% Input
%
%  design - a structure containing the parameters of the machine or system
%   to be evaluated. The specific fields of this structure depend on the
%   system.
%
%  simoptions - structure used to control how the system is evaluated.
%   The simoptions structure must contain a substructure named ODESim. 
%
%   The ODESim structure can contain the following fields:
%
%   PreProcFcn : optional function handle or string containing the function
%     which will be run to generate data prior to running the simulation.
%     simfun will be passed the design and simoptions structure, and can
%     also be supplied with additional arguments by using the
%     'ExtraPreProcFcnArgs' parameter-value pair. The extra arguments must be
%     placed in a cell array. It must return two arguments which will
%     overwrite the design and simoptions variables.
%
%     The supplied function must have the calling syntax
%
%     [design, simoptions] = thefunction(design, simoptions, simfunarg1, simfunarg2, ...)
%
%     Where simfunarg1, simfunarg2, ... if supplied are the elements of the
%     a cell array, passed in using the Parameter-value pair
%     ExtraPreProcFcnArgs, e.g.
%
%     simulatemachine_AM(design, simoptions, 'ExtraPreProcFcnArgs', {1, 'another', [1,2;3,4]})
%
%   PostPreProcFcn : optional function handle or string containing a
%     function which will be run after simfun. finfun will also be passed
%     the design and simoptions structure, and can also be supplied with
%     additional arguments by using the 'ExtraPostPreProcFcnArgs'
%     parameter-value pair. The extra arguments must be placed in a cell
%     array. It must return two arguments which will overwrite the design
%     and simoptions variables.
%
%     The supplied function must have the calling syntax
%
%     [design, simoptions] = thefunction(design, simoptions, finfunarg1, finfunarg2, ...)
%   
%     Where finfunarg1, finfunarg2, ... if supplied are the elements of the
%     a cell array, passed in using the Parameter-value pair
%     ExtraPostPreProcFcnArgs, e.g.
%
%     simulatemachine_AM(design, simoptions, 'ExtraPostPreProcFcnArgs', {1, 'another', [1,2;3,4]})
%
%   EvalFcn : function handle or string containing the function which will
%     be evaluated by the ode solver routines to solve the system of
%     equations. see the ode solvers (e.g. ode45, ode15s) for further
%     information on how to create a suitible function.
%
%   SolutionComponents : Structure containing information on the various
%     components of the system to be solved. This structure will contain
%     fields with names corresponding to the solution components. Each of
%     these fields is then also a structure with information about each of
%     these solution components. The information is contained in the
%     following fields:
%
%     SolutionIndices : these are the indices of the ode system
%       coresponding to this component.
%
%     InitialConditions : this is the initial conditions for variables
%       associated with this component.
%
%     AbsTol : the absolute tolerances of the variables associated with
%       this  component. If AbsTol is not supplied for every solution
%       component it will be ignored for all components.
%
%     The initial conditions for the ODE solution and possibly the absolute
%     tolerances on all variables are assembled from the information in
%     this structure.
%
%     In addition, the field OutputFcn may be present. This contains a
%     function handle or string representing a function to be run on the
%     variables associated with the solution component after each
%     successful time step.
%
%     SolutionComponents can also contain a field named 'NestedSim' for the
%     purposes of setting up a multi-rate solution. If present, it contains
%     solution component information intended for a nested, higher time
%     step rate simulation (implemented using the ode.odesolver class),
%     running between the steps of the 'outer' or top level simulation. The
%     NestedSim sturcture is the same format as the top level
%     SolutionComponents structure and contains the solution component
%     information for the variables in the lower level simulation. The
%     NestedSim structure can also contain another NestedSim structure and
%     any number of levels are supported.
%
%   PostAssemblyFcn : function handle or string containing a function which
%     will be run after the solution components have been read and fully
%     assembled.
%
%   PostSimFcn : function handle or string containing a function which will
%     be run after the simulation has been completed by the ode solver.
%     resfun must take the T and Y matrices, as generated by the ode
%     solver, and the design and simoptions arguments in that order.
%     PostSimFcn can also be supplied with additional arguments by using
%     the 'ExtraPostSimFcnArgs' parameter-value pair. The extra arguments
%     must be placed in a cell array. It must return two arguments, one of
%     which is a results variable containing results of interest to the
%     user, the other of which overwrites the design variable.
%
%     The supplied function must have the calling syntax
%
%     [results, design] = htefunction(T, Y, design, simoptions, odearg1, odearg2, ..., resfunarg1, resfunarg1, ...);
%
%     Where odearg1, odearg2, ..., resfunarg1, resfunarg1, ... if supplied are the elements of the
%     two cell arrays, passed in using the Parameter-value pairs
%     ExtraEvalFcnArgs and ExtraPostSimFcnArgs respectively, e.g.
%
%     simulatemachine_AM ( design, simoptions, ...
%                          'ExtraEvalFcnArgs', {1, true}, ...
%                          'ExtraPostSimFcnArgs', {2, false} )
%
%   Solver : function handle or string specifying the ode solver to use,
%     if not supplied, for Matlab the default is 'ode15s', and for Octave
%     'ode2r'.
%
%   Split : (scalar integer) if this field is present in the structure
%     it indicates that the evaluation of the system of differential
%     equations should be split into manageable chunks, useful for
%     long-running simulations which use substantial memory. The value of
%     ODESim.Split is the desired initial number of chunks into which
%     evaluation will be split. If the system runs out of memory during
%     simulation, the number of blocks will be increased and simulation
%     reattempted this will be attempted at most 4 times. If ODESim.Split
%     is present, the field SplitPointFcn must also be supplied, documented
%     below.
%
%   SplitPointFcn : (string|function handle) if ODESim.Split is provided
%     this field must also be present which should contain a string or
%     function handle. This function will be called at each break in the
%     integration and must have the following syntax:
%
%     results = spfcn (flag, results, sol, design, simoptions)
%   
%     For further information on creating a split point function, see the
%     help for 'odesplit.m'.
%
%   OutputFcn : Top level OutputFcn to be run after each successful time
%     step. If not supplied, the function odesimoutputfcns_AM is used. This
%     function looks for a SolutionComponents structure in
%     simoptions.ODESim and runs the output function (if present) for each
%     solution component with the appropriate inputs.
%
%   Events : Event function to run after every time step to determine if a
%     termination event has occured. See the help for the ode solvers (e.g.
%     ode45) to learn more about the event function.
%
%   Vectorized : true/false flag indicating if the ODE solution function is
%     vectorized
%
%   MaxStep : The maximum allowed time step size
%
%   InitialStep : A suggested initial time step size for the ode solution
%
% Output
%
%  design - the input design structure returned with any modifications made
%   by the performed simulations. The design structure will generally have
%   a number of fields added, containing the main results from the
%   evaluation.
%
%  simoptions - the input design structure returned with any modifications made
%   by the performed simulations. The modification can include the addition
%   of fields containing default options used which were not provided by
%   the user. The exact fields added are simulaiton dependant.
%
%  T - output time vector as produced by ode solver functions, e.g. ode45
%
%  Y - output solution vector as produced by ode solver functions, e.g.
%   ode45
% 
%  results - the results as produced by the supplied function in
%   ODESim.PostSimFcn. Typically a structure containig fields which are the
%   outputs of various quantities at each time step in the simulation.
%
% See also: simulatemachine_AM
%

% Copyright Richard Crozier 2012

    
    % simulate the machine
    [T, Y, results, design, simoptions] = simulatemachine_AM (design, simoptions);
    
%     % evaluate the structure, unless we are told to skip it
%     if ~simoptions.Evaluation.SkipStructural
%         [maxzdef, maxstress, design] = evaluatestructure_RADIAL(design, simoptions);
%     else
%         maxzdef = 0; 
%         maxstress = 0;
%     end
% 
%     % copy some results into the design structure
%     design.MaxDeflection = maxzdef;
%     
%     design.MaxStress = maxstress;
    
end