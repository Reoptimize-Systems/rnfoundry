function [design, simoptions] = circuitprops_linear(design, simoptions)
% circuitprops_linear: calculates and completes various circuit parameters
% for a linear machine for evaluation in a ode simulation
%
% Syntax
% 
% [design, simoptions] = circuitprops_linear(design, simoptions)
%
% Input
%
% design is a structure containing various design parameters of the linear
% machine necessary for calculating some circuit properties. The design
% structure must contain the following fields:
% 
%   CoilResistance - resistance of a single coil winding
%   RgVRc - ratio of coil resistance to grid/load resistance
%
% design can also optionally contain the following fields:
%
%   CoilInductance - Inductance of a single coil, if not present the
%     phase inductance will be set to 1e-4
%   LgVLc - ratio of grid inductance to col inductance, ignored if
%     CoilInductance not present
%   Branches - Number of parallel branches of coils per phase
%   CoilsPerBranch - Number of series coils in each parallel branch
%     NB: (Both CoilsPerBranch and Branches must be present for either to
%     be used) If neither are specified, all coils are assumed to be in
%     parallel. In this case the field 'PowerPoles' must be present which
%     will be the number of parallel coils per phase)
%   
%
% Output
%
% The design matrix will be populated with new fields depending on it's
% inital contents.
% 
% In all cases the field 'R' will be added which contains the phases
% resistance.
%
% If not present previously the fields 'Branches' and 'CoilsPerBranch' will
% be added.

    % we now call the common function
    [design, simoptions] = circuitprops_AM(design, simoptions);

end