function varargout = evaluatefestructure_ACPMSM(IVars, design, options)
% evaluatefestructure_ACPMSM: evaluates whether a given beam can withstand
% the forces in an air-cored linear permanent magnet synchronous machine
% design using the FAESOR finite element toolbox
%
%                                                  
    % We will calculate the deflection of the I-Beams, stopping when we reach
    % an I-Beam that meets the specifications, as this will be the smallest
    % volume with which this can be achieved

%   IVarsCell - (1 x 2) cell array of values, each of which contains a matrix
%           of values for calculating the second moment of area of the
%           stator and translator sections repectively, the the values in
%           the second cell are used for both stator sections if a double
%           sided machine is being investigated.
%
%           First cell of IVars must be the beams supporting the inner and
%           outer 
%

    % complete the IVars cell array by adding the appropriate information
    % from the design structure
    IVarsCell = {IVars, design.InnerStructureBeamVars};
          
    % Determine what the air-gap is with the structure loaded by the maxwell
    % stresses
    new_g = feairgapclosure_ACPMSM(design, options, IVarsCell);

    % Check if beam is successful or a failure
    if min(new_g(:)) >= options.gfactor * design.g
        varargout{1} = 1;
        varargout{2} = IVars;
    else
        varargout{1} = 0;
        varargout{2} = IVars;
    end

end