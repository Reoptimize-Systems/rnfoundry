function modifiedIvars = combinedsheatandcoilI_TM(Icoil, IVars)
% Determines the appropriate inner coil radius to give a cross-sectional
% area in a tubular machine which has the equivalent second moment of area
% of a beam made of the same material as the outer sheath surrounding the
% coils
%
% Input:
%  
%   Icoil - scalar value of the desired second moment of area of the coil
%           section, used to determine the correct inner radius
%
%   IVars - (2 x 4) matrix of values for calculating the second moment of
%           area of the stator and translator sections, the second row of
%           values is used for both stator sections if a double sided
%           machine is being investigated.
%
%           IVars(1,1): Ro, outer coil radius
%           IVars(1,2): Ri, inner coil radius
%
%           IVars(1,1): Ra, outer sheath radius
%           IVars(1,2): Ro, inner sheath radius
%
% Output:
%
%   modifiedIvars - (2 x 2) matrix of values descirbiing the new combined
%                   section
%
%                   IVars(1,1): Ra
%                   IVars(1,2): Ri
%

    % First we determine what the inner radius of the coil section would be
    % if it were made of the same material as the sheath
    Ri = ((IVars(1,1).^4) - (4*Icoil / pi)) .^ (0.25);
    
    % The IVars for the new combined beam are the outer radius of the
    % sheath and the inner radius of the modified coil section
    modifiedIvars = [IVars(2,1) Ri];

end