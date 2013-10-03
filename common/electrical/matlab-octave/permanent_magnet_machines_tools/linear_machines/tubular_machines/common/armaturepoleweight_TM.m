function varargout = armaturepoleweight_TM(WpVRm, RoVRm, Rm, g, WcVWp, cufill, copperDensity, RaVRo, steelDensity)
% Calculates the weight of one pole of the armature of a tubular machine.
% If information about the steel is omitted, it is assumed the machine is
% air-cored.
%
% Arguments: (input)
%
%   WmVWp - scalar value of Wm/Wp Ratio for machine to be evaluated
%
%   WpVRm - scalar value of Wp/Rm Ratio for machine to be evaluated
%
%   RsiVRso - scalar value of Rsi / Rso, the ration of the shaft inner
%             diameter to the shaft outer diameter
%
%   RsoVRm - sclar value of Rso/Rm, the ratio of the shaft outer diameter
%            to the translator radius
%
%   Rm - Radius of translator
%
%   g - the air gap of the machine
%
%   WcVWp - The coil width to pole width ratio, it is assumed this is <=
%           1/3
%
%   cufill - the copper fill factor in the windings
%
%   copperDensity - the density of the copper windings
%
%   RaVRo - The outer atrmature sheath radius to inner sheath radius ratio,
%           if not supplied, it is assumed it is an air-cored machine
%
%   steelDensity - density of steel in the machine, if not supplied, it is
%                  assumed it is an air-cored machine
%
% Output:
%
%   totalWeight - the total weight of copper and steel (if present) in a
%                 pole. If an air-cored machine this will simply be the
%                 weight of copper and no other return values will be
%                 supplied. If it is a slotless machine with a
%                 ferromagnetic sheath the following two return values will
%                 also be supplied
%
%   copperWeight - the total wieght of copper making up one pole
%
%   steelWeight - the total weight of steel in the sheath of one pole

    Wp = WpVRm * Rm;
    Wc = WcVWp * Wp;
    Ro = RoVRm * Rm;
    Ri = Rm + g;
    
    G = 9.81;
    
    Vcopper = 3 * Wc * pi * (Ro^2 - Ri^2) * cufill;
    
    copperWeight = Vcopper * copperDensity * G;
    
    varargout{1} = copperWeight;

    if nargin >= 7
        
        Ra = RaVRo * Ro;
        Vsteel = Wp * pi * (Ra^2 - Ro^2);
        steelWeight = Vsteel * steelDensity * G;
        
        varargout{1} = copperWeight + steelWeight;
        varargout{2} = copperWeight;
        varargout{3} = steelWeight;

    end

    
end