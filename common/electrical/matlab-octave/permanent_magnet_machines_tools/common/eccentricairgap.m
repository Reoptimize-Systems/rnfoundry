function g = eccentricairgap (theta, Ri, Ro, dg, thetadisp)
% calculates the air gap length around an eccentric rotor-stator
%
% Syntax
%
% g = eccentricairgap (theta, Ri, Ro, dg)
% g = eccentricairgap (..., thetadisp)
%
% Description
%
% eccentricairgap calculates the new value of the air gap length around the
% circumference of a machine with an eccentrically mounted rotor/stator.
% All inputs can be scalar values or vectors/matrices the same size as each
% other.
%
% Input
%
%  theta - angular positions(s) around the stator in radians. 
%
%  Ri - inner component radius
%
%  Ro - outer component radius
%
%  dg - the magnitude of the relative displacement of the two components
%
%  thetadisp - (optional) angle of displacement in radians. This is the
%    direction of the displacement causing the eccentricity. Default is
%    zero if not supplied, meaning the displacement is in the positive 'x'
%    direction.
%
% Output
%
%  g -the new air gap size at the corresponding value(s) of theta
%
% See Also: analyticaleccentricump.m
%

    if nargin < 5
        thetadisp = 0;
    end
    
    check.isNumericScalar (thetadisp, true, 'thetadisp');
    
    check.multicheck (@isnumeric, 'all inputs to eccentricairgap must be numeric matrices', '', ...
        theta, Ri, Ro, dg )
    
    assert ( samesizeorscalars (theta, Ri, Ro, dg), ...
        'theta, Ri, Ro and dg must all be the same size or scalar values' );
    
    
    g = Ro - sqrt ( realpow (dg + Ri.*cos (theta - thetadisp), 2) ...
                     + realpow (Ri.*sin (theta - thetadisp), 2) ) ;

end