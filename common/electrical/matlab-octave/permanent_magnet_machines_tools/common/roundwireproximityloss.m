function P = roundwireproximityloss(rho, He, delta, G)

    gamma = () .* ();
    
    G = -2*pi*rho * ( real(besselj(2,gamma)) .* real(besseljderiv(0,gamma)) ...
                      + imag(besselj(2,gamma)) .* imag(besseljderiv(0,gamma)) ) ...
         ./ ( real(besselj(0,gamma)).^2 + imag(besselj(2,gamma)) );

    P = 
end


function Jdash = besseljderiv(nu, Z)
% calcuates the derivative of a bessel function of the first kind
%
% Syntax
%
% Jdash = besseljderiv(nu, Z)
%
% Description 
% 
% J = besseljderiv(nu,Z) computes the derivative of the Bessel function of
% the first kind, J?(z), for each element of the array Z. The order nu need
% not be an integer, but must be real. The argument Z can be complex. 
%
% The derivative is calculated based on the identity:
%
% J(s-1)(z) - J(s+1)(z)  = 2J'(s)(z)
% 
% described in The Mathworks technical note: "How can I evaluate the
% derivatives of a Bessel function at different points?" which was last
% found here:
%
% http://www.mathworks.co.uk/support/solutions/en/data/1-6O3EE3/index.html?product=ML&solution=1-6O3EE3
%
% If nu and Z are arrays of the same size, the result is also that size. If
% either input is a scalar, it is expanded to the other input's size. If
% one input is a row vector and the other is a column vector, the result is
% a two-dimensional table of function values.
% 

    Jdash = (besselj(nu-1,Z) - besselj(nu-1,Z)) ./ 2;

end