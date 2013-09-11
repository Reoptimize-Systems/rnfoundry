function [histloss, eddyloss, excessloss] = ...
    softferrolossrectregionvarxpartcalc(Bx, By, Bz, Hx, Hy, Hz, kc, ke, beta, xstep, dx, dy, dz)
% part calculates the values of the hysteresis, eddy current and excess
% loss in a cuboidal region of soft ferromagnetic material
%
% Syntax
%
% [histloss, eddyloss, excessloss] = softferrolossrectregpartcalc(Bx, By, Bz, Hx, Hy, Hz, dx, dy, dz, kc, ke, beta)
%
% Input
%
%  Bx - 
% 
%  By - 
% 
%  Bz -  
% 
%  Hx -  
% 
%  Hy -  
% 
%  Hz -  
% 
%  kc -  
% 
%  ke -  
% 
%  beta - 
%
%  dx -  
% 
%  dy -  
% 
%  dz
%
% Output
%
%  histloss -  
% 
%  eddyloss -  
% 
%  excessloss - 
%
% Description
%
% Calculates part of the formula required to perform loss calculations
% using the time-domain method laid out in [1]. The resulting formula
% assumes the direction of motion is constrained so that dB/dt can be
% reasonably represented as:
%
% dB/dx * dx/dt = dB/dx * velocity
%
% [1] D. Lin, P. Zhou, W. N. Fu, Z. Badics, and Z. J. Cendes, “A Dynamic
% Core Loss Model for Soft Ferromagnetic and Power Ferrite Materials in
% Transient Finite Element Analysis,�? IEEE Transactions on Magnetics, vol.
% 40, no. 2, pp. 1318--1321, Mar. 2004.

    if nargin == 11
        dV = dx;
    elseif nargin == 13
        dV = dx .* dy .* dz;
    else
        error('Incorrect number of arguments.')
    end

    % the direction of motion is assumed to be in the x direction. Therfore
    % we first find the derivatives of each B component with respect to
    % this direction. dx is assumed to be a scalar, therefore Bx, By and Bz
    % must all be sampled at
%     dBxVdx = bgrad(Bx, xstep);
%     dByVdx = bgrad(By, xstep);
%     dBzVdx = bgrad(Bz, xstep);
    dBxVdx = gradient(Bx, xstep);
    dByVdx = gradient(By, xstep);
    dBzVdx = gradient(Bz, xstep);
    
    % generate the part calculation of the hysteresis losses for the
    % region, these must be multiplied by the velocity to get the actual
    % losses
    histloss = abs(Hx .* dBxVdx).^(2/beta) ...
               + abs(Hy .* dByVdx).^(2/beta) ...
               + abs(Hz .* dBzVdx).^(2/beta);
           
	% replace infinite values with realmax
	histloss(isinf(histloss)) = realmax;
        
    histloss = dV .* (histloss .^ (beta / 2));
    
    histloss = sum(histloss, 3);
    
    histloss = sum(histloss, 1);
    
    % generate the part calculation of the eddy current losses for the
    % region, these must be multiplied by the square of the velocity (v^2)
    % to get the actual losses
    Bderivsquares = realpow(dBxVdx, 2) + realpow(dByVdx, 2) + realpow(dBzVdx, 2);
    
    eddyloss =  dV .* (kc / (2 * pi^2)) .* Bderivsquares;
    
    eddyloss = sum(eddyloss, 3);
    
    eddyloss = sum(eddyloss, 1);
    
    % generate the part calculation of the excess (sometimes innacurately
    % known as anomalous) losses for the region, these must be multiplied
    % by the cube of the square root of the velocity (sqrt(v)^3) to get the
    % actual losses
    
    % Ce is a constant found by the numerical solution of 
    % <latex>
    % \[ Ce = (2*\pi)^{1.5} \frac{2}{\pi} \int_{0}^{\frac{\pi}{2}} \cos^{1.5}(\theta) d \theta \]
    % <latex>
    % This is described in [1]
    Ce = 8.763363;
    
    excessloss = dV .* (ke / Ce) .*  realpow(Bderivsquares, 0.75);
    
    excessloss = sum(excessloss, 3);
    
    excessloss = sum(excessloss, 1);

end



function Bgrad = bgrad(Bmat, xstep)

    if isscalar(xstep)
        xstep = 0:xstep:(xstep*(size(Bmat, 2)-1));
    end
    
    % switch the rows with the columns as interp1, which we will use to
    % generate a piecewise polynomial fitted to the data, interpolates down
    % the first dimension of the matrix (the rows), but the data is
    % provided with x along the 2nd dimenaion, (the columns). This is
    % identical to doing the transpose of each slice of the 3D matrix in
    % the 3rd dimension
    Bmat = permute(Bmat,[2,1,3]);
    
    % generate the piecewise polynomial fits to the data using interp1
    pp = interp1(xstep.', Bmat, 'pchip', 'pp');
    
    % now construct a new set of polynomials which are the derivatives of
    % the originals
    [breaks,coefs,npieces,order,targetdim] = unmkpp(pp);
    
    % make the polynomials that describe the derivatives
    dpp = mkpp(breaks, bsxfun(@times, coefs(:,1:order-1), order-1:-1:1), targetdim);
    
    % record that the input data was oriented according to INTERP1's rules.
    % Thus PPVAL will return its values oriented according to INTERP1's
    % convention
    dpp.orient = 'first';
    
    % calculate the gradients using these new polynomials
    Bgrad = ppval(dpp, xstep);
    
    % inverse permute to get the original direction
    Bgrad = ipermute(Bgrad, [2,1,3]);

end


