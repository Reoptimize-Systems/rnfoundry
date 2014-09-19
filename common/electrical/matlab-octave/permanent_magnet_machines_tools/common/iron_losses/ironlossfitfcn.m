function res = ironlossfitfcn(fitvars,f,B,coreloss)
% computes the residuals of a fit to experimental epstein sample test iron
% loss data to the standard iron loss equation
%
% Syntax
%
% res = ironlossfitfcn(fitvars,f,B,coreloss)
%
% Description
%
% The main purpose of this function is as an objective function for the
% function LMFnlsq.m to produce a fit to loss data for a material generated
% using epstein sample tests with sinusoidal applied magnetic fields of
% various frequencies and magnitudes.
%
% ironlossfitfcn fits the coefficients kh, kc, ke and beta in a function of
% the form:
%
%             kh*f*B^beta + kc*(f*B)^2 + ke*(f*B)^1.5
%
% Where kh is the hysteresis loss coefficient, kc is the eddy current loss
% coefficient and ke is the excess, or anomalous loss coefficient. See [1]
% for more information on this function and its use.
%
% [1] G. Bertotti, General properties of power losses in soft
% ferromagnetic materials, IEEE Transactions on Magnetics, vol. 24, no. 1,
% pp. 621-630, 1988.
%
% Input
%
%  firvars - a four element vector containing values of kh, kc, ke and beta
%    to be tested, i.e. the vector [kh, kc, ke, beta].
%
%  f - a matrix of frequency values at which the tests were performed.
%
%  B - a matrix of peak flux density values at which the tests were
%    performed.
%
%  coreloss - a vector of values of the losses at the corresponding
%    positions in f and B
%
% Output
%
%  res - the residuals between the function and the actual core loss data
%
% Example
%
% % choose a starting point for the fit variables
% Xo = [0.02, 0.0001, 0.001, 1.5];
% 
% options = LMFnlsq('default');
% options = LMFnlsq(options, 'Display', 1, 'XTol', 1e-9);
% 
% % fit values of kh, kc, ke and beta
% [xf, SS, cnt, res, XY] = LMFnlsq(@(fitvars) ironlossfitfcn(fitvars,f,B,coreloss), Xo, options);
%

% Copyright Richard Crozier 2012

    fitvars = abs(fitvars);

    kh = fitvars(1);
    kc = fitvars(2);
    ke = fitvars(3);
    beta = fitvars(4);
    
    % res = FUN(x) - y
    fitcoreloss = kh .* f(:) .* B(:).^beta + kc .* (f(:).*B(:)).^2 + ke .* (f(:).*B(:)).^1.5;
    
    res = real(fitcoreloss - coreloss(:));
    
end
