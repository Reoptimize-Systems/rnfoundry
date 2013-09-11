function emf = peakemfest_ROTARY(peakfl, omega, pp)
% peakemfest_ROTARY: estimates the peak voltage of a rotary machine at a
% specified angular velocity given the peak flux linkage and the number of
% pole pairs
%
% Syntax
%
% emf = peakemfest_linear(peakfl, omega, pp)
%
%
% Input
%
%  peakfl - peak flux linkage in a coil or phase
%
%  omega - machine angular velocity
%
%  pp - the number of pole pairs in the machine
%
% Output
%
%  emf - the estimated peak emf at the specified angular velocity in a coil
%   or phase
%

    % v = d/t, v = (2Wp)/t, v/(2*Wp) = 1/t = f
    % elec omega = omega * pole pairs
    omegaE = omega * pp;
    
    % emf = d lambda / d t = (d/dt)(lambda sin(omega t))
    %     = lambda omega cos(omega t)
    %
    % peak emf = lambda omega
    emf = peakfl * omegaE;

end