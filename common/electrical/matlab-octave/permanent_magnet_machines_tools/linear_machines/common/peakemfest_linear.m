function emf = peakemfest_linear(peakfl, velocity, polewidth)
% peakemfest_linear: estimates the peak voltage of a linear machine at a
% specified velocity given the peak flux linkage and pole width
%
% Syntax
%
% emf = peakemfest_linear(peakfl, velocity, polewidth)
%


    % v = d/t, v = (2Wp)/t, v/(2*Wp) = 1/t = f
    % omega = 2 pi f = 2 pi v / (2*Wp)
    omega = 2 * pi * velocity / (2 * polewidth);
    
    % emf = d Psi / d t = (d/dt)(Psi sin(omega t))
    %     = Psi omega cos(omega t)
    %
    % peak emf = Psi omega
    emf = peakfl * omega;


end