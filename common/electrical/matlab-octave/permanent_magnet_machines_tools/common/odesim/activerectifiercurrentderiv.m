function idot = activerectifiercurrentderiv(I, E, R, M, D, Vo, VnN)
% calculates the current derivatives in a 3 phase rectifier and motor
%
% Syntax
% 
% idot = activerectifiercurrentderiv(I, E, R, M, D, Vo, VnN)
%
% Description
%
% Solves the right hand side of the differential equations describing the
% average behaviour of a three-phase active rectifier.
%                              
% For an explanation, see Suntio, T., Messo, T and Puukko, J., 2018.
% 'Power Electronic Converters, Dynamics and Control in Conventional and
% Renewable Energy Applications', First Edition, Chapter 17: Dynamic
% Modeling of Thre-Phase Active Rectifiers, Wiley-VHC Verlag GmbH & Co.
% KGaA.                    
%
% Input
%
%  I - is a vector of values of the currentin the the machine Phases at the
%   current time.
%
%  E - is a vector of values of the emfs in each phase at the current time
%
%  R - is an (n x n) matrix with diagonals the phase resistances of the
%   machine, i.e.
%
%   R = [ Ra  0   0;
%         0   Rb  0;
%         0   0   Rc ];
%
%  M - is an inductance matrix with diagonal terms the main inductances of the
%   Phases, and off-diagonal terms the mutual inductances between Phases
%
%   M = [ La   Lba  Lca;
%         Lba  Lb   Lcb;
%         Lba  Lcb  Lc  ];
% 
%   Often it is assumed that: La = Lb = Lc = Lp, and Lba = Lca = Lcb = Mpp
%   such that
%
%   M = [ Lp   Mpp  Mpp;
%         Mpp  Lp   Mpp;
%         Mpp  Mpp  Lp ];
%
%
%  D - (3 x 1) vector od duty ratios for the switches for each phase, i.e.
%    [da; db; dc]
%
%  Vo - Desired DC output voltage
%
%  VnN - Voltage between ground and base rail of rectifier
%
%

    % Calculate di/dt
    idot = ( (E(:) - R * I(:) - D*Vo + VnN)' / M )';

end