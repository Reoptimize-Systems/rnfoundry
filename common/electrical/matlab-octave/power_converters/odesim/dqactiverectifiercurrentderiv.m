function dy = dqactiverectifiercurrentderiv(Idq, Vdq, R, L, D, VC, omegas)
% calculates the dq state derivatives in a 3 phase rectifier and motor
%
% Syntax
% 
% dy = dqactiverectifiercurrentderiv(Idq, Vdq, R, L, D, VC)
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
%  Idq - is a vector of values of the dq current in the the machine Phases
%   at the current time.
%
%  Vdq - is a vector of values of the voltages applied to each phase at the 
%   current time
%
%  R - is the scalar value of the phase resistances
%
%  L - 
%
%  D - (3 x 1) vector od duty ratios for the switches for each phase, i.e.
%    [da; db; dc]
%
%  VC - Desired DC output voltage
%
%

    dy = [ ( Vdq(1) / L ) - (R / L)*Idq(1) + omegas*Idq(2) - (D(1)/L) * VC;
           ...
           ( Vdq(2) / L ) - (R / L)*Idq(2) - omegas*Idq(1) - (D(2)/L) * VC;
         ];
    
end