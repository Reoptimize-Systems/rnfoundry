function Z = phaseimpedance_AM (design, omega)
% calculates the machine phase impedance seen by an external 3 phase load
%
% Syntax
%
% Z = phaseimpedance_AM (design, omega)
%
% Description
%
% calculates the machine phase impedance seen by an external three phase
% load connected to the machine terminals. Uses the cyclic inductance which
% incorporates the mutual inductance between phases. The cylcic inductance
% is calculated based on [1].
%
% [1] Estanislao Juan Pablo Echenique Subiabre, "Improving the performance
% of a wind energy system", PhD thesis, The University of Edinburgh, 2014,
% pp 108-109.
%
% Input
%
%  design - structure which must contain the following fields:
%
%   PhaseResistance : containing the scalar value of the design phase
%    resistance.
%
%   PhaseInductance : vector of one or two values. The first value must be
%    the machine pahse self-inductance. If the optional second value is
%    present, this is the mutual inductance between phases. If a secon
%    value is not present the mutual inductance is assumet to be zero.
%
%  omega - electrical frequency in rad/s at which the phase impedance is to
%   be calculated. Can be a matrix of any size to calculate the impedance
%   at mutliple freqencies.
%
% Output
%
%  Z - complex valued impedance of the phase impdeance seen by an external
%   three phase load at the frequencies given in omega.
%
%
%
% See Also: cyclicinductance_AM
%

    check.isNumeric (omega, true, 'omega', 1);
    
    Ls = cyclicinductance_AM (design);
        
    Z = design.PhaseResistance + 1j .* omega .* Ls;

end

