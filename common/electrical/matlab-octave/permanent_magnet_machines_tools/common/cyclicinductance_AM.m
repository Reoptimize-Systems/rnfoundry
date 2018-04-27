function Ls = cyclicinductance_AM (design)
% calculates the machine inductance seen by an external 3 phase load
%
% Syntax
%
% Z = phaseimpedance_AM (design, omega)
%
% Description
%
% calculates the machine cyclic inductance seen by an external three phase
% load connected to the machine terminals. Taken from [1].
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
% Output
%
%  Ls - cyclic phase inductance seen by an external three phase load,
%   taking into account mutual inductance.
%
%
% See Also: phaseimpedance_AM
%

    Laa = design.PhaseInductance(1);
    
    if numel (design.PhaseInductance) > 1
        M = design.PhaseInductance(2);
    else
        M = 0;
    end

    Lm = M / -0.5;

    Lsigma = Laa - Lm;

    Ls = Lsigma + 1.5 * Lm;
    
end