function Pr = powerConverterRating(Vppeak, Ippeak, omegaPeak, Lp)
% powerConverterRating: determines the rating for the power converter given
% a peak phase voltage, phase current and machine inductance
%
% Input:
%
%   Vppeak - peak phase voltage of the machine
%
%   Ippeak - peak phase current of the machine
%
%   omegaPeak - peak electrical frequency of the machine
%
%   Lp - the machine phase inductance
%
% Output:
%
%   Pr - RELATIVE Power rating for the back-to-back AC/AC power converter
%

    Pr = 2 .* sqrt(1 + (omegaPeak .* Lp .* Ippeak ./ Vppeak)^2);
    
end