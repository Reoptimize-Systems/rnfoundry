function Pv = steinmetzloss (f, Bm, kh, kc, ke, beta)
% calculates iron loss per volume using steinmetz loss equation
%
% Syntax
%
% 
% Pv = steinmetzloss (f, Bm, kh, kc, ke, beta)
%
% Inputs
%
%  f - frequency of applied induction
%
%  Bm - magnitude of aplied induction
%
%  kh - hysteresis loss coefficient
%
%  kc - eddy current loss coefficient
%
%  ke - excess loss coefficient
%
%  beta - hysteresis loss exponent
%
% Outputs
%
%  Pv - power loss per unit volume for each combination of inputs
%
%

    Pv = kh .* f .* (Bm.^beta) ...
        + kc .* realpow (f .* Bm, 2) ...
        + ke .* realpow (f .* Bm, 1.5);
    
end