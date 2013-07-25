function [ Ptot, Pi, Pd, Pon, Poff, Prr ]= igbtpwminverterlossses(I_CN,V_CEN,V_CE0,V_FN,V_F0,t_rN,t_fN,t_rrN,Q_rrN,I_CM,V_cc,F_s,M)
% calculates losses in the IGBT switches of PWM inverters
%
% Syntax 
%
% [ Ptot, Pi, Pd, Pon, Poff, Prr ]= igbtpwmlossses()
%
% Description
%
% calculates the losses in a single switch of a PWM inverter employing
% IGBTs from information usually available in manufacturors catalogs. The
% system will be able to calculate losses for different modulation methods,
% provided the current output is sinusoidal. This function implements the
% methods described in [1].
% 
% Input
%
%   All inputs may be vectors or matrices of the same size, the outputs
%   will then be the same size.
%
%   I_CN - rated collector current
%
%   V_CEN - rated collector-to-emitter voltage
%
%   V_CE0 - threshold V_CE
%
%   V_FN - rated diode forward voltage
%
%   V_F0 - diode threshold voltage
%
%   t_rN - rated rise time
%
%   t_fN - rated fall time
%
%   t_rrN - rated recovery time
%
%   Q_rrN - rated recovery charge
%
%   I_CM - maximum collector current
%
%   V_cc - DC bus voltage
%
%   F_s - switching frequency
%
% Output
%
%   Ptot - Total power loss in a single switch
% 
%   Pi - The IGBT forward conduction power loss
% 
%   Pd - The diode conduction power loss
% 
%   Pon - The IGBT gate turn-on power loss
% 
%   Poff - The IGBT gate turn-off power loss
% 
%   Prr - The IGBT recovery loss
%
% [1] Casanellas, F., "Losses in PWM inverters using IGBTs",
% IEE-Proc.-Elecrr. Power Appl., Vol. 141, No. 5, September 1994
%
%

% Created by Richard Crozier 2012

    % The typical voltage/current graph V_CE / I_CE is approximated by
    % the following linear equation, where:
    % 
    % I_C = collector current
    % I_CN = rated I_C
    % V_CE = collector-to-emitter voltage
    % V_CN = rated V_CE
    % V_CE0 = threshold V_CE
    
%     V_CE = ( (V_CEN - V_CE0) .* I_C ) / I_CN + V_CE0;
%     
%     % The diode forward voltage follows an exponential law. In the working
%     % range, we may simplify the equation, approximating it to a linear law
%     % with the origin at V_F0. 
%     if fdafd
%         V_F = ( (V_FN - V_F0) * I_C ) / I_CN + V_F0;
%     else
%         % The threshold voltage is taken as 0.7 V
%         V_F = 0.7;
%     end
    
    % calculate the IGBT forward losses
    Pi = ( 0.125 + M./(3*pi)) .* ( (V_CEN - V_CE0) ./ I_CN) .* I_CM.^2 ...
         + (1/(2*pi) + 0.125.*M.*cos(theta)) .* V_CE0 .* I_CM;
    
    % the diode forward losses
    Pd = ( 0.125 - M./(3*pi)) .* ( (V_FN - V_F0) ./ I_CN) .* I_CM.^2 ...
         + (1/(2*pi) - 0.125.*M.*cos(theta)) .* V_CE0 .* I_CM;
     
    % IGBT turn-on losses
    % 
    % t_rN = rated rise time
    % V_cc = DC bus voltage
    % F_s = switching frequency
    % I_CM = maximum collector current
    % I_CN = rated collector current
    Pon = 0.125 .* V_cc .* t_rN .* (I_CM.^2 / I_CN) .* F_s;
    
    % IGBT turn-off losses
    % 
    % t_fN = rated fall time
    % V_cc = DC bus voltage
    % F_s = switching frequency
    % I_CM = maximum collector current
    % I_CN = rated collector current
    Poff = V_cc .* I_CM .* t_fN .* F_s .* ( (1/(3*pi)) + (1/24).*(I_CM/I_CN));
    
    % Recovery losses
    %
    % Q_rrN = rated recovery charge
    % t_rrN = rated recovery time
    %
    Prr = F_s .* V_cc .* ( (0.28 + (0.38.*I_CM)./(pi.*I_CN) + 0.015 .*(I_CM./I_CN).^2) .* Q_rrN ...
                          + ( 0.8/pi + 0.05.*(I_CM./ICN) ).*I_CM.*t_rrN  );
    
    
	% get the total power loss
    Ptot = Pi + Pd + Pon + Poff + Prr;
    
end

