function [ Ptot, Pi, Pd, Pon, Poff, Prr ]= ...
    igbtpwminverterlossses(I_CN, V_CEN, V_CE0, V_FN, V_F0, t_rN, t_fN, t_rrN, Q_rrN, I_CM, V_cc, F_s, M, pf, mtype)
% calculates losses in the IGBT switches of PWM inverters
%
% Syntax 
%
% [ Ptot, Pi, Pd, Pon, Poff, Prr ]= ...
%       igbtpwmlossses(I_CN,V_CEN,V_CE0,V_FN,V_F0,t_rN,t_fN,t_rrN,Q_rrN,I_CM,V_cc,F_s,M)
%
% Description
%
% calculates the losses in a single switch of a PWM inverter employing
% IGBTs from information usually available in manufacturors catalogs. The
% system will be able to calculate losses for different modulation methods,
% provided the current output is sinusoidal. This function implements the
% methods described in [1].
%
% Note that the results here are defined for 2nd generation 1200V IGBTs,
% and depend on the
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
%   pf - desired power factor
%
%   mtype - (optional) scalar specifying modulation type, can be 0 or 1, if
%     0 sine modulation is used, if 1, sine modulation with the third
%     harmonic is used. Default is 1 if not supplied. An mtype of 1 can be
%     used for other modulation systems which also give 100% output
%     voltage, such as bus clamping and vector modulation with good
%     accuracy.
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

% Created by Richard Crozier 2013

    if nargin < 15
        mtype = 1;
    end
    
    if mtype == 0
        % modulation is simple sine wave
        
        % calculate the IGBT forward losses
        Pi = ( 0.125 + M./(3*pi)) .* ( (V_CEN - V_CE0) ./ I_CN ) .* I_CM.^2 ...
             + (1/(2*pi) + 0.125.*M.*pf) .* V_CE0 .* I_CM;

        % the diode forward losses
        Pd = ( 0.125 - M./(3*pi)) .* ( (V_FN - V_F0) ./ I_CN) .* I_CM.^2 ...
             + (1/(2*pi) - 0.125.*M.*pf) .* V_CE0 .* I_CM;        

    elseif mtype == 1
        % modulation is sine wave with third harmonic
        
        % do some subcalculations for speed
        a = 2 * sqrt(3) ./ (9 * pi);
        b = sqrt(3) / (45 * pi);
        c = sqrt(3) / 12;
        theta = acos(pf);
        
        % calculate the IGBT forward losses
        Pi = ( 0.125 + a.*M.*pf - b.*M.*cos(3.*theta) ) .* ( (V_CEN - V_CE0) ./ I_CN ) .* I_CM.^2 ...
             + (1/(2*pi) + c.*M.*pf) .* V_CE0 .* I_CM;

        % the diode forward losses
        Pd = ( 0.125 - a.*M.*pf + b.*M.*cos(3.*theta) ) .* ( (V_FN - V_F0) ./ I_CN) .* I_CM.^2 ...
             + (1/(2*pi) - c.*M.*pf) .* V_CE0 .* I_CM;
         
    end
    
    % IGBT turn-on losses
    % 
    % t_rN = rated rise time
    % V_cc = DC bus voltage
    % F_s = switching frequency
    % I_CM = maximum collector current
    % I_CN = rated collector current
    Pon = 0.125 .* V_cc .* t_rN .* (I_CM.^2 ./ I_CN) .* F_s;
    
    % IGBT turn-off losses
    % 
    % t_fN = rated fall time
    % V_cc = DC bus voltage
    % F_s = switching frequency
    % I_CM = maximum collector current
    % I_CN = rated collector current
    Poff = V_cc .* I_CM .* t_fN .* F_s .* ( (1/(3*pi)) + (1/24).*(I_CM./I_CN));
    
    % Recovery losses
    %
    % Q_rrN = rated recovery charge
    % t_rrN = rated recovery time
    %
    Prr = F_s .* V_cc .* ( (0.28 + (0.38.*I_CM)./(pi.*I_CN) + 0.015 .*(I_CM./I_CN).^2) .* Q_rrN ...
                          + ( 0.8/pi + 0.05.*(I_CM./I_CN) ).*I_CM.*t_rrN  );
    
    
	% get the total power loss
    Ptot = Pi + Pd + Pon + Poff + Prr;
    
end

