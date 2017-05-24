function design = machineodesim_common_AM (design, simoptions, vCoilRelToField)

    % Now modify the resistance matrix to account for temperature
    
    if simoptions.TempDependantResistivity
        % get the temperature dependent resistivity
        rho = tempdepresistivity ( design.WireResistivityBase, ...
                                   design.AlphaResistivity, ...
                                   design.TemperatureBase, ...
                                   simoptions.Temperature );
                               
        % modify the resistance matrix to account for temperature
        design.RPhase = design.RDCPhase .* rho ./ design.WireResistivityBase;
    end
    
    if simoptions.FreqDependantResistance
        % calculate the electrical frequency
        fe = abs(vCoilRelToField) ./ (2 .* design.PoleWidth);  %velocity2electricalfreq (vCoilRelToField, design.PoleWidth);

        % Then modify to account for skin depth
        design.RPhase = design.NStrands .* ...
            roundwirefreqdepresistance (design.WireStrandDiameter./2, ...
                                        design.RPhase./design.NStrands, ...
                                        rho, 1, fe);
    end
                                
end

function rho = tempdepresistivity(rhobase, alpha, Tbase, T)
% calculates the temperature dependent resistivity from base values and a
% coefficient
%
% Syntax
%
% rho = tempdepresistivity(rhobase, alpha, Tbase, T)
%
% Input
%
%  rhobase - base value of the resistivity at the temperature supplied in
%    Tbase
%
%  alpha - temperature coefficient of resistivity of the material
%
%  Tbase - temperature at which the material has the resistivity supplied
%    in rhobase
%
%  T - the temperature for which the actual resistivity is to be calculated
%
% Output
%
%  rho - the resistivity at the temperature supplied in T
%
% 

% Copyright Richard Crozier 2012 - 2012

    % Rho is the resistivity of copper at 'design.temperature'. We
    % calculate this based on the known resistivity at a reference
    % simparams.temperature (design.rho20wire). In this case the reference
    % temperature is 20 degrees celcius
    rho = rhobase * (1 + alpha * (T - Tbase));
            
end

function Rac = roundwirefreqdepresistance(Dc, Rdc, rho, mu_r, freq)
% calcuates the AC resitance of a wire of round cross-section due to the
% skin effect
%
% Syntax
%
% Rac = roundwirefreqdepresistance(Dc, Rdc, rho, mu_r, freq)
%
% Input
%
%   Dc - wire diameter
%
%   Rdc - the DC resistance of the wire
% 
%   rho - the resistivity of the wire material
% 
%   mu_r - the relative permeability of the wire material
% 
%   freq - the frequency of the current waveform
%
% Output
%
%   Rac - the AC wire resistance including the skin effect
%
% Description
%
% The AC winding resistance is calculated according to the formulas
% presented in 'The Analysis of Eddy Currents', Richard L Stoll, Chapter 2,
% Section 2.8, page 25
%
% 

% Copyright Richard Crozier 2012 - 2012

    % determine the skin depth
    delta = skindepth (rho, mu_r, freq);
    
    % calculate the AC resistance, this is dependent on the ratio of the
    % wire radius to the skin depth, as described in 'The Analysis of Eddy
    % Currents', Richard L Stoll, Chapter 2, Section 2.8, page 25
    if Dc > 7 * delta
        
        Rac = Rdc .* ( (Dc / (2*delta)) + 0.25 + (3*delta / (32*Dc)) );
        
    else
        
        Rac = Rdc .* ( 1 + realpow (Dc,2) / ( 4 .* realpow (delta,2) ) );
        
    end


end
            