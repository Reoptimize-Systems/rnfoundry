function [design, simoptions] = materialmasses_TM (design, simoptions)
% calculates the masses of common components of tubular linear machines
%
% Syntax
%
% [design, simoptions] = materialmasses_TM(design, simoptions)
%
% 

    [design, simoptions] = materialmasses_arm_TM (design, simoptions);

    [design, simoptions] = materialmasses_field_TM (design, simoptions);
    
    % set the epoxy mass to zero for the moment
    design.EpoxyVol = 0;
    design.EpoxyMass = design.EpoxyVol * simoptions.Evaluation.EpoxyDensity;
    
    % TODO calculate the support mass
    design.ArmatureSupportMass = 0;
    design.FieldSupportMass = 0;
    
    % calculate the component and module masses
    design.FieldMass = design.FieldIronMass ...
                       + design.MagnetMass ...
                       + design.FieldSupportMass;
    
    design.ArmatureMass = design.CopperMass ...
                        + design.EpoxyMass ...
                        + design.ArmatureSupportMass ...
                        + design.ArmatureIronMass;
                    
	design.TotalMass = design.FieldMass + design.ArmatureMass;
    
end

function [design, simoptions] = materialmasses_arm_TM (design, simoptions)
% calculates the masses of the field component of tubular linear machines
%
% Syntax
%
% [design, simoptions] = materialmasses_arm_TM(design, simoptions)
%
% 

    if ~isfield (design, 'MTL')
        design.MTL = 2 * pi * design.Rcm;
    end
    
    design.WireLengthTotal = design.CoilTurns * design.MTL;
    
    % Determine the copper cost based on the length of wire in a coil
    % multiplied by the number of Phases and armature Poles and the wire
    % cross-sectional area etc.
    design.CopperVol = design.Qc ...
                        * design.WireLengthTotal ...
                        * (pi * (design.Dc/2)^2);
    
    if (numel (design.Poles)) > 1 && (design.Poles(2) > design.Poles(1))
        % scale up by the amount the armature overlaps field in terms of
        % number of poles
        design.CopperVol = design.CopperVol * (design.Poles(2) / design.Poles(1));
        
        design.ArmatureIronVol = design.ArmatureIronVol ...
                                  * (design.Poles(2) / design.Poles(1));
    end
    
    design.CopperMass = design.CopperVol ...
                        * simoptions.Evaluation.CopperDensity;
    
    design.ArmatureIronMass = design.ArmatureIronVol ...
                              * simoptions.Evaluation.ArmatureIronDensity;

end

function [design, simoptions] = materialmasses_field_TM (design, simoptions)
% calculates the masses of the field component of tubular linear machines
%
% Syntax
%
% [design, simoptions] = materialmasses_field_TM(design, simoptions)
%
% 

    if ~isfield (design, 'MagnetVol')
        
        design.MagnetVol = design.Poles(1) ...
                            * design.Wm ...
                            * pi * (design.Rm.^2 - design.Rso^2);

    end
    
    % Now calculate the total mass of the magnets
    design.MagnetMass = design.MagnetVol ...
                        * simoptions.Evaluation.MagnetDensity;

    if ~isfield (design, 'FieldSpacerVol')

        % Calculate the cost of the back iron on the field, i.e. the steel
        % discs in the translator
        if design.mode == 1 || design.mode == 3
            cavityvol = steelcavityvol_TM (design);
        else
            cavityvol = 0;
        end

        design.FieldSpacerVol = (design.zs * pi * (design.Rfo^2 - design.Rso^2) - cavityvol) ...
            * design.Poles(1);

    end
    
    design.FieldIronMass = design.FieldSpacerVol * simoptions.Evaluation.FieldIronDensity;
    
end

