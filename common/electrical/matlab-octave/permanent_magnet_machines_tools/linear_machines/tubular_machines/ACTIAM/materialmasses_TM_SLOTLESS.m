function [design, simoptions] = materialmasses_TM_SLOTLESS (design, simoptions)

    % Calculate the mass of the laminated iron in the armature
    if ~isfield (design, 'ArmatureIronVol')
        design.ArmatureIronVol = design.Poles(2) ...
                                  * design.zp ...
                                  * pi ...
                                  * (design.Ryi.^2 - design.Ryo.^2);
    end
    
    [design, simoptions] = materialmasses_TM (design, simoptions);
    
end