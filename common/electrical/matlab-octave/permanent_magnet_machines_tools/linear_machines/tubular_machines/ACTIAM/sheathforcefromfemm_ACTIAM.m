function [force, stress] = sheathforcefromfemm_ACTIAM(design)
% determines the force on the sheath of an ACTIAM using FEMM
%
% Syntax
%
% [force, stress] = sheathforcefromfemm_ACTIAM(design)
%
% 

    RunStructFEMMSimNew_ACTIAM(design);
    
    mi_analyse(1); mi_loadsolution;
    
    %[sheathForce, Area, NormSheathForce] = FEMMGetSheathForce(design.WpVRm, design.RoVRm, design.Rm, 'm');
    
    % Extract the total energy stored in the magnetic field
    fEnergy = mo_TotalFieldEnergy;
    
    mo_close; mi_close;
    
    % Reduce the air gap and perform the sim again
    displ = design.Ro * 0.01;
    
    design.Ro = design.Ro - displ;
    design.Ra = design.Ra - displ;
    
    design = dimensions2ratios_ACTIAM(design);
    
    RunStructFEMMSimNew_ACTIAM(design);
    
    mi_analyse(1); mi_loadsolution;
    
    fEnergy(1,3:4) = mo_TotalFieldEnergy;
    
    mo_close; mi_close;
    
    force = (fEnergy(1,4)-fEnergy(1,2)) ./ displ;
    
    stress = force / Area;

end