function [force, stress] = sheathforcefromfemm_ACTIAM (design, displ)
% determines the force on the sheath of an ACTIAM using FEMM
%
% Syntax
%
% [force, stress] = sheathforcefromfemm_ACTIAM(design)
%
% 

    design.FemmProblem = femmprob_ACTIAM(design);

    % analyse
    [ansfilename, femfilename] = analyse_mfemm (design.FemmProblem, ...
                                                false, ...
                                                true);
    % load solution
    solution = fpproc(ansfilename);
    
    delete (femfilename);
    
    %[sheathForce, Area, NormSheathForce] = FEMMGetSheathForce(design.WpVRm, design.RoVRm, design.Rm, 'm');
    
    % Extract the total energy stored in the magnetic field
    fEnergy = [ solution.totalfieldenergy(), solution.totalfieldcoenergy()];
    
    % Reduce the air gap and perform the sim again
    if nargin < 2
        displ = design.Ro * 0.01;
    end
    
    design.Ro = design.Ro - displ;
    design.Ra = design.Ra - displ;
    
    design = dimensions2ratios_ACTIAM(design);
    
    design.FemmProblem = femmprob_ACTIAM(design);

    % analyse
    [ansfilename, femfilename] = analyse_mfemm (design.FemmProblem, ...
                                                false, ...
                                                true);
    % load solution
    solution = fpproc(ansfilename);
    
    delete (femfilename);
    
    fEnergy(1,3:4) = [ solution.totalfieldenergy(), solution.totalfieldcoenergy()];
    
    force = (fEnergy(1,4)-fEnergy(1,2)) ./ displ;
    
    Area = 2 * pi() * design.Ro * design.Wp;
    
    stress = force / Area;

end