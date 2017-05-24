function design = test_design_ACTM ()

    % Number of phases in machine
    design.Phases = 3;  
    % Magnet outer radius
    design.Rm = 0.1;
    % air gap
    design.g = 5/1000;
    % Coil inner radius
    design.Ri = design.Rm + design.g;
    % ratio of magnet width to pole width
    design.WmVWp = 0.75;
    % ratio of pole width to magnet outer radius
    design.WpVRm = 0.5;
    % ratio of coil inner radius to magnet outer radius
    design.RiVRm = design.Ri / design.Rm;
    % ratio of coil outer radius to magnet outer radius
    design.RoVRm = 1.2;
    % ratio of coil sheath outer radius to coil outer radius
    design.RaVRo = 1.025;
    % ratio of magnet stack shaft outer radius to magnet outer radius
    design.RsoVRm = 0.1;
    % ratio of magnet stack shaft inner radius to magnet stack shaft outer
    % radius
    design.RsiVRso = 0;
    % ratio of coil axial height to pole height
    design.WcVWp = 1/3;
    % the coil wire configuration
    design.CoilFillFactor = 0.65;
    design.CoilTurns = 50;

    design.mode = 2; 

    % poles on the translator and stator respectively
    design.Poles = [10 30];
    % series/parallel branch configuration in phase winding
    design.CoilsPerBranch = 10;
    design.Branches = 1;
    
    % ratio of load resistance to phase resistance
    design.RlVRp = 10;
    % ratio of load inductance to phase inductance
    design.LgVLc = 0;

end