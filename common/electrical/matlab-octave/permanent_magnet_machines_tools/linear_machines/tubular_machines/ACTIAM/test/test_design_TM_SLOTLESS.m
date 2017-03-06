function design = test_design_TM_SLOTLESS (spectype)

    if nargin < 1
        spectype = 'dims';
    end
    
    design.ArmatureType = 'external';
    
    % Number of phases in machine
    design.Phases = 3;
    % poles on the translator and stator respectively
    design.Poles = 10;
    % series/parallel branch configuration in phase winding
    design.CoilsPerBranch = 10;
    design.Branches = 1;
    design.CoilLayers = 1;
    
    design.Qc = design.Phases * design.Poles;
%     if design.CoilLayers == 1
%         design.Qs = 2 * design.Qc;
%     elseif design.CoilLayers == 2
        design.Qs = design.Qc;
%     end
    design.qc = fr (design.Qc, design.Poles * design.Phases);
    design.yd = 1;
    % the coil wire configuration
    design.CoilFillFactor = 0.65;
    design.CoilTurns = 50;
    
    switch spectype
        case 'ratios'
            % RyiVRyo, RtsbVRyi, RagVRtsb, RfoVRag, RsoVRfo, tsgVtsb, zmVzp, zcgVzs, zcyVzs, zsgVzcg

%             % Magnet outer radius
%             design.Rfo = 0.1;
%             % air gap
%             design.g = 5e-3;
%             % ratio of magnet width to pole width
%             design.zmVzp = 0.75;
%             % ratio of pole width to magnet outer radius
%             design.zpVRfo = 0.5;
%             % ratio of coil inner radius to magnet outer radius
%             design.RfoVRag = design.Rfo / design.Rag;
%             % ratio of coil outer radius to magnet outer radius
%             design.RcyVRfo = 1.2;
%             % ratio of coil sheath outer radius to coil outer radius
%             design.RaVRcy = 1.025;
%             % ratio of magnet stack shaft outer radius to magnet outer radius
%             design.RsoVRfo = 0.1;
%             % ratio of magnet stack shaft inner radius to magnet stack shaft outer
%             % radius
%             design.RsiVRso = 0;
%             % ratio of coil axial height to pole height
%             design.zcyVzs = 1/3;
%             design.zcgVzs = 1/3;
    
        case 'dims'
            
            % Ryo, g, ry, rm, rc, rsb, rsg, zm, zcg, zcy, zsg
            design.Ryo = 0.5;
            design.g = 3e-3;
            design.ry = 10e-3;
            design.rm = 120e-3;
            design.rc = 30e-3;
            design.rsb = 10e-3;
            design.rsg = 5e-3;
            design.zp = 0.2;
            design.zm = 0.8 * design.zp;
            zs = design.zp / 3;
            design.zcg = zs * 0.8;
            design.zcy = zs * 0.9;
            design.zsg = design.zcg * 0.5;
            
            
        case 'rdims'
            
            % Ryo, Ryi, Rtsb, Rag, Rso, Rfo, rsg, zm, zcg, zcy, zsg
%             design.Ryo = 
%             design.Ryi = 
%             design.Rtsb = 
%             design.Rag = 
%             design.Rso = 
%             design.Rfo = 
%             design.rsg = 
%             design.zm = 
%             design.zcg = 
%             design.zcy = 
%             design.zsg = 
            
            
    end
    

    design.mode = 2;
    
    
    
    % ratio of load resistance to phase resistance
    design.RlVRp = 10;
    % ratio of load inductance to phase inductance
    design.LgVLc = 0;
    
    design.MagFEASimMaterials.Magnet = 'NdFeB 40 MGOe';
    design.MagFEASimMaterials.FieldSheath = '1117 Steel';
    design.MagFEASimMaterials.MagnetSpacer = '1117 Steel';
    design.MagFEASimMaterials.ArmatureYoke = 'M-27 Steel';
    design.MagFEASimMaterials.ArmatureCoil = '36 AWG';
    design.MagFEASimMaterials.AirGap = 'Air';
    design.MagFEASimMaterials.CoilInsulation = 'Air';
    
    design.HeatFEASimMaterials.Magnet = 'Iron, Pure';
    design.HeatFEASimMaterials.FieldSheath = 'Iron, Pure';
    design.HeatFEASimMaterials.MagnetSpacer = 'Iron, Pure';
    design.HeatFEASimMaterials.ArmatureYoke = 'M-27 Steel';
    design.HeatFEASimMaterials.ArmatureCoil = 'Copper, Pure';
    design.HeatFEASimMaterials.AirGap = 'Water';
    design.HeatFEASimMaterials.CoilInsulation = 'Water';

end