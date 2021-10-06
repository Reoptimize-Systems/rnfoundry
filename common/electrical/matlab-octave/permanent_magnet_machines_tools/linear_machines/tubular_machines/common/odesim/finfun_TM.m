function [design, simoptions] = finfun_TM (design, simoptions)
% finalises common aspects of tubular linear generator models prior to a
% dynamic simulation
%
% Syntax
%
% [design, simoptions] = finfun_TM(design, simoptions)
%

    design = setfieldifabsent (design, 'PowerPoles', design.Poles(1));
    
    % check and set common linear machine parameters
    [design, simoptions] = finfun_linear(design, simoptions);
    
    if simoptions.DoPreLinSim
        
        evaloptions = designandevaloptions_TM(simoptions.Evaluation);
        
        simoptions.FieldIronDensity = evaloptions.FieldIronDensity;
        simoptions.MagnetDensity = evaloptions.MagnetDensity;
        simoptions.ShaftDensity = evaloptions.StructMaterialDensity;
        simoptions.CopperDensity = evaloptions.CopperDensity;
        simoptions.ArmatureIronDensity = evaloptions.ArmatureIronDensity;

        design.PoleWeight = fieldpoleweight_TM(design.WmVWp, ...
            design.WpVRm, design.RsiVRso, design.RsoVRm, ...
            design.Rm, simoptions.FieldIronDensity, simoptions.MagnetDensity, ...
            simoptions.ShaftDensity, design.Rs2VHmag, design.Rs1VHmag, ...
            design.Ws2VhalfWs, design.Ws1VhalfWs);
        
        % armaturepoleweight_TM(WpVRm, RoVRm, Rm, g, WcVWp, cufill, copperDensity,
        % RaVRo, steelDensity)
        design.PoleWeight = design.PoleWeight + armaturepoleweight_TM(design.WmVWp, ...
            design.RoVRm, design.Rm, design.g, design.WcVWp, ...
            design.CoilFillFactor, simoptions.CopperDensity, design.RaVRo, ...
            simoptions.ArmatureIronDensity);
        
    end
    
end
