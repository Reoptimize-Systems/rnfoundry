function design = costestimate_AM(design, simoptions)
% calculates the cost of materials present in a design structure, and the
% total cost based on these costs
%
% Syntax
% 
% design = costestimate_AM(design, simoptions)
%
% Description
%
%  simoptions - structure containng unit cost information for material in a
%    in the field 'Evaluation'. 'Evaluation' is itself a structure
%    containing fields with cost information on materials used. What fields
%    must be supplied in the evaloptions structure depends on the contents
%    of the design structure.
%
%  design - a structure containing the following optional fields. Depeding
%   on which are supplied new fields containing cost information are added
%   to the structure.
%
%   CostEstimate - starting cost for the design to which additional
%     material costs will be added. If not supplied, the cost before
%     material costs are added is assumed to be zero.
%
%   CopperMass - total mass of copper in the design. If present the
%     evaloptions structure must contain the field 'CopperCost', which is
%     the cost per kg of copper. The field 'CopperCost' will be added to
%     the design structure containing the total cost of the copper used in
%     the design.
%
%   MagnetMass - total mass of magnet material in the design. If present
%     the evaloptions structure must contain the field 'MagnetCost', which
%     is the cost per kg of the magnets. The field 'MagnetCost' will be added to
%     the design structure containing the total cost of the magnets used in
%     the design.
%
%   FieldIronMass - total mass of field iron material in the design. If
%     present the evaloptions structure must contain the field
%     'FieldIronCost', which is the cost per kg of the material. The field
%     'FieldIronCost' will be added to the design structure containing the
%     total cost of the field iron used in the design.
%
%   ArmatureIronMass - total mass of armature iron material in the design.
%     If present the evaloptions structure must contain the field
%     'ArmatureIronCost', which is the cost per kg of the material. The
%     field 'ArmatureIronCost' will be added to the design structure
%     containing the total cost of the armature iron used in the design.
%
%   StructMaterialMass - total mass of structural material in the design.
%     If present the evaloptions structure must contain the field
%     'StructMaterialCost', which is the cost per kg of the material. The
%     field 'StructMaterialCost' will be added to the design structure
%     containing the total cost of the structural material used in the
%     design.
%
%   EpoxyMass - total mass of epoxy in the design. If present the
%     evaloptions structure must contain the field 'EpoxyCost', which is
%     the cost per kg of the epoxy. The field 'EpoxyCost' will be added to
%     the design structure containing the total cost of the epoxy used in
%     the design.
%

    if ~isfield(design, 'CostEstimate')
        design.CostEstimate = 0;    
    end
    
    if isfield(design, 'CopperMass')
        design.CopperCost = design.CopperMass * simoptions.Evaluation.CopperCost;
        design.CostEstimate = design.CostEstimate + design.CopperCost;
    end
    
    if isfield(design, 'MagnetMass')
        design.MagnetCost = design.MagnetMass * simoptions.Evaluation.MagnetCost;
        design.CostEstimate = design.CostEstimate + design.MagnetCost;
    end
    
    if isfield(design, 'StructMaterialMass')
        design.StructuralCost = design.StructMaterialMass * simoptions.Evaluation.StructMaterialCost;
        design.CostEstimate = design.CostEstimate + design.StructuralCost;
    end

    if isfield(design, 'FieldIronMass')
        design.FieldIronCost = design.FieldIronMass * simoptions.Evaluation.FieldIronCost;
        design.CostEstimate = design.CostEstimate + design.FieldIronCost;
    end
    
    if isfield(design, 'ArmatureIronMass')
        design.ArmatureIronCost = design.ArmatureIronMass * simoptions.Evaluation.ArmatureIronCost;
        design.CostEstimate = design.CostEstimate + design.ArmatureIronCost;
    end
    
    if isfield(design, 'EpoxyMass')
       design.EpoxyCost = design.EpoxyMass * simoptions.Evaluation.EpoxyCost;
       design.CostEstimate = design.CostEstimate + design.EpoxyCost;
    end
    
end