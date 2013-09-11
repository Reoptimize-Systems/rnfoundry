function [g, closingForces, BeamInfo, design] = feairgapclosure_ACPMSM(design, options, BeamInfo)

    if nargin < 3 || isempty(BeamInfo)
        BeamInfo.LastOuterStructureNBeams = NaN;
        BeamInfo.PreviousMeshes = {};
        BeamInfo.PreviousMeshesNBeams = [];
    end
    
    % Calculate the weight acting on the support beams due to the mass of
    % the translator
    OuterPoleWeight = fpoleweight_ACPMSM(design, options);
    
    [g, closingForces, BeamInfo, design] = feairgapclosure_FM(design, options, BeamInfo, OuterPoleWeight, @completebeaminfo_ACPMSM);
    
end