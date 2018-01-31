function [fens,gcells,maglabels,shaftlabels,outersep,simoptions] = rotormodulesmesh_AF(design, simoptions)
% creates a FAESOR finite element mesh of one or more rotor modules of an
% axial flux machine
%
% Syntax
%
% [fens,gcells,maglabels,shaftlabels,outersep,simoptions] = rotormodulesmesh_AF(design, simoptions)
%
% 

    % set the subset of rotor modules to be evaluates to 1 if not 
    simoptions.Evaluation.structmeshoptions = ...
        setfieldifabsent(simoptions.Evaluation.structmeshoptions, 'ModuleSubset', 1);
    
    maglabels = 1:design.NStages+1;

    tdiscsep = (2*design.g) + design.tc;
    
    maxDiscAxialLayers = 3;
    % limit the number of layers in the mesh
    if any(simoptions.Evaluation.structmeshoptions.DiscAxialLayersPerM * design.tbi > maxDiscAxialLayers)
        simoptions.Evaluation.structmeshoptions.DiscAxialLayersPerM = (1 / max(design.tbi)) * maxDiscAxialLayers; 
    end
    
    maxShaftAxialLayers = 3;
    % limit the number of layers in the mesh
    if simoptions.Evaluation.structmeshoptions.ShaftAxialLayersPerM * tdiscsep > maxShaftAxialLayers
        simoptions.Evaluation.structmeshoptions.ShaftAxialLayersPerM = (1 / tdiscsep) * maxShaftAxialLayers; 
    end
    
    maxSupportAxialLayers = 3;
    % limit the number of layers in the mesh
    if simoptions.Evaluation.structmeshoptions.SupportAxialLayersPerM * design.tsuppb > maxSupportAxialLayers
        simoptions.Evaluation.structmeshoptions.SupportAxialLayersPerM = (1 / design.tsuppb) * maxSupportAxialLayers; 
    end
    
    % ensure not too many circumferential points are present in the mesh
    maxCircumPointsFraction = 1 / 20;
    if 1 / (simoptions.Evaluation.structmeshoptions.CircumPointsPerM * (pi * 2 * design.Rbo)) < maxCircumPointsFraction;
        simoptions.Evaluation.structmeshoptions.CircumPointsPerM = 1 / (maxCircumPointsFraction * pi * 2 * design.Rbo);
    end
    
    % also make sure a minimum number of circumferential points are present
    % in the mesh
    minModuleCircumPoints = 5;
    if simoptions.Evaluation.structmeshoptions.CircumPointsPerM * (pi * 2 * design.Rmo * simoptions.Evaluation.structmeshoptions.ModuleSubset) / design.NModules < minModuleCircumPoints
        simoptions.Evaluation.structmeshoptions.CircumPointsPerM = minModuleCircumPoints / ((pi * 2 * design.Rmo * simoptions.Evaluation.structmeshoptions.ModuleSubset) / design.NModules);
    end
    
    maxBackIronRadialPoints = 3;
    if simoptions.Evaluation.structmeshoptions.BackIronRadialPointsPerM * (design.Rmi - design.Rbi) > maxBackIronRadialPoints
        simoptions.Evaluation.structmeshoptions.BackIronRadialPointsPerM = (1 / (design.Rmi - design.Rbi)) * maxBackIronRadialPoints;
    end
    
    maxMagnetRadialPoints = 3;
    if simoptions.Evaluation.structmeshoptions.MagnetRadialPointsPerM * (design.Rmo - design.Rmi) > maxMagnetRadialPoints
        simoptions.Evaluation.structmeshoptions.MagnetRadialPointsPerM = (1 / (design.Rmo - design.Rmi)) * maxMagnetRadialPoints;
    end
    
    [fens,gcells,maglabels,shaftlabels,outersep] = ...
         axialrotormesh(design.Rs, ...
                        design.Rbi, ...
                        design.Rmi, ...
                        design.Rmo, ...
                        design.Rbo, ...
                        design.tbi, ...
                        design.tsuppb, ...
                        tdiscsep, ...
                        design.tausupp, ...
                        design.NModules, ...
                        design.NModuleSupports, ...
                        'ShaftAxialLayersPerM', simoptions.Evaluation.structmeshoptions.ShaftAxialLayersPerM, ...
                        'DiscAxialLayersPerM', simoptions.Evaluation.structmeshoptions.DiscAxialLayersPerM, ...
                        'SupportAxialLayersPerM', simoptions.Evaluation.structmeshoptions.SupportAxialLayersPerM, ...
                        'CircumPointsPerM', simoptions.Evaluation.structmeshoptions.CircumPointsPerM, ...
                        'MagnetLabels', maglabels, ...
                        'NStages', design.NStages, ...
                        'BackIronRadialPointsPerM', simoptions.Evaluation.structmeshoptions.BackIronRadialPointsPerM, ...
                        'MagnetRadialPointsPerM', simoptions.Evaluation.structmeshoptions.MagnetRadialPointsPerM, ...
                        'Subset', simoptions.Evaluation.structmeshoptions.ModuleSubset, ...
                        'MeshType', simoptions.Evaluation.structmeshoptions.MeshType);
                    
end