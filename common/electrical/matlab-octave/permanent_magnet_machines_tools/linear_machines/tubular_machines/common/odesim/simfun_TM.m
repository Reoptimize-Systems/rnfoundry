function [design, simoptions] = simfun_TM (design, simoptions)
% simulation preprocessing function common to all tubular linear machines
%
% Syntax
%
% [design, simoptions] = simfun_TM(design, simoptions)
%
% Input
%
% design - a tubular machine design structure 
% 

    % for a tubular machine the MTL is the radial distance to the coil
    % 2D shape centroid, use it if present, othersie use middle of coil
    if isfield (design, 'CoilCentroid')
        design.MTL = design.CoilCentroid(1);
    else
        % if not avaialble use the 
        design.MTL = 2 * pi * design.Rcm;
    end
                          
	% do common linear machine setup tasks
    [design, simoptions] = simfun_linear(design, simoptions);
    
    % set some default mesh size options
    simoptions.MagFEASim = setfieldifabsent (simoptions.MagFEASim, ...
        'MagnetRegionMeshSize', choosemesharea_mfemm(design.rm, design.zm, 1/10));
    
    simoptions.MagFEASim = setfieldifabsent (simoptions.MagFEASim, ...
        'MagnetSpacerRegionMeshSize', choosemesharea_mfemm(design.rm, design.zsd, 1/10));
    
    simoptions.MagFEASim = setfieldifabsent (simoptions.MagFEASim, ...
        'AirGapMeshSize', choosemesharea_mfemm(design.g, design.zp, 1/10));
    
    simoptions.MagFEASim = setfieldifabsent (simoptions.MagFEASim, ...
        'OuterRegionsMeshSize', [choosemesharea_mfemm(design.rm, design.zp, 1/5), -1]);
    
    simoptions.MagFEASim = setfieldifabsent (simoptions.MagFEASim, 'YokeRegionMeshSize', ...
                                       mean( [choosemesharea_mfemm(design.ry, 2*design.zp, 1/10), ...
                                        choosemesharea_mfemm(design.rc(1), design.zs-max(design.zc), 1/10)] ) );
                                    
    simoptions.MagFEASim = setfieldifabsent (simoptions.MagFEASim, 'CoilRegionMeshSize', ...
                                    choosemesharea_mfemm(design.rc(1), mean(design.zc)) );
    
end