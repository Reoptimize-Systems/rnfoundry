function [flux_linkage, FEF, lossinfo] = feaode_RADIAL_SLOTTED (design, simoptions, theta, Icoils)
% generate fea information for a slotted radial flux machine
%
% Syntax
%
% [flux_linkage, FEF, lossinfo] = feaode_RADIAL_SLOTTED (design, simoptions, theta, Icoils)
%
% Inputs
%
%  design - design structure for slotted radial flux machine
%
%  simoptions - simulation options structure
%
%  theta - angular position of the rotor
%
%  Icoils - Currents in the machine coils
%
% Output
%
%  flux_linkage - (n x 1) vector of flux linkage values in the coil
%    circuits
%
%  FEF - machine torque
%
%  lossinfo - loss information
%
%

    FEF = nan * ones (1, size (Icoils, 2));
    flux_linkage = nan * ones (size (Icoils));
    
    if simoptions.MagFEASim.UseParFor
        
        parfor ind = 1:size (Icoils, 2)
            [ FEF(ind), ...
              BxCoreLossData(:,:,ind), ...
              ByCoreLossData(:,:,ind), ...
              ~, ...
              flux_linkage(:,ind) ] = feasim_RADIAL_SLOTTED (design, simoptions, theta, ...
                                                      'IsInitialisation', false, ...
                                                      'PhaseCurrents', Icoils(:,ind), ...
                                                      'GatherIronLossData', false);
        end
        
    else
        
        for ind = 1:size (Icoils, 2)
            [ FEF(ind), ...
              BxCoreLossData(:,:,ind), ...
              ByCoreLossData(:,:,ind), ...
              ~, ...
              flux_linkage(:,ind) ] = feasim_RADIAL_SLOTTED (design, simoptions, theta, ...
                                                      'IsInitialisation', false, ...
                                                      'PhaseCurrents', Icoils(:,ind), ...
                                                      'GatherIronLossData', false);
        end
        
    end
    
    lossinfo = struct ( 'BxCoreLossData', BxCoreLossData, ...
                        'ByCoreLossData', ByCoreLossData );
    
end