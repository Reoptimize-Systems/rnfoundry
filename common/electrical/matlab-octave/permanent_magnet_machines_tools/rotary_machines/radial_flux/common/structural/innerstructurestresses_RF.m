function [Kmat, Fmat] = innerstructurestresses_RF(faeprob, Rsoo, ls, shearforce, radialforce, structdensity, omega, labels)
% applies forces to the finite element block representing the internal
% structure of a radial flux machine
%
% Syntax
%
% [Kmat, Fmat] = inerstructurestresses_RF(fens, gcells, feb, geom, u, ...
%                   Rsoo, ls, shearforce, radialforce, structdensity, omega, labels)
%
% 

    
    % Assemble the system matrix
    ems = stiffness(faeprob.feb, faeprob.geom, faeprob.u);
    K = start (sparse_sysmat, get(faeprob.u, 'neqns'));
    K = assemble (K, ems);

    % Loads
    
    % initialise the force matrix
    F = start (sysvec, get(faeprob.u, 'neqns'));

    outersurfgcells = getoutersurface(faeprob, labels.outersurface, Rsoo, ls);
    
    % the shear stress on the surface
    shearfeb = feblock_defor_ss(struct ('mater', getmater(faeprob.feb), ...
                                        'gcells', outersurfgcells,...
                                        'integration_rule', gauss_rule(2, 2)));

    % generate the shear force intensity object
    shearfi = force_intensity( struct ('magn', @(xyz) circtangentforcevec(xyz, shearforce)') );
    
    F = assemble (F, distrib_loads(shearfeb, faeprob.geom, faeprob.u, shearfi, 2));
    
    % the radial stress on the surface
    radfeb = feblock_defor_ss(struct ('mater', getmater(faeprob.feb), ...
                                        'gcells', outersurfgcells,...
                                        'integration_rule', gauss_rule(2, 2)));

    % generate the radial force intensity object
    radfi = force_intensity( struct ('magn', @(xyz) radialforcevec(xyz, radialforce)') );
    
    F = assemble (F, distrib_loads(radfeb, faeprob.geom, faeprob.u, radfi, 2));
    
%                                   
%     %   then the axial forces intensity object
%     axialfi = force_intensity( struct ('magn', [0, 0, radialforces(1)]) );                           
%     
%     F = assemble (F, distrib_loads(magfeb, geom, u, axialfi, 2));
%     
%     [ botmagsurfgcells ] = getmagsurface(fens, gcells, maglabels(end), Rmo);
%     
%     magfeb = feblock_defor_ss(struct ('mater', getmater(feb), ...
%                                       'gcells', botmagsurfgcells,...
%                                       'integration_rule', gauss_rule(2, 2)));
%                                   
%     %   then the axial forces intensity object
%     axialfi = force_intensity( struct ('magn', [0, 0, radialforces(end)]) );                           
%     
%     F = assemble (F, distrib_loads(magfeb, geom, u, axialfi, 2));
%      
%     % do the end rotor forces 
%     magregionvol = (pi * (realpow(Rmo, 2) - realpow(Rmi, 2)) * tbi(1));
%     
%     % calculate the force per unit volume in the magnet region
%     forcemag = shearforce(1) / magregionvol;
%     
%     % generate the shear force intensity object
%     shearfi = force_intensity( struct ('magn', @(xyz) circtangentforcevec(xyz, forcemag)') );
% %     % then the axial forces intensity object
% %     axialfi = force_intensity( struct ('magn', [0, 0, axialforces(1)/magregionvol]) );
%     
%     % shear forces
%     magcelllist = gcell_select(fens, gcells, struct ('label', maglabels(1)) );
%     
%     magfeb = feblock_defor_ss(struct ('mater', getmater(feb), ...
%                                       'gcells', subset(gcells,magcelllist),...
%                                       'integration_rule', gauss_rule(3, 2)));
% 
%     F = assemble (F, distrib_loads(magfeb, geom, u, shearfi, 3));
% %     F = assemble (F, distrib_loads(magfeb, geom, u, axialfi, 3));
% %     
%     F = applymagcentrifugalload(F, fens, gcells, feb, geom, u, maglabels(1), [0, magnetmass], Rmo, Rmi, nmags, omega);
% %                                   
%     % other end
%     % calculate the force per unit volume in the magnet region
%     forcemag = shearforce(end) / magregionvol;
% %     
%     % first generate the shear force intensity object
%     shearfi = force_intensity( struct ('magn', @(xyz) circtangentforcevec(xyz, forcemag)') );
% %     % then the axial forces intensity object
% %     axialfi = force_intensity( struct ('magn', [0, 0, axialforces(end)/magregionvol]) );
% %     
%     magcelllist = gcell_select(fens, gcells, struct ('label', maglabels(end)) );
% 
%     magfeb = feblock_defor_ss(struct ('mater', getmater(feb), ...
%                                       'gcells', subset(gcells,magcelllist),...
%                                       'integration_rule', gauss_rule(3, 2)));
% 
%     F = assemble (F, distrib_loads(magfeb, geom, u, shearfi, 3));
% %     F = assemble (F, distrib_loads(magfeb, geom, u, axialfi, 3));
% %     
%     F = applymagcentrifugalload(F, fens, gcells, feb, geom, u, maglabels(end), [magnetmass, 0], Rmo, Rmi, nmags, omega);
%     
%     % now do the rest of the forces
%     magregionvol = (pi * (realpow(Rmo, 2) - realpow(Rmi, 2)) * tbi(2));
%     
%     for i = 2:numel(maglabels)-1
%     
%         if radialforces(i)  > 0
%             
%             [ botmagsurfgcells, topmagsurfgcells ] = getmagsurface(fens, gcells, maglabels(1), Rmo);
%             
%             magfeb = feblock_defor_ss(struct ('mater', getmater(feb), ...
%                                               'gcells', topmagsurfgcells,...
%                                               'integration_rule', gauss_rule(2, 2)));
%             
%             %   then the axial forces intensity object
%             axialfi = force_intensity( struct ('magn', [0, 0, radialforces(i)]) );
%             
%             F = assemble (F, distrib_loads(magfeb, geom, u, axialfi, 2));
%             
%         elseif radialforces < 0
%             
%             [ botmagsurfgcells ] = getmagsurface(fens, gcells, maglabels(end), Rmo);
%             
%             magfeb = feblock_defor_ss(struct ('mater', getmater(feb), ...
%                                               'gcells', botmagsurfgcells,...
%                                               'integration_rule', gauss_rule(2, 2)));
%             
%             %   then the axial forces intensity object
%             axialfi = force_intensity( struct ('magn', [0, 0, radialforces(i)]) );
%             
%             F = assemble (F, distrib_loads(magfeb, geom, u, axialfi, 2));
%             
% %             
% %             % calculate the force per unit volume in the magnet region
% %             forcemag = shearforce(i) / magregionvol;
% % 
% %             % first generate the shear force intetensity object
% %             shearfi = force_intensity( struct ('magn', @(xyz) circtangentforcevec(xyz, forcemag)') );
% %             % then the axial forces intensity object
% %             axialfi = force_intensity( struct ('magn', [0, 0, axialforces(i)/magregionvol]) );
% % 
% %             magcelllist = gcell_select(fens, gcells, struct ('label', maglabels(end)) );
% % 
% %             magfeb = feblock_defor_ss(struct ('mater', getmater(feb), ...
% %                                               'gcells', subset(gcells,magcelllist),...
% %                                               'integration_rule', gauss_rule(3, 2)));
% % 
% %             F = assemble (F, distrib_loads(magfeb, geom, u, shearfi, 3));
% %             F = assemble (F, distrib_loads(magfeb, geom, u, axialfi, 3));
% %             
%         end
%         
%         F = applymagcentrifugalload(F, fens, gcells, feb, geom, u, maglabels(i), [magnetmass, magnetmass], Rmo, Rmi, nmags, omega);
%         
%     end
    
    % apply radial cetrifugal forces to the entire rotor structure
    % the radial body force applied everywhere in the structure
    radialfi = force_intensity( struct ('magn', @(xyz) centrifugalbodyforce(xyz, omega, structdensity)') );
    
    F = assemble(F, distrib_loads(faeprob.feb, faeprob.geom, faeprob.u, radialfi, 3));
    
    % apply gravity in x direction
    gravityfi = force_intensity( struct ('magn', [1;0;0] * structdensity * 9.81) );
    
    F = assemble(F, distrib_loads(faeprob.feb, faeprob.geom, faeprob.u, gravityfi, 3));
    
    F = finish(F);

    Kmat = get(K,'mat');
    Fmat = get(F,'vec');

end


function [F] = applymagcentrifugalload(F, faeprob, maglabel, magnetmass, Rmo, Rmi, Nm, omega)

    if any(magnetmass > 0) 
        [ botmagsurfgcells, topmagsurfgcells ] = getmagsurface(faeprob.fens, faeprob.gcells, maglabel, Rmo);
    end
    
    % calculate the force per m^2 on the plate face due to the magnet
    % centrifugal force, at the mean magnet radius
    magcentrifugalfpermsq = magnetmass * omega * Nm * (Rmo + Rmi) / (2 * pi * (Rmo^2 - Rmi^2));
    
    if magnetmass(1) > 0
        
        % the magnet suface centrifugal force intensity object
        magradfi = force_intensity( struct ('magn', @(xyz) radialforcevec(xyz, magcentrifugalfpermsq(1))') );

        % apply the centrifugal surface loads due to the magnets
        magfeb = feblock_defor_ss(struct ('mater', getmater(faeprob.feb), ...
                                          'gcells', botmagsurfgcells,...
                                          'integration_rule', gauss_rule(2, 2)));

        F = assemble (F, distrib_loads(magfeb, faeprob.geom, faeprob.u, magradfi, 2));

    end
    
    if magnetmass(2) > 0
        
        % the magnet suface centrifugal force intensity object
        magradfi = force_intensity( struct ('magn', @(xyz) radialforcevec(xyz, magcentrifugalfpermsq(2))') );
        
        % apply the centrifugal surface loads due to the magnets
        magfeb = feblock_defor_ss(struct ('mater', getmater(faeprob.feb), ...
                                          'gcells', topmagsurfgcells,...
                                          'integration_rule', gauss_rule(2, 2)));

        F = assemble (F, distrib_loads(magfeb, faeprob.geom, faeprob.u, magradfi, 2));
        
    end
    
end

function [ outersurfgcells ] = getoutersurface(faeprob, outerlabel, Rsoo, ls)

    % first get the 3D gcells making up the outer surface of the inner
    % structure
    outercells = subset( faeprob.gcells, gcell_select(faeprob.fens, faeprob.gcells, struct('label', outerlabel)) );
    
    % then get the boundary cells of this region
    bdrygcells = mesh_bdry(outercells);

    tol = 10 * eps;

    % then extract only parts of the surface of interest
    outersurfgcells = subset(bdrygcells, ...
                             gcell_select(faeprob.fens, bdrygcells, ...
                                struct('cylinder', [0, 0, 0, Rsoo - tol, 2 * ls], ...
                                       'invert', true)));

	% strip the rim cells on either edge
    outersurfgcells = subset(outersurfgcells, ...
                             gcell_select(faeprob.fens, outersurfgcells, ...
                                struct('box', [Rsoo, -Rsoo, Rsoo, -Rsoo, 0, 0], ...
                                       'inflate', tol, ...
                                       'invert', true)));
                                   
    outersurfgcells = subset(outersurfgcells, ...
                             gcell_select(faeprob.fens, outersurfgcells, ...
                                struct('box', [Rsoo, -Rsoo, Rsoo, -Rsoo, ls, ls], ...
                                       'inflate', tol, ...
                                       'invert', true)));

end


