function [Kmat, Fmat] = axialrotorstresssys(Rmi, Rmo, tbi, shearforce, axialforces, magnetmass, nmags, structdensity, omega, maglabels, fens, gcells, feb, geom, u)

    if isscalar(shearforce)
        shearforce = repmat(shearforce, size(maglabels));
    end
    
    if isscalar(axialforces)
        axialforces = [axialforces, zeros(1, numel(maglabels)-2), -axialforces];
    end
    
    % Assemble the system matrix
    ems = stiffness(feb, geom, u);
    K = start (sparse_sysmat, get(u, 'neqns'));
    K = assemble (K, ems);

    % Loads
    
    % initialise the force matrix
    F = start (sysvec, get(u, 'neqns'));

    [ botmagsurfgcells, topmagsurfgcells ] = getmagsurface(fens, gcells, maglabels(1), Rmo);
    
    magfeb = feblock_defor_ss(struct ('mater', getmater(feb), ...
                                      'gcells', topmagsurfgcells,...
                                      'integration_rule', gauss_rule(2, 2)));
                                  
    %   then the axial forces intensity object
    axialfi = force_intensity( struct ('magn', [0, 0, axialforces(1)]) );                           
    
    F = assemble (F, distrib_loads(magfeb, geom, u, axialfi, 2));
    
    [ botmagsurfgcells ] = getmagsurface(fens, gcells, maglabels(end), Rmo);
    
    magfeb = feblock_defor_ss(struct ('mater', getmater(feb), ...
                                      'gcells', botmagsurfgcells,...
                                      'integration_rule', gauss_rule(2, 2)));
                                  
    %   then the axial forces intensity object
    axialfi = force_intensity( struct ('magn', [0, 0, axialforces(end)]) );                           
    
    F = assemble (F, distrib_loads(magfeb, geom, u, axialfi, 2));
     
    % do the end rotor forces 
    magregionvol = (pi * (realpow(Rmo, 2) - realpow(Rmi, 2)) * tbi(1));
    
    % calculate the force per unit volume in the magnet region
    forcemag = shearforce(1) / magregionvol;
    
    % generate the shear force intensity object
    shearfi = force_intensity( struct ('magn', @(xyz) circtangentforcevec(xyz, forcemag)') );
%     % then the axial forces intensity object
%     axialfi = force_intensity( struct ('magn', [0, 0, axialforces(1)/magregionvol]) );
    
    % shear forces
    magcelllist = gcell_select(fens, gcells, struct ('label', maglabels(1)) );
    
    magfeb = feblock_defor_ss(struct ('mater', getmater(feb), ...
                                      'gcells', subset(gcells,magcelllist),...
                                      'integration_rule', gauss_rule(3, 2)));

    F = assemble (F, distrib_loads(magfeb, geom, u, shearfi, 3));
%     F = assemble (F, distrib_loads(magfeb, geom, u, axialfi, 3));
%     
    F = applymagcentrifugalload(F, fens, gcells, feb, geom, u, maglabels(1), [0, magnetmass], Rmo, Rmi, nmags, omega);
%                                   
    % other end
    % calculate the force per unit volume in the magnet region
    forcemag = shearforce(end) / magregionvol;
%     
    % first generate the shear force intensity object
    shearfi = force_intensity( struct ('magn', @(xyz) circtangentforcevec(xyz, forcemag)') );
%     % then the axial forces intensity object
%     axialfi = force_intensity( struct ('magn', [0, 0, axialforces(end)/magregionvol]) );
%     
    magcelllist = gcell_select(fens, gcells, struct ('label', maglabels(end)) );

    magfeb = feblock_defor_ss(struct ('mater', getmater(feb), ...
                                      'gcells', subset(gcells,magcelllist),...
                                      'integration_rule', gauss_rule(3, 2)));

    F = assemble (F, distrib_loads(magfeb, geom, u, shearfi, 3));
%     F = assemble (F, distrib_loads(magfeb, geom, u, axialfi, 3));
%     
    F = applymagcentrifugalload(F, fens, gcells, feb, geom, u, maglabels(end), [magnetmass, 0], Rmo, Rmi, nmags, omega);
    
    % now do the rest of the forces
    magregionvol = (pi * (realpow(Rmo, 2) - realpow(Rmi, 2)) * tbi(2));
    
    for i = 2:numel(maglabels)-1
    
        if axialforces(i)  > 0
            
            [ botmagsurfgcells, topmagsurfgcells ] = getmagsurface(fens, gcells, maglabels(1), Rmo);
            
            magfeb = feblock_defor_ss(struct ('mater', getmater(feb), ...
                                              'gcells', topmagsurfgcells,...
                                              'integration_rule', gauss_rule(2, 2)));
            
            %   then the axial forces intensity object
            axialfi = force_intensity( struct ('magn', [0, 0, axialforces(i)]) );
            
            F = assemble (F, distrib_loads(magfeb, geom, u, axialfi, 2));
            
        elseif axialforces < 0
            
            [ botmagsurfgcells ] = getmagsurface(fens, gcells, maglabels(end), Rmo);
            
            magfeb = feblock_defor_ss(struct ('mater', getmater(feb), ...
                                              'gcells', botmagsurfgcells,...
                                              'integration_rule', gauss_rule(2, 2)));
            
            %   then the axial forces intensity object
            axialfi = force_intensity( struct ('magn', [0, 0, axialforces(i)]) );
            
            F = assemble (F, distrib_loads(magfeb, geom, u, axialfi, 2));
            
%             
%             % calculate the force per unit volume in the magnet region
%             forcemag = shearforce(i) / magregionvol;
% 
%             % first generate the shear force intetensity object
%             shearfi = force_intensity( struct ('magn', @(xyz) circtangentforcevec(xyz, forcemag)') );
%             % then the axial forces intensity object
%             axialfi = force_intensity( struct ('magn', [0, 0, axialforces(i)/magregionvol]) );
% 
%             magcelllist = gcell_select(fens, gcells, struct ('label', maglabels(end)) );
% 
%             magfeb = feblock_defor_ss(struct ('mater', getmater(feb), ...
%                                               'gcells', subset(gcells,magcelllist),...
%                                               'integration_rule', gauss_rule(3, 2)));
% 
%             F = assemble (F, distrib_loads(magfeb, geom, u, shearfi, 3));
%             F = assemble (F, distrib_loads(magfeb, geom, u, axialfi, 3));
%             
        end
        
        F = applymagcentrifugalload(F, fens, gcells, feb, geom, u, maglabels(i), [magnetmass, magnetmass], Rmo, Rmi, nmags, omega);
        
    end
    
    % now apply radial cetrifugal forces to the entire rotor structure
    % the radial body force applied everywhere in the machine
    radialfi = force_intensity( struct ('magn', @(xyz) centrifugalbodyforce(xyz, omega, structdensity)') );
    
    F = assemble(F, distrib_loads(feb, geom, u, radialfi, 3));
    
    F = finish(F);

    Kmat = get(K,'mat');
    Fmat = get(F,'vec');

end


function magn = centrifugalbodyforce(xyz, omega, density)
% calculates the centrifugal force per unit volume throughout a mass
% rotating about the z axis
%
% Syntax
% 
% magn = centrifugalbodyforce(xyz, v, density)
%
% Input
%
%   xyz - coordinate at which the force is to be calculated
%
%   omega - velocity in radians/s
%
%   density - density of the mrotating material
%
% Output
%
%   magn - a vector of forces per unit mass acting at the point xyz in the
%          mass
%

    % get the radial distance from the center
    [theta,R] = cart2pol(xyz(1), xyz(2));
    
    % calculate the acceleration at this distance
    a = omega^2 * R;
    
    % the force per unit volume is density times the acceleration (the
    % total force would be F = ma = density * V * a, so force per unit
    % volume is just F = density * a)
    Fcentrifugal = density * a;
    
    magn = radialforcevec(xyz, Fcentrifugal);

end


function [F] = applymagcentrifugalload(F, fens, gcells, feb, geom, u, maglabel, magnetmass, Rmo, Rmi, Nm, omega)

    if any(magnetmass > 0) 
        [ botmagsurfgcells, topmagsurfgcells ] = getmagsurface(fens, gcells, maglabel, Rmo);
    end
    
    % calculate the force per m^2 on the plate face due to the magnet
    % centrifugal force, at the mean magnet radius
    magcentrifugalfpermsq = magnetmass * omega * Nm * (Rmo + Rmi) / (2 * pi * (Rmo^2 - Rmi^2));
    
    if magnetmass(1) > 0
        
        % the magnet suface centrifugal force intensity object
        magradfi = force_intensity( struct ('magn', @(xyz) radialforcevec(xyz, magcentrifugalfpermsq(1))') );

        % apply the centrifugal surface loads due to the magnets
        magfeb = feblock_defor_ss(struct ('mater', getmater(feb), ...
                                          'gcells', botmagsurfgcells,...
                                          'integration_rule', gauss_rule(2, 2)));

        F = assemble (F, distrib_loads(magfeb, geom, u, magradfi, 2));

    end
    
    if magnetmass(2) > 0
        
        % the magnet suface centrifugal force intensity object
        magradfi = force_intensity( struct ('magn', @(xyz) radialforcevec(xyz, magcentrifugalfpermsq(2))') );
        
        % apply the centrifugal surface loads due to the magnets
        magfeb = feblock_defor_ss(struct ('mater', getmater(feb), ...
                                          'gcells', topmagsurfgcells,...
                                          'integration_rule', gauss_rule(2, 2)));

        F = assemble (F, distrib_loads(magfeb, geom, u, magradfi, 2));
        
    end
    
end

function [ botmagsurfgcells, topmagsurfgcells ] = getmagsurface(fens, gcells, maglabel, Rmo)

    % first get the 3D gcells making up the magnet region in the fin
    maggcells = subset( gcells, gcell_select(fens, gcells, struct('label', maglabel)) );
    
    % then get the boundary cells of this region
    bdrygcells = mesh_bdry(maggcells);
    
    nodes = unique(getconn(bdrygcells));
    
    xyz = getxyz(subset(fens, nodes));
    
    minzpos = min(xyz(:,3));
    maxzpos = max(xyz(:,3));
    
    tol = (maxzpos - minzpos) / 1000;
    
    % then extract only parts of the surface of interest
    botmagsurfgcells = subset(bdrygcells, gcell_select(fens, bdrygcells, struct('box', [Rmo, -Rmo, Rmo, -Rmo, minzpos, minzpos], 'inflate', tol)));

    topmagsurfgcells = subset(bdrygcells, gcell_select(fens, bdrygcells, struct('box', [Rmo, -Rmo, Rmo, -Rmo, maxzpos, maxzpos], 'inflate', tol)));

end


