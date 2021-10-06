function [feb,geom,u,mater] = axialrotorfeb(fens, gcells, Rbi, outersep, tbi, tsuppb, E, nu)
    
    % define a material with a young's modulus, E and poisson's ratio nu,
    % both specified above
    prop = property_linel_iso( struct('E', E, 'nu', nu) );

    mater = mater_defor_ss_linel_triax( struct('property', prop) );

    % Create a finite element block made of the gcells, all having material
    % mater and integration_rule gauss_rule
    feb = feblock_defor_ss( struct('mater', mater, ...
                                   'gcells', gcells,...
                                   'integration_rule', gauss_rule(3, 3)) );

    % Geometry: initialize a 3D geometry field (i.e. a field where the
    % values are the coordinates of the nodes)
    geom = field(struct ('name', 'geom', 'dim', 3, 'fens', fens));

    % Define the displacement field, initially all zeros. This field will
    % be the displacements of each of the points in the geom field
    u = geom * 0; % zero out

    % Apply EBC's (element boundary conditions)

    % first select all the nodes in a rectangular region. In this case
    % these are the nodes on
    maxz = (outersep/2 + tsuppb + tbi(1));
    ebc_fenids = fenode_select ( fens, struct('box', [-Rbi, Rbi, -Rbi, Rbi, maxz, maxz], 'inflate', 100*eps) );
    
    ebc_fenids = [ ebc_fenids, fenode_select( fens, ...
                                 struct('box', [-Rbi, Rbi, -Rbi, Rbi, -maxz, -maxz], 'inflate', 100*eps) ) ];

    % create a set of prescribed ebcs the same size as the selected nodes
    ebc_prescribed = ebc_fenids * 0 + 1;

    ebc_comp = [];

    ebc_val = ebc_fenids * 0;

    % apply the conditions
    u = set_ebc(u, ebc_fenids, ebc_prescribed, ebc_comp, ebc_val);

    u = apply_ebc(u);

    % Number equations
    u = numbereqns (u);

end