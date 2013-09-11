function [faeprob,mater] = outerstructurefeb_RF(faeprob, Rsii, ls, lst, lp, E, nu, fixed)
    
    % define a material with a young's modulus, E and poisson's ratio nu,
    % both specified above
    prop = property_linel_iso( struct('E', E, 'nu', nu) );

    mater = mater_defor_ss_linel_triax( struct('property', prop) );

    % Create a finite element block made of the gcells, all having material
    % mater and integration_rule gauss_rule
    faeprob.feb = feblock_defor_ss( struct('mater', mater, ...
                                   'gcells', faeprob.gcells,...
                                   'integration_rule', gauss_rule(3, 3)) );

    % Geometry: initialize a 3D geometry field (i.e. a field where the
    % values are the coordinates of the nodes)
    faeprob.geom = field(struct ('name', 'geom', 'dim', 3, 'fens', faeprob.fens));

    % Define the displacement field, initially all zeros. This field will
    % be the displacements of each of the points in the geom field
    faeprob.u = faeprob.geom * 0; % zero out

    % Apply EBC's (element boundary conditions)

%     % first select some nodes on one end of the outer structure cylinder
%     hpnts = hullpnts(Rsii, ls);
%     ebc_fenids = fenode_select ( faeprob.fens, struct('convhull', hpnts) );

    if fixed(1)
        ebc_fenids1 = fenode_select ( faeprob.fens, struct('cylinder', [0,0,-(lp + lst)/2,Rsii,(lp + lst)], ...
                                                           'inflate', 100*eps(Rsii) ) );
    else
        ebc_fenids1 = [];
    end
    
    if fixed(2)
        ebc_fenids2 = fenode_select( faeprob.fens, struct('cylinder', [0,0,ls+(lp+lst)/2,Rsii, (lst+lp) ], ...
                                                          'inflate', 100*eps(Rsii) ) );
    else
        ebc_fenids2 = [];
    end
      
	ebc_fenids = [ebc_fenids1, ebc_fenids2];
    
%     % first select some nodes on either end of the inner shaft
%     ebc_fenids = fenode_select ( faeprob.fens, struct('box', [0, 0, -Rsii, Rsii, 0, 0], 'inflate', 100*eps(Rsii)) );
%     
%     ebc_fenids = [ ebc_fenids, fenode_select( faeprob.fens, ...
%                                  struct('box', [0, 0, -Rsii, Rsii, ls, ls], 'inflate', 100*eps(Rsii)) ) ];

    % create a set of prescribed ebcs the same size as the selected nodes
    ebc_prescribed = ebc_fenids * 0 + 1;

    ebc_comp = [];

    ebc_val = ebc_fenids * 0;

    % apply the conditions
    faeprob.u = set_ebc(faeprob.u, ebc_fenids, ebc_prescribed, ebc_comp, ebc_val);

    faeprob.u = apply_ebc(faeprob.u);

    % Number equations
    faeprob.u = numbereqns(faeprob.u);

end

