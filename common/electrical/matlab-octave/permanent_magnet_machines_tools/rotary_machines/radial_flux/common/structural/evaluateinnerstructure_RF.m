function [maxzdef, maxstress, design] = evaluateinnerstructure_RF(design, simoptions, omega)
% finds the largest deflections in the structure of a radial flux machine
% rotor
%
% Syntax
%
% [maxzdef, design] = evaluateinnerstructure_RF(design, simoptions)
% 
% 

    if nargin < 2
        omega = 0;
    end
    
    if strcmp(design.StatorType, 'so')
        
        % create the inner structure mesh
        [fens,gcells,labels] = innerstructuremesh_RF( design.Rsii, ...
                                                      design.Rsio, ...
                                                      design.Ryi, ...
                                                      design.Ryo, ...
                                                      design.thetasp, ...
                                                      design.ls, ...
                                                      design.NStructSpokes, ...
                                                      layers, ...
                                                      ncircumpoints, ...
                                                      ncircumsppoints, ...
                                                      nradsipoints, ...
                                                      nradsppoints, ...
                                                      nradsopoints, ...
                                                      labels );

        %% Finite Element Block

        nu = 0.31;

        % tic
        % create the finite element  %fens, gcells, Rsio, E, nu
        [feb,geom,u,mater] = innerstructurefeb_RF(fens, gcells, design.Rbi, simoptions.evaloptions.E, nu);
        % toc

        %% Add the forces

        % calculate the shear forces
        shearforce = (design.TorquePtoPeak / design.Ryo) / (2 * pi * design.Ryo * design.ls);
        
        % get the outer radial force
        radialforce = polyvaln(design.p_gforce, design.g) / (2 * pi * design.Ryo * design.thetap * design.ls);
    
        % apply the forces to the finite element block
        [Kmat, Fmat] = inerstructurestresses_RF( fens, gcells, feb, geom, u, ...
                                                 design.Ryo, design.ls, ...
                                                 shearforce, ...
                                                 radialforce, ...
                                                 structdensity, ...
                                                 omega, ...
                                                 labels );
                                             
    else
        
        % calculate the shear forces
        shearforce = (design.TorquePtoPeak / design.Rbo) / (2 * pi * design.Ryo * design.ls);
        
        % get the outer radial force
        radialforce = polyvaln(design.p_gforce, design.g) / (2 * pi * design.Rbo * design.thetap * design.ls);
        
        % apply the forces to the finite element block
        [Kmat, Fmat] = inerstructurestresses_RF( fens, gcells, feb, geom, u, ...
                                                 design.Rbo, design.ls, ...
                                                 shearforce, ...
                                                 radialforce, ...
                                                 structdensity, ...
                                                 omega, ...
                                                 labels );
        
    end
    % toc
%     disp('done system setup');

    %%
    % Solve
    %
    % tic
    x = Kmat \ Fmat;
    % toc
    % tic
    % % use an iterative method to save memory
    % x = bicgstab(Kmat, Fmat, [], 5000);
    % toc
    u = scatter_sysvec(u, x);
    
    design.StructMaterialVolume = measure(feb, geom, inline('1')) * (design.NModules / simoptions.evaloptions.structmeshoptions.ModuleSubset);
    
    design.StructMaterialMass = design.StructMaterialVolume * simoptions.evaloptions.StructMaterialDensity;

%     disp('done system solving');
    
    % get the maximum deflection in the z direction
    maxdef = max(getvals(u));
    
    maxzdef = maxdef(3);
    
    % obtain the maximum stress experienced anywhere
    
    % Obtain the continuous (Cauchy) stress field
    stressfield  = field_from_integration_points(feb, geom, u, [], 'Cauchy', 6);
    maxstress = max(abs(get(stressfield,'values')));

end


