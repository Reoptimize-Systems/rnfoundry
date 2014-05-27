function [mingap, maxstress, design] = evaluateinnerstructure_RADIAL_SLOTTED(design, simoptions, omega)
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
    
    faeprob = mechfaesorprob;
    
    if strcmp(design.StatorType, 'so')
        design.InnerStructure.Rsoi = design.Ryi;
	    design.InnerStructure.Rsoo = design.Ryo;                                     
    else
        design.InnerStructure.Rsoi = design.Rbi;
	    design.InnerStructure.Rsoo = design.Rbo;
    end

    % create the inner structure mesh
    [faeprob,labels] = innerstructuremesh_RF( faeprob, ...
                                              design.InnerStructure.Rsii, ...
                                              design.InnerStructure.Rsio, ...
                                              design.InnerStructure.Rsoi, ...
                                              design.InnerStructure.Rsoo, ...
                                              design.InnerStructure.thetasp, ...
                                              design.ls, ...
                                              design.InnerStructure.NSpokes );

    %% Finite Element Block

    nu = 0.31;

    % tic
    % create the finite element  %fens, gcells, Rsio, E, nu
    [faeprob,mater] = innerstructurefeb_RF( faeprob, design.InnerStructure.Rsii, design.ls, ...
                            simoptions.evaloptions.StructModulusOfElasticity, ...
                            simoptions.evaloptions.StructPoissonRatio );
    % toc

    %% Add the forces

    % calculate the shear forces
    shearforce = (design.TorquePtoPeak / design.InnerStructure.Rsoo) / (2 * pi * design.InnerStructure.Rsoo * design.ls);

    % get the outer radial force
    radialforce = -polyvaln(design.p_gforce, design.g) / (design.InnerStructure.Rsoo * design.thetap * design.ls);

    % apply the forces to the finite element block
    [Kmat, Fmat] = innerstructurestresses_RF( faeprob, ...
                                              design.InnerStructure.Rsoo, design.ls, ...
                                              shearforce, ...
                                              radialforce, ...
                                              simoptions.evaloptions.StructMaterialDensity, ...
                                              omega, ...
                                              labels );

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
    faeprob.u = scatter_sysvec(faeprob.u, x);
    
    design.StructMaterialVolume = measure(faeprob.feb, faeprob.geom, inline('1')); 
    
    design.StructMaterialMass = design.StructMaterialVolume * simoptions.evaloptions.StructMaterialDensity;
    
    % get the minimim air gap by finding the minimum radial distance from a position displaced by
    
    [mingap, maxstress, design] = evaluateinnerstructure_RADIAL(design, simoptions, faeprob);

end


function [mingap, maxstress, design] = evaluateinnerstructure_RADIAL(design, simoptions, faeprob)

    cartdispxyz = getvals(faeprob.u) + getvals(faeprob.geom);
    
%     poldispxyz = zeros(size(cartdispxyz));
    
    [~, R, ~] = cart2pol(cartdispxyz(:,1), cartdispxyz(:,2), cartdispxyz(:,3));
    
    mingap = min(design.InnerStructure.Rsoo + design.g - R); 
    
    % obtain the maximum stress experienced anywhere
    
    % Obtain the continuous (Cauchy) stress field
    stressfield  = field_from_integration_points(faeprob.feb, faeprob.geom, faeprob.u, [], 'Cauchy', 6);
    
    maxstress = max(abs(get(stressfield,'values')));

end