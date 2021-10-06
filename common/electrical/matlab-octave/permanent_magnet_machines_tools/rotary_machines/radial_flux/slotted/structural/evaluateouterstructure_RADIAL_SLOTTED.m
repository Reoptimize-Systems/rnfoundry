function [mingap, maxstress, design] = evaluateouterstructure_RADIAL_SLOTTED(design, simoptions, omega)
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

    % constrain only one end by default
    design.OuterStructure = setfieldifabsent(design.OuterStructure, 'OuterConstraint', [1, 0]);
    
    faeprob = mechfaesorprob;
    % Rr
    if strcmp(design.ArmatureType, 'internal')
        design.OuterStructure.Rsoi = design.Rbi;
        design.OuterStructure.Rsoo = design.Rbo;
    else
        design.OuterStructure.Rsoi = design.Ryi;
	    design.OuterStructure.Rsoo = design.Ryo;
    end

    % create the outer structure mesh
%     [faeprob,innersurfacelabel] = outerstructuremesh_RF( faeprob, ...
%                                                          design.OuterStructure.Rsii, ...
%                                                          design.OuterStructure.Rsoi, ...
%                                                          design.OuterStructure.Rsoo, ...
%                                                          design.ls, ...
%                                                          design.OuterStructure.lst, ...
%                                                          design.OuterStructure.OuterConstraint );
%                                                      
    [faeprob,rimlabel] = outerstructuremesh_RF_2( faeprob, ...
                                                  design.OuterStructure.Rsii, ...
                                                  design.OuterStructure.Rsio, ...
                                                  design.OuterStructure.Rsoi, ...
                                                  design.OuterStructure.Rsoo, ...
                                                  design.OuterStructure.thetasp, ...
                                                  design.ls, ...
                                                  design.OuterStructure.lst, ...
                                                  design.OuterStructure.lp, ...
                                                  design.OuterStructure.NSpokes, ...
                                                  design.OuterStructure.OuterConstraint );
                                          
    %% Finite Element Block

    nu = 0.31; 

    % tic
    % create the finite element  %fens, gcells, Rsio, E, nu
    [faeprob,mater] = outerstructurefeb_RF( faeprob, ...
                                            design.OuterStructure.Rsii, ...
                                            design.ls, ...
                                            design.OuterStructure.lst, ...
                                            design.OuterStructure.lp, ...
                                            simoptions.Evaluation.StructModulusOfElasticity, ...
                                            simoptions.Evaluation.StructPoissonRatio, ...
                                            design.OuterStructure.OuterConstraint );
    % toc

    %% Add the forces

    % calculate the shear forces
    shearforce = -(design.TorquePtoPeak / design.OuterStructure.Rsoi) / (2 * pi * design.OuterStructure.Rsoi * design.ls);

    % get the outer radial force
    radialforce = polyvaln(design.p_gforce, design.g) / (design.OuterStructure.Rsoi * design.thetap * design.ls);

    % apply the forces to the finite element block
    [Kmat, Fmat] = outerstructurestresses_RF( faeprob, ...
                                              design.OuterStructure.Rsoi, design.ls, ...
                                              shearforce, ...
                                              radialforce, ...
                                              simoptions.Evaluation.StructMaterialDensity, ...
                                              omega, ...
                                              rimlabel );

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
    
    design.OuterStructMaterialVolume = measure(faeprob.feb, faeprob.geom, inline('1')); 
    
    design.OuterStructMaterialMass = design.StructMaterialVolume * simoptions.Evaluation.StructMaterialDensity;
    
    % get the minimim air gap by finding the minimum radial distance from a position displaced by
    
    [mingap, maxstress, design] = evaluateouterstructure_RADIAL(design, simoptions, faeprob);

end



function [mingap, maxstress, design] = evaluateouterstructure_RADIAL(design, simoptions, faeprob)

    xyz = getvals(faeprob.geom);
    
    cartdispxyz = getvals(faeprob.u) + getvals(faeprob.geom);
    
    cartdispxyz = cartdispxyz(xyz(:,3) > 0 & xyz(:,3) < design.ls, :);
    
%     poldispxyz = zeros(size(cartdispxyz));
    
    [~, R, ~] = cart2pol(cartdispxyz(:,1), cartdispxyz(:,2), cartdispxyz(:,3));
    
    mingap = min(R - (design.OuterStructure.Rsoi - design.g)); 
    
    % obtain the maximum stress experienced anywhere
    
    % Obtain the continuous (Cauchy) stress field
    stressfield  = field_from_integration_points(faeprob.feb, faeprob.geom, faeprob.u, [], 'Cauchy', 6);
    
    maxstress = max(abs(get(stressfield,'values')));

end