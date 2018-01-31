function [maxzdef, maxstress, design] = evaluatestructure_AF(design, simoptions)
% finds the largest deflections in the structure of an axial flux machine
% rotor
%
% Syntax
%
% [maxzdef, design] = evaluatestructure_AF(design, simoptions)
% 
% 

    simoptions.Evaluation.structmeshoptions.ModuleSubset = 1;
    
    [fens,gcells,maglabels,shaftlabels,outersep,simoptions] = rotormodulesmesh_AF(design, simoptions);

%     disp('done meshing');                

    %% Finite Element Block
    
    nu = 0.31;

    % tic
    % create the finite element block
    [feb,geom,u,mater] = axialrotorfeb(fens, gcells, design.Rbi, outersep, ...
                            design.tbi, design.tsuppb, simoptions.Evaluation.E, nu);
    % toc
    
    %%
    
    omega = max(simoptions.vT) / design.Rmm;
    
    magnetmass = simoptions.Evaluation.MagnetDensity * (pi * (realpow(design.Rmo, 2) - realpow(design.Rmi, 2)) * design.tm * (design.taumm/design.taupm));
    
    magtheta = design.taumm / design.Rmm;
    magarea = (0.5 * magtheta * (design.Rmo.^2 - design.Rmi.^2));
    
    % get the outer axial forces
    axialforces = polyvaln(design.p_gforce, design.g) / magarea;
    
    % invent some forces to ensure the inner back iron must have some
    % thickness
    axialforces = [axialforces, repmat(0.05*axialforces, 1, numel(maglabels)-2), -axialforces];
    
    [Kmat, Fmat] = axialrotorstresssys(design.Rmi, design.Rmo, design.tbi, ...
                                       design.MaxFpto, ...
                                       axialforces, ...
                                       magnetmass, ...
                                       design.Poles, ...
                                       simoptions.Evaluation.StructMaterialDensity, ...
                                       omega, ...
                                       maglabels, ...
                                       fens, gcells, feb, geom, u);
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
    
    design.StructMaterialVolume = measure(feb, geom, inline('1')) * (design.NModules / simoptions.Evaluation.structmeshoptions.ModuleSubset);
    
    design.StructMaterialMass = design.StructMaterialVolume * simoptions.Evaluation.StructMaterialDensity;

%     disp('done system solving');
    
    % get the maximum deflection in the z direction
    maxdef = max(getvals(u));
    
    maxzdef = maxdef(3);
    
    % obtain the maximum stress experienced anywhere
    
    % Obtain the continuous (Cauchy) stress field
    stressfield  = field_from_integration_points(feb, geom, u, [], 'Cauchy', 6);
    maxstress = max(abs(get(stressfield,'values')));

end


