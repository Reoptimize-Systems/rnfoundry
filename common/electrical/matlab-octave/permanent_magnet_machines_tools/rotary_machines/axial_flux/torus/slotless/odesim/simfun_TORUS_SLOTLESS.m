function [design, simoptions] = simfun_TORUS_SLOTLESS(design, simoptions)
% gathers simulation data for the slotless axial flux TORUS type machine in
% preparation for a dynamic simulation
%
% Syntax
%
% [design, simoptions] = simfun_TORUS_SLOTLESS(design, simoptions)
%
% 

    design.Hc = design.tc;
    design.Wc = design.tauco;
    
    % set up variables for calculation of the iron losses, if they have not
    % been supplied
    if ~isfield(design, 'CoreLoss')
        % Coreloss will be the armature back iron data
        [design.CoreLoss.fq, ...
         design.CoreLoss.Bq, ...
         design.CoreLoss.Pq ] = m36assheared26gagecorelossdata(false);
    end

    [design, simoptions] = simfun_TORUS(design, simoptions);
    
    design.gap = design.g + design.tc;
    
    if ~simoptions.SkipFEA

        % Draw the main sim to extract the vector potential
        [design.FemmProblem, design.outermagsep] = slotlessfemmprob_torus(design, ...
                                                'NStages', simoptions.ndrawnstages, ...
                                                'DrawCoils', true, ...
                                                'Position', 0, ...
                                                'MagnetRegionMeshSize', simoptions.MagFEASim.MagnetRegionMeshSize, ...
                                                'BackIronRegionMeshSize', simoptions.MagFEASim.BackIronRegionMeshSize, ...
                                                'OuterRegionsMeshSize', simoptions.MagFEASim.OuterRegionsMeshSize, ...
                                                'AirGapMeshSize', simoptions.MagFEASim.AirGapMeshSize);

        % get the material properties of the back Iron
        backironmatstruct = matstr2matstruct_mfemm(design.MagFEASimMaterials.FieldBackIron);
        
        % copy the permeability and resistivity to the appropriate field in
        % the design structure for calculation of losses later
        design.CoreMuR = backironmatstruct.Mu_x;
        design.CoreResistivity = 1 / backironmatstruct.Sigma / 1e6;
        
        % make an appropriate name for the .fem file from the base name created
        % by simfun_ROTARY
        femfilename = [simoptions.filenamebase, '.fem'];

        if isfield(design, 'HcMag')
            design.FemmProblem.Materials(2).H_c = design.HcMag;
    %         design.FemmProblem.Materials(2).Name = regexprep(design.FemmProblem.Materials(2).Name, '\d+', sprintf('%3.0d', hc2mgoe(design.HcMag)));
        end

    %     femfilename = 'temp_simfun_TORUS_CORELESS.fem';

        % write the fem file to disk
        writefemmfile(femfilename, design.FemmProblem);
        ansfilename = analyse_mfemm(femfilename);

        % now, get the necessary data to perform a simulation

        % we will fit two polynomials to the vector potential in the gap either
        % side of the yoke
        npoints = 45;
        xmin = -(design.outermagsep/2) + design.FEMMTol;
        xmax = xmin + design.gap - (2*design.FEMMTol);
        xrange = linspace(xmin, xmax, npoints);
        yrange = linspace(0, 1, npoints);
        [design.X, design.Y] = meshgrid(xrange,yrange);

        gridsize = size(design.X);

        xycoords = [design.X(:), design.Y(:)];

        if exist('fpproc_interface_mex', 'file') == 3
            % open the solution
            solution = fpproc(ansfilename);
            % Extract various pieces of information from the simulation
            p = solution.getpointvalues(xycoords(:,1), xycoords(:,2) .* design.taupm);
            
            design.A = reshape(p(1,:), gridsize);
            design.Bx = reshape(p(2,:), gridsize);
            design.By = reshape(p(3,:), gridsize);
        else
            opendocument(ansfilename);
            % Extract various pieces of information from the simulation
            p = mo_getpointvalues(xycoords(:,1), xycoords(:,2) .* design.taupm);
            
            design.A = reshape(p(:,1), gridsize);
            design.Bx = reshape(p(:,2), gridsize);
            design.By = reshape(p(:,3), gridsize);
        end
        
        design.MaxCoreSurfaceField = max(abs(design.Bx(:,end)));

        % Fit a polynomial to the vector potential in the gap between the left
        % hand torus rotor and the left hand side of the yoke
        design.APoly = machineApoly(xycoords, 8, design.A(:));

        [design.p_Bx, design.p_By] = machineBpolys(xycoords, 8, [design.Bx(:), design.By(:)]);

        xmin = -(design.outermagsep/2) + design.g + design.tc + design.ty + design.FEMMTol;
        xmax = xmin + design.gap - (2*design.FEMMTol);
        xrange = linspace(xmin, xmax, npoints);

        [design.X(:,:,2), design.Y(:,:,2)] = meshgrid(xrange,yrange);

        xycoords = [reshape(design.X(:,:,2), numel(design.X(:,:,2)), []),  ...
                    reshape(design.Y(:,:,2), numel(design.Y(:,:,2)), []) ];

        if exist('fpproc_interface_mex', 'file') == 3
            % Extract various pieces of information from the simulation
            p = solution.getpointvalues(xycoords(:,1), xycoords(:,2) .* design.taupm);
            
            design.A(:,:,2) = reshape(p(1,:), gridsize);
            design.Bx(:,:,2) = reshape(p(2,:), gridsize);
            design.By(:,:,2) = reshape(p(3,:), gridsize);
        else
            % Extract various pieces of information from the simulation
            p = mo_getpointvalues(xycoords(:,1), xycoords(:,2) .* design.taupm);
            
            design.A(:,:,2) = reshape(p(:,1), gridsize);
            design.Bx(:,:,2) = reshape(p(:,2), gridsize);
            design.By(:,:,2) = reshape(p(:,3), gridsize);
        end

        % Fit a polynomial to the vector potential in the gap between the torus
        % rotors
        design.APoly(2) = machineApoly(xycoords, 8, reshape(design.A(:,:,2),[],1));

        if exist('fpproc_interface_mex', 'file') == 3
            % get the forces
            solution.clearblock();
            solution.groupselectblock(simoptions.ndrawnstages+1)
            design.gforce = solution.blockintegral(18)/2;
        else
            % get the forces
            mo_clearblock();
            mo_groupselectblock(simoptions.ndrawnstages+1)
            design.gforce = mo_blockintegral(18)/2;
        end
        design.gvar = design.g;
        
        % extract the information necessary to calculate the losses in the
        % core material, note that the directions of the field are changed
        % when stored, as the FEA sim y-direction is not the direction of
        % motion of the translator. 
        design.CoreLoss.dy = min([0.01 / 5, design.ty / 4]);
        design.CoreLoss.corexcoords = -(design.outermagsep/2) + design.g + design.tc + ((design.CoreLoss.dy/2):design.CoreLoss.dy:(design.ty-design.CoreLoss.dy/2))';
        design.CoreLoss.dx = min([0.01 / 5, design.taupm / 4]);
        design.CoreLoss.coreycoords = ((design.CoreLoss.dx/2):design.CoreLoss.dx:(design.taupm-design.CoreLoss.dx/2))';
        [meshx, meshy] = meshgrid(design.CoreLoss.corexcoords,design.CoreLoss.coreycoords);
        
        if exist('fpproc_interface_mex', 'file') == 3
            p = solution.getpointvalues(meshx(:), meshy(:));

            design.CoreLoss.By = reshape(p(2,:), size(meshx))';
            design.CoreLoss.Bx = reshape(p(3,:), size(meshx))';
            design.CoreLoss.Hy = reshape(p(6,:), size(meshx))';
            design.CoreLoss.Hx = reshape(p(7,:), size(meshx))';
            clear solution;
        else
            p = mo_getpointvalues(meshx(:), meshy(:));

            design.CoreLoss.By = reshape(p(:,2), size(meshx))';
            design.CoreLoss.Bx = reshape(p(:,3), size(meshx))';
            design.CoreLoss.Hy = reshape(p(:,6), size(meshx))';
            design.CoreLoss.Hx = reshape(p(:,7), size(meshx))';
            mo_close;
            mi_close;
        end
        
        design.CoreLoss.Bz = zeros(size(design.CoreLoss.Bx));
        design.CoreLoss.Hz = zeros(size(design.CoreLoss.Hx));
        % make dz value equal to the depth of the simulation as it is used
        % to calculate the volume of material generating the loss (volume
        % will be dx X dy X dz)
        design.CoreLoss.dz = design.Rmo - design.Rmi;
        
        pos = linspace(design.FEMMTol - design.g, 0, simoptions.NForcePoints-1);
        pos(end) = 2*design.g;

        if simoptions.GetVariableGapForce
            design.gforce = [design.gforce, rotorforces_TORUS_SLOTLESS(design, simoptions.ndrawnstages, 2, pos)];
        else
            design.gforce = [design.gforce, repmat(design.gforce, 1, numel(pos))]; 
        end
        design.gvar = [design.gvar, design.g + pos];

        % clean up by removing the FEA files
        delete(femfilename);
        delete(ansfilename);
        
        % make an inductance sim
        Lcurrent = inductancesimcurrent(design.CoilArea, design.CoilTurns);
        InductanceFemmProblem = slotlessLfemmprob_torus(design, ...
                                    'NStages', simoptions.ndrawnstages, ...
                                    'CoilCurrent', Lcurrent, ...
                                    'MagnetRegionMeshSize', simoptions.MagFEASim.MagnetRegionMeshSize, ...
                                    'BackIronRegionMeshSize', simoptions.MagFEASim.BackIronRegionMeshSize, ...
                                    'OuterRegionsMeshSize', simoptions.MagFEASim.OuterRegionsMeshSize, ...
                                    'AirGapMeshSize', simoptions.MagFEASim.AirGapMeshSize);
                                
        % write the fem file to disk
        writefemmfile(femfilename, InductanceFemmProblem);
        ansfilename = analyse_mfemm(femfilename);

        if exist('fpproc_interface_mex', 'file') == 3
            solution = fpproc(ansfilename);
            [design.CoilResistance, design.CoilInductance] = solution.circuitRL('1');
            % get the mutual inductance by dividing the flux in a neighbouring
            % phase by the applied current in the first phase
            temp = solution.getcircuitprops('2');
        else
            opendocument(ansfilename);
            % Now get the resistance and inductance of the machine coils
            [design.CoilResistance, design.CoilInductance] = RandLfromFEMMcircuit('1');
            % get the mutual inductance by dividing the flux in a neighbouring
            % phase by the applied current in the first phase
            temp = mo_getcircuitproperties('2');
            mo_close;
            mi_close;
        end
        design.CoilInductance(2) = abs(temp(3)) / Lcurrent;
        % clean up the files
        delete(femfilename);
        delete(ansfilename);
    end
    
    % recalculate the resistance to account for end-windings etc
    design.MTL = rectcoilmtl(design.ty, design.Rmo + design.Rmi, design.tc);
    
    design.CoilResistance = 1.7e-8 * design.MTL ...
                            * design.CoilTurns ./ (pi * (design.Dc/2)^2);
             
end
