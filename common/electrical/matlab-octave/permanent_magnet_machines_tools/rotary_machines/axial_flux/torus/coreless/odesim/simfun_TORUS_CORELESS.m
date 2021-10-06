function [design, simoptions] = simfun_TORUS_CORELESS(design, simoptions)
% gathers data on a coreless torus type axial flux permanent magnet machine
% in preperation for a dynamic simulation.
%
% Syntax
%
% [design, simoptions] = simfun_TORUS_CORELESS(design, simoptions)
%
% 

    design.Hc = design.tc;
    design.Wc = (design.tauco - design.tauci) / 2;
    
    [design, simoptions] = simfun_TORUS(design, simoptions);
    
    design.gap = 2*design.g + design.tc;
    
    if ~simoptions.SkipFEA
        
        simoptions.MagFEASim.MagnetRegionMeshSize = -1;
        simoptions.MagFEASim.BackIronRegionMeshSize = -1;
        simoptions.MagFEASim.OuterRegionsMeshSize = [-1, -1];
        simoptions.MagFEASim.AirGapMeshSize = -1;
        
        if ~simoptions.SkipMainFEA 

            % Draw the main sim to extract the vector potential
            [design.FemmProblem, design.outermagsep] = corelessfemmprob_torus(design, ...
                                                    'NStages', simoptions.ndrawnstages, ...
                                                    'DrawCoils', false, ...
                                                    'Position', 0, ...
                                                    'MagnetRegionMeshSize', simoptions.MagFEASim.MagnetRegionMeshSize, ...
                                                    'BackIronRegionMeshSize', simoptions.MagFEASim.BackIronRegionMeshSize, ...
                                                    'OuterRegionsMeshSize', simoptions.MagFEASim.OuterRegionsMeshSize, ...
                                                    'AirGapMeshSize', simoptions.MagFEASim.AirGapMeshSize);

            % make an appropriate name for the .fem file from the base name
            % created by simfun_ROTARY, this will probably be a temporary file
            % location on the local system
            femfilename = [simoptions.filenamebase, '.fem'];

            if isfield(design, 'HcMag')
                design.FemmProblem.Materials(2).H_c = design.HcMag;
        %         design.FemmProblem.Materials(2).Name = regexprep(design.FemmProblem.Materials(2).Name, '\d+', sprintf('%3.0d', hc2mgoe(design.HcMag)));
            end

        %     femfilename = 'temp_simfun_TORUS_CORELESS.fem';

            % write the fem file to disk
            writefemmfile(femfilename, design.FemmProblem);
            
            ansfilename = analyse_mfemm(femfilename);
            
            % set up some variabls for the data extraction
            npoints = 45;
            xmin = -(design.outermagsep/2) + (0.00001*design.gap);
            xmax = xmin + design.gap - (2*0.00001*design.gap);
            xrange = linspace(xmin, xmax, npoints);
            yrange = linspace(0, 1, npoints);
            [design.X, design.Y] = meshgrid(xrange,yrange);
            xycoords = [design.X(:), design.Y(:)];
            
            if (exist('fpproc_interface_mex')==3)
                % using xfemm interface
                femmsolution = fpproc(ansfilename);
                
                % Extract various pieces of information from the simulation
                p = femmsolution.getpointvalues(design.X(:), design.Y(:) * design.taupm);

                design.A = reshape(p(1,:), size(design.X));
                design.Bx = reshape(p(2,:), size(design.X));
                design.By = reshape(p(3,:), size(design.X));

                % get the forces
                femmsolution.clearblock();
                femmsolution.groupselectblock(simoptions.ndrawnstages+1)
                design.gforce = -femmsolution.blockintegral(18)/simoptions.ndrawnstages;
                design.gvar = design.g;
                
                clear femmsolution;
            else
                % using original femm interface
                
                % open the solution in femm
                opendocument(ansfilename);
                
                % Extract various pieces of information from the simulation
                p = mo_getpointvalues(design.X(:), design.Y(:) * design.taupm);
                
                design.A = reshape(p(:,1), size(design.X));
                design.Bx = reshape(p(:,2), size(design.X));
                design.By = reshape(p(:,3), size(design.X));

                % get the forces
                mo_clearblock();
                mo_groupselectblock(simoptions.ndrawnstages+1)
                design.gforce = -mo_blockintegral(18)/simoptions.ndrawnstages;
                design.gvar = design.g;
                
                mo_close;
            
            end
            
            % clean up by deleting the .fem and .ans files
            delete(femfilename);
            delete(ansfilename);
            
            % Fit a polynomial to the vector potential in the gap between the torus
            % rotors
            design.APoly = machineApoly(xycoords, 8, design.A(:));
            
            % Fit a polynomial to the x and y directed flux density in the air gap
            [design.p_Bx, design.p_By] = machineBpolys(xycoords, 8, [design.Bx(:), design.By(:)]);
            
            % now get some more data points for forces
            simoptions = setfieldifabsent(simoptions, 'NForcePoints', 4);

            pos = linspace(design.FEMMTol - design.g, 0, simoptions.NForcePoints-1);
            pos(end) = 2*design.g;

            if simoptions.GetVariableGapForce
                design.gforce = [design.gforce, rotorforces_TORUS_CORELESS(design, simoptions.ndrawnstages, 2, pos)];
            else
                design.gforce = [design.gforce, repmat(design.gforce, 1, numel(pos))];
            end
            design.gvar = [design.gvar, design.g + pos];

        end
        
        if ~simoptions.SkipInductanceFEA
            
            % make an inductance sim
            Lcurrent = inductancesimcurrent(design.CoilArea, design.CoilTurns);
            
            InductanceFemmProblem = corelessLfemmprob_torus(design, ...
                                        'NStages', simoptions.ndrawnstages, ...
                                        'CoilCurrent', Lcurrent, ...
                                        'MagnetRegionMeshSize', simoptions.MagFEASim.MagnetRegionMeshSize, ...
                                        'BackIronRegionMeshSize', simoptions.MagFEASim.BackIronRegionMeshSize, ...
                                        'OuterRegionsMeshSize', simoptions.MagFEASim.OuterRegionsMeshSize, ...
                                        'AirGapMeshSize', simoptions.MagFEASim.AirGapMeshSize);

            % write the fem file to disk
            writefemmfile(femfilename, InductanceFemmProblem);

            ansfilename = analyse_mfemm(femfilename);
            
            % open the solution in femm
            opendocument([femfilename(1:end-4), '.ans']);

            % Now get the resistance and inductance of the machine coils
            [design.CoilResistance, design.CoilInductance] = RandLfromFEMMcircuit('1');

            % get the mutual inductance by dividing the flux in a neighbouring
            % phase by the applied current in the first phase
            temp = mo_getcircuitproperties('2');
            design.CoilInductance(2) = abs(temp(3)) / Lcurrent;

            mo_close;
            delete(femfilename);
            delete(ansfilename);
        
        end
    
    end
    
end


% 
% function LAint = mutualinductance(design)
% 
%         % now, get the necessary data to perform a simulation
%         npoints = 45;
%         xmin = -(design.outermagsep/2) + (0.00001*design.gap);
%         xmax = xmin + design.gap - (2*0.00001*design.gap);
%         xrange = linspace(xmin, xmax, npoints);
%         yrange = linspace(0, 1, npoints);
%         [design.X, design.Y] = meshgrid(xrange,yrange);
% 
%         % Extract various pieces of information from the simulation
%         p = mo_getpointvalues(design.X(:), design.Y(:) * design.taupm);
% 
%         design.Bx = reshape(p(:,2), size(design.X));
% 
%         design.By = reshape(p(:,3), size(design.X));
% 
%         xycoords = [design.X(:), design.Y(:)];
% 
% end