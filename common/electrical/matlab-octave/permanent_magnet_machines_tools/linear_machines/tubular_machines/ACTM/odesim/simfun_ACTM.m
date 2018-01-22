function [design, simoptions] = simfun_ACTM(design, simoptions)
% simfun_ACTM: performs FEA to gather data for the simulation of an ACTM
% machine design
%
% Syntax
%
% [design, simoptions] = simfun_ACTM(design, simoptions)
%

    design = ratios2dimensions_ACTM(design);
    
    design.Rcm = mean([design.Ro, design.Ri]);
    
    [design, simoptions] = simfun_TM(design, simoptions);

    design = setfieldifabsent (design, 'MagFEASimMaterials', struct ());
    design.MagFEASimMaterials = setfieldifabsent (design.MagFEASimMaterials, 'ArmatureIron', 'Air');
    design.MagFEASimMaterials = setfieldifabsent (design.MagFEASimMaterials, 'Magnet', 'NdFeB 40 MGOe');
    design.MagFEASimMaterials = setfieldifabsent (design.MagFEASimMaterials, 'ArmatureCoil', '10 AWG');
    design.MagFEASimMaterials = setfieldifabsent (design.MagFEASimMaterials, 'FieldIron', '1117 Steel');
    design.MagFEASimMaterials = setfieldifabsent (design.MagFEASimMaterials, 'Shaft', '304 Stainless Steel');
    design.MagFEASimMaterials = setfieldifabsent (design.MagFEASimMaterials, 'LibraryLocation', '');
    
    if ~simoptions.SkipMainFEA
        
        % Draw and analyse the machine
        design.FemmProblem = femmprob_ACTM(design);
        femmfilename = [tempname, '.fem'];
        % generate the file
        writefemmfile(femmfilename, design.FemmProblem);
        % analyse
        ansfilename = analyse_mfemm(femmfilename, ...
                                    simoptions.MagFEASim.UseFemm, ...
                                    simoptions.MagFEASim.QuietFemm);
        % load solution
        solution = fpproc(ansfilename);

        npoints = 45;
        xrange = linspace(design.Rm+design.FEMMTol, design.Ro + 2*design.g, npoints);
        yrange = linspace(0, 1, npoints);
        [design.X, design.Y] = meshgrid(xrange,yrange);

        p = solution.getpointvalues(design.X(:), design.Y(:) * design.Wp);

        xycoords = [design.X(:), design.Y(:)];

        % store the data in the design structure
        design.A = reshape(p(1,:)', size(design.X));
        design.Bx = reshape(p(2,:)', size(design.X));
        design.By = reshape(p(3,:)', size(design.X));

        design.APoly = machineApoly(xycoords, 8, design.A(:));
        [design.p_Bx, design.p_By] = machineBpolys(xycoords, 8, [design.Bx(:) design.By(:)]);

        % clean up
        delete(femmfilename);
        delete(ansfilename);
    
    end
    
    design = checkcoilprops_AM (design);
    
    if ~simoptions.SkipInductanceFEA
    
        % Now get the inductances etc.
        I = inductancesimcurrent(design.CoilArea, design.CoilTurns, 0.1e6);
        Ldesign = design;
        % set the magnet coercivity to zero for inductance sim
        Ldesign.HcMag = 0;
        Ldesign.FemmProblem = femmprob_ACTIAM(Ldesign, 'CoilCurrents', [0 I 0]);
        femmfilename = [tempname, '.fem'];
        % generate the file
        writefemmfile(femmfilename, Ldesign.FemmProblem);
        % analyse
        ansfilename = analyse_mfemm(femmfilename, ...
                                    simoptions.MagFEASim.UseFemm, ...
                                    simoptions.MagFEASim.QuietFemm);
        % load solution
        solution = fpproc(ansfilename);

        % Now get the resistance and inductance of the machine coils
        [design.CoilResistance, design.CoilInductance] = solution.circuitRL('Coil B');

        delete(ansfilename);
        delete(femmfilename);
        
    end
    
end
