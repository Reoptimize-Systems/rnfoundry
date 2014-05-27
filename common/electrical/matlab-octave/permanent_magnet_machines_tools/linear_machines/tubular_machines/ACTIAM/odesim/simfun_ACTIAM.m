function [design, simoptions] = simfun_ACTIAM(design, simoptions)
% generates electromagnetic data for the slotless tubular permanent magnet
% synchronous machine
%
% Syntax
%
% [design, simoptions] = simfun_ACTIAM(design, simoptions)
%
% 

    design = ratios2dimensions_ACTIAM (design);
    
    % set up variables for calculation of the iron losses, if they have not
    % been supplied
    if ~isfield(design, 'CoreLoss')
        % Coreloss will be the armature back iron data
        [design.CoreLoss.fq, ...
         design.CoreLoss.Bq, ...
         design.CoreLoss.Pq ] = m36assheared26gagecorelossdata(false);
    end
    
    [design, simoptions] = simfun_TM(design, simoptions);
    
    if ~any(design.mode == [2,3])
        error('Design mode must be 2 or 3 with steel in middle of sim to obtain flux linkage.');
    end
    
    % Draw and analyse the machine
    design.FemmProblem = femmprob_ACTIAM(design);
    femmfilename = [tempname, '.fem'];
    % generate the file
    writefemmfile(femmfilename, design.FemmProblem);
    % analyse
    ansfilename = analyse_mfemm(femmfilename, ...
                                simoptions.usefemm, ...
                                simoptions.quietfemm);
    % load solution
    solution = fpproc(ansfilename);

    % select locations at which the field quantities will be sampled
    npoints = 45;
    xrange = linspace(design.Rm+design.FEMMTol, design.Ro-design.FEMMTol, npoints);
    yrange = linspace(0, 1, npoints);
    [design.X, design.Y] = meshgrid(xrange,yrange);
    xycoords = [design.X(:), design.Y(:)];
    
    % ensure all points are within the bounds of the solution
    boundinds = xycoords(:,2) == 1;
    xycoords(boundinds,2) = xycoords(boundinds,2) - eps(1);
    
    % get the values from the solution
    p = solution.getpointvalues(xycoords(:,1), xycoords(:,2) * design.Wp);
    
    % reset the boundary points to actually be at the boundary
    xycoords(boundinds,2) = 1;
    
    % stroe the indformation in the design structure
    design.A = reshape(p(1,:)', size(design.X));
    design.Bx = reshape(p(2,:)', size(design.X));
    design.By = reshape(p(3,:)', size(design.X));
    
    design.APoly = machineApoly(xycoords, 8, p(1,:)');
    
    [design.p_Bx, design.p_By] = machineBpolys(xycoords, 8, [design.Bx(:), design.By(:)]);
    
    % extract the information necessary to calculate the losses in the
    % back iron material, note that the directions of the field are changed
    % when stored, as the FEA sim y-direction is not the direction of
    % motion of the translator.
    hbi = (design.Ra - design.Ro);
    design.CoreLoss(1).dy = hbi / 10;
    design.CoreLoss(1).coreycoords = ((design.Ro + design.CoreLoss(1).dy/2):design.CoreLoss(1).dy:(design.Ra-design.CoreLoss(1).dy/2))';
    
    design.CoreLoss(1).dx = min([0.01 / 5, design.Wp / 4]);
    design.CoreLoss(1).corexcoords = ((design.CoreLoss(1).dx/2):design.CoreLoss(1).dx:(design.Wp-design.CoreLoss(1).dx/2))';
    [meshx, meshy] = meshgrid(design.CoreLoss(1).coreycoords,design.CoreLoss(1).corexcoords);
    
    p = solution.getpointvalues(meshx(:), meshy(:));
    
    design.CoreLoss(1).By = reshape(p(2,:)', size(meshx))';
    design.CoreLoss(1).Bx = reshape(p(3,:)', size(meshx))';
    design.CoreLoss(1).Bz = zeros(size(design.CoreLoss(1).Bx));
    design.CoreLoss(1).Hy = reshape(p(6,:)', size(meshx))';
    design.CoreLoss(1).Hx = reshape(p(7,:)', size(meshx))';
    design.CoreLoss(1).Hz = zeros(size(design.CoreLoss(1).Hx));
    % make dz value equal to the circumference at each position
    design.CoreLoss(1).dz = pi * 2 * meshx';
    
    delete(ansfilename);
    delete(femmfilename);

    % Now get the inductances etc.
    
    % Draw and analyse the machine
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
                                simoptions.usefemm, ...
                                simoptions.quietfemm);
    % load solution
    solution = fpproc(ansfilename);
    
    % Now get the resistance and inductance of the machine coils
    [design.CoilResistance, design.CoilInductance] = solution.circuitRL('Coil B');
    
    delete(ansfilename);
    delete(femmfilename);

end