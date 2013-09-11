% Test_createmeshfebeamdef_FM

clear

load PMSM_Test_Design.mat

design.poles(2) = 2 * design.poles(2);

design.BeamSpreadFactor = 1;


design.sides = 2;

design.tols = [0;0;0];
design.AngleFromHorizontal = 0.5 * pi/2;

design.OuterWebs = 3;
design.sides = 2;
design.GuideRailIVars = [0.1, 0.1, 0.095, 0.095];
design.GuideRailIMethod = '1.3';
design.InnerStructureBeamVars = [];

options = designandevaloptions_PMSM;

options.nu = 0.31;

options.IMethod = '1.3';
options.gfactor = 0.5;

options.sections = 10;

% options.alphab = 1.05;

IVars = beamvars(options.IMethod);
    
% Now calculate the total volume of beam material that would be
% required if that beam were used. We do this by determining the number
% of beams required, as they are not distributed with regard for pole
% width.
n_beams = zeros(size(IVars, 1), 1);
totalVol = n_beams;

design.OuterStructureBeamIMethod = options.IMethod;

for i = 1:size(IVars, 1)

    n_beams(i,1) = numbeams(span1(IVars(i,:), options.IMethod), ...
        design.PoleWidth, ...
        design.poles(2) * design.PoleWidth, ...
        design.BeamSpreadFactor);

    design.OuterStructureBeamVars = IVars(i,:);

    % Determine the total volume of material that would be required in each
    % case
    totalVol(i,1) = structvol_PMSM(design, options);

end

beams = [IVars, totalVol, n_beams];

beams = sortrows(beams, size(IVars, 2)+1);
    
    
design.OuterStructureBeamVars = [beams(1, 1:size(IVars, 2)); beams(1, 1:size(IVars, 2))];

BeamInfo = completebeaminfo_PMSM(design, options);

[ x, y, z, n, ...
  znodes, ...
  zguidecon, ...
  zsupp, ...
  zframe, ...
  zextreme, ...
  guideyshift, ...
  tolerance ... 
 ] = dimsfebeamdef_FM(BeamInfo);

[fens, gcells, BeamInfo, OuterSupportCells] = createmeshfebeamdef2_FM(x, y, z, n, znodes, zguidecon, zsupp, zframe, zextreme, guideyshift, tolerance, BeamInfo)

%%




% Finite element block containing all the elements (useful for drawing)
feb = feblock_defor_ss_beam3 (struct ('gcells', gcells));

geom = field(struct ('name', 'geom', ...
                             'dim', 3, ...
                             'fens', fens));

ur  = field(struct ('name', 'ur', ...
                            'dim', 6, ...
                            'data', zeros(get(geom,'nfens'),6)));
% Graphics display
gv = graphic_viewer;

gv = reset (gv,[]);

% get the displacements
u = slice(ur, (1:3), 'u');

% Get the rotations
rot = slice(ur, (4:6), 'rot');

% Now draw the results
draw(feb, gv, struct ('x', geom, 'u', 0*u, 'rot', 0*rot, 'facecolor', 'none', 'drawscale', 0.2));

