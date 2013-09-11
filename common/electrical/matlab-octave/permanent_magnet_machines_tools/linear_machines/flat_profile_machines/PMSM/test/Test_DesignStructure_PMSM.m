% Test_DesignStructure_PMSM
clear

load PMSM_Test_Design.mat

design.poles(2) = 3 * design.poles(2);

design.BeamSpreadFactor = 0;


design.sides = 2;

design.tols = [0;0;0];
design.AngleFromHorizontal = 0.5 * pi/2;

design.OuterWebs = 2;
design.sides = 2;
design.GuideRailIVars = [0.1, 0.1, 0.095, 0.095];
design.GuideRailIMethod = '1.3';
design.InnerStructureBeamVars = [];

options = designandevaloptions_PMSM;

options.nu = 0.31;

options.IMethod = '1.3';
options.gfactor = 0.5;

options.sections = 10;

design.alphab = 1.1;

beaminfofcn = @completebeaminfo_PMSM;
poleweightfcn = @fpoleweight_PMSM;
airgapclosurefcn = @feairgapclosure_PMSM;
structvolfcn = @structvol_PMSM;


design.OuterStructureBeamIMethod = options.IMethod;

% Load the database of available beams base don the beam type stored in
% the options structure
IVars = beamvars(options.IMethod);

IVarsLength = size(IVars, 2);

beams = zeros(size(IVars, 1), IVarsLength + 4);

beams(:, 1:IVarsLength) = IVars;

beams(:, IVarsLength + 1) = MomentOfInertiaY1(IVars, design.OuterStructureBeamIMethod);

% sort by moment of inertia
beams = sortrows(beams, IVarsLength + 1);

% create an index to beam with next biggest Moment of inertia
beams(1:end-1, IVarsLength + 2) = 2:size(beams,1);

% No stronger beam than last beam, so it should point to itself
beams(end, IVarsLength + 2) = size(beams,1);
IVarsLinkCol = size(IVars, 2) + 2;

for i = 1:size(IVars, 1)

    % get the number of pole support beams used in this design
    beams(i, IVarsLength + 4) = numbeams(span1(beams(i, 1:size(IVars, 2)), options.IMethod), ...
                                         design.PoleWidth, ...
                                         design.poles(2) * design.PoleWidth, ...
                                         design.BeamSpreadFactor);

    % Set the outer structure beam variables, the first row of these
    % will be the outer pole support variables, while the second is the
    % variables for the members to which these outer pole supports are
    % fixed (typically referrred to with the variable name
    % 'SupportBeams' elsewhere. The support beams are the next
    % strongest available beam to the beam used by the pole supports
    design.OuterStructureBeamVars = [beams(i, 1:IVarsLength); 
                                     beams( beams(i,IVarsLinkCol), 1:IVarsLength )];                            

    % Determine the total volume of material that would be required in
    % each case
    beams(i, IVarsLength + 3) = structvolfcn(design, options);

end

% sort by the total volume of the structural material
beams = sortrows(beams, IVarsLength + 3);
    
design.OuterStructureBeamVars = [beams(1, 1:IVarsLength); 
                                 beams( beams(1,IVarsLinkCol), 1:IVarsLength )];  

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

[design.fens, design.gcells] = createmeshfebeamdef2_FM(x, y, z, n, znodes, ...
    zguidecon, zsupp, zframe, zextreme, guideyshift, tolerance, BeamInfo);


design = designstructure_FM(design, options, structvolfcn, beaminfofcn, airgapclosurefcn);
