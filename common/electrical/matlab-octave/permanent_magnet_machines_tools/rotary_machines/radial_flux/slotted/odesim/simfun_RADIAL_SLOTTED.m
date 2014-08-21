function [design, simoptions] = simfun_RADIAL_SLOTTED(design, simoptions)
% generates simulation data for a radial flux slotted machine in
% preparation for a time series ode simulation
%
% Syntax
%
% [design, simoptions] = simfun_RADIAL_SLOTTED(design, simoptions)
%
% Description
%
% simfun_RADIAL_SLOTTED takes a slotted radial flux machine design, specified in 
% the fields of a structure and performs a series of electromagnetic simulation
% to obtain  data on the machine performance. In addition, various other 
% calculations are  performed such as estimating the resistance and inductance 
% etc. The results of the calcualtions are added as fields to the design 
% structure, which is returned with the appended fields. Also expected to be
% provided is a simoptions structure containing fields which specify various
% options and settings determining what simulation data will be gathered. More
% information on the simoptions fields is given later in this help.
%
% The design can be for either an internal or external armature. This is
% specified in the field 'ArmatureType' which must contain a string, either
% 'internal' or 'external'. Any unambiguous substring is also acceptable,
% e.g 'e' or'i', or 'ext', or 'in'.
%
% Both internal and external armature designs must have the following fields:
%
%    tsg - thickness of shoe in radial direction at shoe gap
%
%    thetam - angular pitch of magnet in radians
%
%    thetac - angular pitch of coil slot, this can be a single value or two 
%      values. If two values, it is the pitch at the base (close to armature 
%      yoke) and top (close to slot opening) of the slot respectively. If this 
%      is a single value, the same pitch is used for both the base and the top.
%
%    thetasg - angular pitch of the coil slot opening between shoes.
%
%    ls - stack length (depth 'into the page' of simulation)
% 
% In general, all dimensions which refer to a radial measurement from the
% center of the machine are prefixed with the capital letter 'R'.
%
% For an INTERNAL ARMATURE machine, the design structure must also contain
% all the fields:
%
%    Rbo - radial distance to outer back iron surface
%
%    Rmo - radial distance to outer magnet surface
%
%    Rmi - radial distance to inner magnet surface
%
%    Rao - radial distance to armature outer surface (surface of teeth or coils)
%
%    Rtsb - radial distance to tooth shoe base
%
%    Ryi - radial distance to armature yoke inner surface
%
%    Ryo - radial distance to armature yoke outer surface
%
%
% For an EXTERNAL ARMATURE machine, the design structure must also contain the 
% fields:
%
%    Ryo - radial distance to armature yoke outer surface
%
%    Ryi - radial distance to armature yoke inner surface
%
%    Rtsb - radial distance to tooth shoe base
%
%    Rai - radial distance to armature inner surface (surface of teeth or coils)
%
%    Rmi - radial distance to inner magnet surface
%
%    Rmo - radial distance to outer magnet surface
%
%    Rbi - radial distance to back iron inner surface
%
%
% This completes the specification of the physical dimentions of the
% armature and field. Alternative methods of specifying the dimensions (such as 
% as a set of dimensionless ratios) are possible which can be easily converted 
% to this form using the helper function completedesign_RADIAL_SLOTTED.m. See 
% the help for  completedesign_RADIAL_SLOTTED for more details.
%
% In addition, a winding specification must be supplied. The winding is
% described using the following variables:
%
%  yp - Average coil pitch as defined by (Qs/Poles)
%  yd - Actual coil pitch as defined by round(yp) +/- k
%  Qs  -  total number of stator slots in machine
%  Qc  -  total number of winding coils in machine 
%  q  -  number of slots per pole and phase
%  qn  -  numerator of q
%  qd  -  denominator of q
%  qc - number of coils per pole and phase
%  qcn  -  numerator of qc
%  qcd  -  denominator of qc
%  Qcb - basic winding (the minimum number of coils required to make up a 
%    repetitive segment of the machine that can be modelled using symmetry)
%  pb - the number of poles corresponding to the basic winding in Qcb
%  CoilLayers - number of layers in the winding
%
% This pole/slot/coil/winding terminology is based on that presented in
% [1].
%
% Machine windings can be single or double layered, in which case:
%
% Single layer
%   q = 2qc
%   Qs = 2Qc
% Double layer
%   q = qc
%   Qs = Qc
%
% To specify a winding, the 'minimum spec' that must be provided is based on
% a combination of some or all of the following fields:
%
%  Phases - The number of phases in the machine
%  Poles - The number of magnetic poles in the machine
%  NBasicWindings - the number of basic winding segments in the machine
%  qc - number of coils per pole and phase (as a fraction object)
%  Qc - total number of coils (in all phases) in the machine
%  CoilLayers - the number of layers in the coil, either 1 or 2
%
% Any of the following combinations may be supplied to specify the winding:
%
%   Poles, Phases, Qc, CoilLayers
%   Poles, Phases, qc, CoilLayers
%   qc, Phases, NBasicWindings, CoilLayers
%
% These variables must be provided as fields in the design structure. If
% 'qc' is supplied, it must be an object of the class 'fr'. This is a class
% designed for handling fractions. See the help for the ''fr'' class for
% further information.
%
% In addition some further design option fields may be specified if desired:
%
%   NStrands               1
%
%   NStages                1
%
%   MagnetSkew             0
%
%   FEMMTol - tolerance used in FEA drawings to determine if components are
%    separate. Defaults to 1e-5 if not supplied.
%
%   WireResistivityBase    1.7e-8
%
%   AlphaResistivity       3.93e-3
%
%   TemperatureBase        20
%
%
% SIMOPTIONS 
%
% Also expected to be provided is a structure containing various simulation
% options.
%
%
%  usefemm - optional flag to allow forcing the use of FEMM for FEA simulation 
%   rather than mfemm. If true FEMM will be used, otherwise, mfemm will be used
%   if it is available. Defaults to false if not supplied.
%
%  quietfemm - optional flag determining if output to the command line from
%   mfemm should be suppressed. Defaults to true if not supplied.
%
%  GetVariableGapForce - optional flag determining if a series of force 
%   simulations with various air gap sizes will be performed. Must be supplied.
%
%  NForcePoints - optional number of gap sizes to use if GetVariableGapForce is
%   true. Defaults to 4 if not supplied. The force positions will be linearly
%   spaced from zero air gap (plus a small toperance value) to twice the air gap
%   length.
%
%  femmmeshoptions - 
%
%  SkipMainFEA - optional flag which allows skipping the main fea simulations.
%   This can be desirable when only coil winding configurations change etc.
%   Defaults to false if not supplied.
%
%  SkipInductanceFEA - optional flag which allows skipping the inductance fea
%   simulations. This can be desirable when only coil winding configurations
%   change etc. Defaults to false if not supplied.
%
%  
%
% [1] J. J. Germishuizen and M. J. Kamper, "Classification of symmetrical
% non-overlapping three-phase windings," in The XIX International
% Conference on Electrical Machines - ICEM 2010, 2010, pp. 1-6.
%
%
% See also: completedesign_RADIAL_SLOTTED.m, simfun_RADIAL.m, simfun_ROTARY.m
%           simfun_AM.m, fr.m
%
%
    
    if design.ypd ~= 1 && design.ypd ~= 2
    	error('denominator of slots per pole must be 1 or 2, other values not yet supported')
    end

    if ~isfield(design, 'CoreLoss')
        % CoreLoss will be the armature back iron data
        [ design.CoreLoss.kh, ...
          design.CoreLoss.kc, ...
          design.CoreLoss.ke, ...
          design.CoreLoss.beta ] = corelosscoeffs ('M-36', '26');
    end
    
    % We don't check the coil turns etc at this stage (done by default in 
    % simfun_RADIAL) as we don't yet know the coil cross-sectional area, this is
    % done later below
    simoptions.SkipCheckCoilProps = true;
    
    % call the common radial simulation function
    [design, simoptions] = simfun_RADIAL(design, simoptions);
    
    if design.tsg > 1e-5
        if design.tsb > 1e-5
            simoptions.femmmeshoptions = setfieldifabsent(simoptions.femmmeshoptions, 'ShoeGapRegionMeshSize', choosemesharea_mfemm(max(design.tsg, design.tsb), (design.Rmo*design.thetasg), 1/20));
        else
            simoptions.femmmeshoptions = setfieldifabsent(simoptions.femmmeshoptions, 'ShoeGapRegionMeshSize', choosemesharea_mfemm(design.tsb, (design.Rmo*design.thetasg), 1/20));
        end
    else
        if design.tsb > 1e-5
            simoptions.femmmeshoptions = setfieldifabsent(simoptions.femmmeshoptions, 'ShoeGapRegionMeshSize', choosemesharea_mfemm(design.tsb, (design.Rmo*design.thetasg), 1/20));
        else
            simoptions.femmmeshoptions = setfieldifabsent(simoptions.femmmeshoptions, 'ShoeGapRegionMeshSize', -1);
        end
    end
    
    % now get the flux linkage in the coils, we will do this by getting the
    % integral of the vector potential in a slot over a half pole of
    % displacement, giving a quarter of the wave period
    
    if ~simoptions.SkipMainFEA
    
        nfeapos = 10;
    %     design.FirstSlotCenter = design.thetap + (design.thetap / slotsperpole / 2);
        design.FirstSlotCenter = 0;
        design.feapos = linspace(0, 1, nfeapos);

        design.intAdata.slotPos = [];
        design.intAdata.slotIntA = [];
        design.intBdata.slotPos = [];
        design.intBdata.slotIntB = [];
        
        femfilename = [tempname, '_simfun_RADIAL_SLOTTED.fem'];
        
        % determine a unit vector pointing in the direction normal to the air
        % gap half way between the Poles for the purpose of extracting the
        % air-gap closing forces
        [gvector(1), gvector(2)] = pol2cart(design.thetap,1);
        gvector = unit(gvector);
        
        % determine a unit vector pointing in the direction tangential to the
        % radius half way between the Poles for the purpose of extracting
        % the cogging forces
        %
        % The dot product of orthogonal vectors is zero. In two dimensions the
        % slopes of perpendicular lines are negative reciprocals. Switch the x
        % and y coefficients and change the sign of one of them.
        % 
        % For vector V1 = <a, b>, the reciprocal V2 = <b, -a>. Now make it a
        % unit vector by dividing by the magnitude.
        coggingvector = unit([gvector(2), -gvector(1)]);
        
        armirongroup = 2;
            
        for i = 1:numel(design.feapos)

            % Draw the sim, i.e by creating the FemmProblem structure
            [design.FemmProblem, design.coillabellocs] = ...
                                slottedfemmprob_radial(design, ...
                                    'ArmatureType', design.ArmatureType, ...
                                    'NWindingLayers', design.CoilLayers, ...
                                    'Position', (design.feapos(i)+1) * design.thetap + design.FirstSlotCenter, ...
                                    'ArmatureBackIronGroup', armirongroup, ...
                                    'MagnetRegionMeshSize', simoptions.femmmeshoptions.MagnetRegionMeshSize, ...
                                    'BackIronRegionMeshSize', simoptions.femmmeshoptions.BackIronRegionMeshSize, ...
                                    'OuterRegionsMeshSize', simoptions.femmmeshoptions.OuterRegionsMeshSize, ...
                                    'AirGapMeshSize', simoptions.femmmeshoptions.AirGapMeshSize, ...
                                    'ShoeGapRegionMeshSize', simoptions.femmmeshoptions.ShoeGapRegionMeshSize, ...
                                    'YokeRegionMeshSize', simoptions.femmmeshoptions.YokeRegionMeshSize, ...
                                    'CoilRegionMeshSize', simoptions.femmmeshoptions.CoilRegionMeshSize);

            % write the fem file to disk
            writefemmfile(femfilename, design.FemmProblem);
            % analyse the problem
            ansfilename = analyse_mfemm(femfilename, ...
                                        simoptions.usefemm, ...
                                        simoptions.quietfemm);
                                
            
            if (exist('fpproc_interface_mex', 'file')==3) && ~simoptions.usefemm
                
                solution = fpproc(ansfilename);
                % activate field smoothing so values are interpolated across
                % mesh elements
                solution.smoothon();
                
                if i == 1
                    design = corelosssetup(design, design.feapos, solution);
                end
                
                % get the integral of the vector potential in a slot, if two
                % layers get both layers
                design = slotintAdata(design, simoptions, design.feapos(i), solution);
                
                design = slotintBdata(design, simoptions, design.feapos(i), solution);

                for ii = 1:numel(design.CoreLoss)
                    
                    p = solution.getpointvalues( design.CoreLoss(ii).meshx(:), ...
                                                 design.CoreLoss(ii).meshy(:) );

                    design.CoreLoss(ii).By(:,i,:) = reshape(p(2,:)', size(design.CoreLoss(ii).meshx));
                    design.CoreLoss(ii).Bx(:,i,:) = reshape(p(3,:)', size(design.CoreLoss(ii).meshx));

                end
                
                if i == 1
                    % get the forces
                    solution.clearblock();
                    solution.groupselectblock(1)
                    
                    design.gforce = dot([solution.blockintegral(18)/2, solution.blockintegral(19)/2], ...
                                         gvector);
                                    
                    design.gvar = design.g;
                    
                    % get the cross-sectional area of the armature iron for
                    % calcuation of material masses later
                    solution.clearblock();
                    solution.groupselectblock(armirongroup);
                    design.ArmatureIronAreaPerPole = solution.blockintegral(5)/2;
                    
                    % get the cross-sectional area of the coil winding bundle
                    design.CoilArea = solution.blockintegral ( 5, design.coillabellocs(1,1), design.coillabellocs(1,2) );
                    
                end
                
                % get the cogging forces
                solution.clearblock();
                solution.groupselectblock(1)
                
                design.coggingforce(i) = dot([solution.blockintegral(18)/2, solution.blockintegral(19)/2], ...
                                             coggingvector);

                design.coggingforce(i) = design.coggingforce(i) * design.Poles(1);
                
%                if simoptions.getintermagflux
%                    design.InterMagFlux = getintermagflux (design, solution);
%                end
                
                % explicitly call the delete method on the solution
                delete(solution);
                clear solution;
                
            else
                if i == 1
                    design = corelosssetup(design, design.feapos, ansfilename);
                end
                
                % open the solution in FEMM
                opendocument(ansfilename);

                % get the integral of the vector potential in a slot, if two
                % layers get both layers
                design = slotintAdata(design, simoptions, design.feapos(i));

                design = slotintBdata(design, simoptions, design.feapos(i));

                for ii = 1:numel(design.CoreLoss)
                    
                    p = mo_getpointvalues( design.CoreLoss(ii).meshx(:), ...
                                           design.CoreLoss(ii).meshy(:) );

                    design.CoreLoss(ii).By(:,i,:) = reshape(p(:,2), size(design.CoreLoss(ii).meshx));
                    design.CoreLoss(ii).Bx(:,i,:) = reshape(p(:,3), size(design.CoreLoss(ii).meshx));

                end
                
                if i == 1
                    % get the forces
                    mo_clearblock();
                    mo_groupselectblock(1);
                    design.gforce = dot([mo_blockintegral(18)/2, mo_blockintegral(19)/2], ...
                                        gvector);
                    design.gvar = design.g;
                    
                    % get the cross-sectional area of the armature iron for
                    % calcuation of material masses later
                    mo_clearblock();
                    mo_groupselectblock(armirongroup);
                    design.ArmatureIronAreaPerPole = mo_blockintegral(5)/2;
                    
                end

                mo_close;
                
            end
            
            % clean up the FEA files
            delete(femfilename);
            delete(ansfilename);
        
        end

    end % ~SkipMainFEA
    
    % now get more force points
    pos = linspace(design.FEMMTol - design.g, 0, simoptions.NForcePoints-1);
    pos(end) = 2*design.g;
    
    if simoptions.GetVariableGapForce
        design.gforce = [design.gforce, rotorforces_RADIAL_SLOTTED(design.FemmProblem, 2, pos)];
    else
        design.gforce = [design.gforce, repmat(design.gforce, 1, numel(pos))];
    end
    design.gvar = [design.gvar, design.g + pos];
    
    % make sure the winding properties (number of turns etc.) are up to date
    design = checkcoilprops_AM(design);
    
    if ~simoptions.SkipInductanceFEA
    
        % perform an inductance sim
        Lcurrent = inductancesimcurrent(design.CoilArea, design.CoilTurns);
        [InductanceFemmProblem, Lcoillabellocs] = slottedLfemmprob_radial(design, ...
            'ArmatureType', design.ArmatureType, ...
            'NWindingLayers', design.CoilLayers, ...
            'CoilCurrent', Lcurrent, ...
            'MagnetRegionMeshSize', simoptions.femmmeshoptions.MagnetRegionMeshSize, ...
            'BackIronRegionMeshSize', simoptions.femmmeshoptions.BackIronRegionMeshSize, ...
            'OuterRegionsMeshSize', simoptions.femmmeshoptions.OuterRegionsMeshSize, ...
            'AirGapMeshSize', simoptions.femmmeshoptions.AirGapMeshSize, ...
            'ShoeGapRegionMeshSize', simoptions.femmmeshoptions.ShoeGapRegionMeshSize, ...
            'YokeRegionMeshSize', simoptions.femmmeshoptions.YokeRegionMeshSize, ...
            'CoilRegionMeshSize', simoptions.femmmeshoptions.CoilRegionMeshSize);
        
        % analyse the problem
        [ansfilename, femfilename] = analyse_mfemm(InductanceFemmProblem, ...
                                                   simoptions.usefemm, ...
                                                   simoptions.quietfemm);
        
        if exist('fpproc_interface_mex', 'file') == 3 && ~simoptions.usefemm
            solution = fpproc(ansfilename);
            [design.CoilResistance, design.CoilInductance] = solution.circuitRL('1');
            % get the mutual inductance by dividing the flux in a neighbouring
            % phase by the applied current in the first phase
            temp = solution.getcircuitprops('2');
            % explicitly call the delete method on the solution
            delete(solution);
            clear solution;
        else
            opendocument(ansfilename);
            % Now get the resistance and inductance of the machine coils
            [design.CoilResistance, design.CoilInductance] = RandLfromFEMMcircuit('1');
            % get the mutual inductance by dividing the flux in a neighbouring
            % phase by the applied current in the first phase
            temp = mo_getcircuitproperties('2');
            mo_close;
        end
        design.CoilInductance(2) = abs(temp(3)) / Lcurrent;
        
        % clean up the FEA files
        delete(femfilename);
        delete(ansfilename);
    
    end % ~simoptions.SkipInductanceFEA
    
    % now recalculate coil resistance
    design.MTL = rectcoilmtl( design.ls, ...
                              design.yd * design.thetas * design.Rcm, ...
                              mean(design.thetac) * design.Rcm );

    design.CoilResistance = 1.7e-8 * design.CoilTurns * design.MTL ./ (pi * (design.Dc/2)^2);
    
end


function design = slotintBdata(design, simoptions, feapos, solution)

    % extract the flux integral data from all the slots at the given
    % positions of the magnet relative to coil
    if design.CoilLayers == 1
        
        slotypos = design.coillabellocs(1:(size(design.coillabellocs,1)),2);
        slotxpos = design.coillabellocs(1:(size(design.coillabellocs,1)),1);
                
        for i = 1:numel(slotypos)
            if exist('fpproc_interface_mex', 'file') == 3 && ~simoptions.usefemm
                % x and y directed flux in slot on left hand side
                solution.clearblock();
                solution.selectblock(slotxpos(i,1), slotypos(i));
                design.intBdata.slotIntB(end+1,1,1) = solution.blockintegral(8);
                design.intBdata.slotIntB(end,2,1) = solution.blockintegral(9);
            else
                % x and y directed flux in slot on left hand side
                mo_clearblock();
                mo_selectblock(slotxpos(i,1), slotypos(i));
                design.intBdata.slotIntB(end+1,1,1) = mo_blockintegral(8);
                design.intBdata.slotIntB(end,2,1) = mo_blockintegral(9);
            end
        end
        
    elseif design.CoilLayers == 2
        
        slotypos = [design.coillabellocs(1:2:(size(design.coillabellocs,1)),2), ...
                    design.coillabellocs((1:2:(size(design.coillabellocs,1)))+1,2)];
        slotxpos = [design.coillabellocs(1:2:(size(design.coillabellocs,1)),1), ...
                    design.coillabellocs((1:2:(size(design.coillabellocs,1)))+1,1)];
              
        for i = 1:size(slotypos,1)

            if exist('fpproc_interface_mex', 'file') == 3 && ~simoptions.usefemm
                % x and y directed flux in slot on left hand side outer
                % (leftmost) layer
                solution.clearblock();
                solution.selectblock(slotxpos(i,1), slotypos(i,1));
                design.intBdata.slotIntB(end+1,1,1) = solution.blockintegral(8);
                design.intBdata.slotIntB(end,2,1) = solution.blockintegral(9);
                
                % x and y directed flux in slot on left hand side inner
                % (rightmost) layer
                solution.clearblock();
                solution.selectblock(slotxpos(i,2), slotypos(i,2));
                design.intBdata.slotIntB(end,3,1) = solution.blockintegral(8);
                design.intBdata.slotIntB(end,4,1) = solution.blockintegral(9);
                
            else
                % x and y directed flux in slot on left hand side outer
                % (leftmost) layer
                mo_clearblock();
                mo_selectblock(slotxpos(i,1), slotypos(i,1));
                design.intBdata.slotIntB(end+1,1,1) = mo_blockintegral(8);
                design.intBdata.slotIntB(end,2,1) = mo_blockintegral(9);
                
                % x and y directed flux in slot on left hand side inner
                % (rightmost) layer
                mo_clearblock();
                mo_selectblock(slotxpos(i,2), slotypos(i,2));
                design.intBdata.slotIntB(end,3,1) = mo_blockintegral(8);
                design.intBdata.slotIntB(end,4,1) = mo_blockintegral(9);
            end
        end

    end
    
    % store the relative coil/slot positions, we use -ve slotypos as the
    % direction of sampling is the opposite of the direction of the fea
    % drawing, so choosing a slot in the +ve y direction is the same as
    % the magnets being in the opposite direction
    [thetaslotypos, ~] = cart2pol(slotxpos(:,1), slotypos(:,1));
    design.intBdata.slotPos = [design.intBdata.slotPos; 
                               (-thetaslotypos./design.thetap) + design.FirstSlotCenter + feapos];
    
    % sort the data in ascending position order
    [design.intBdata.slotPos, idx] = sort(design.intBdata.slotPos);
    design.intBdata.slotIntB = design.intBdata.slotIntB(idx,:,:);
    
end


function design = slotintAdata(design, simoptions, feapos, solution)

    % extract the flux integral data from all the slots at the given
    % positions of the magnet relative to coil
    if design.CoilLayers == 1
        
        slotypos = design.coillabellocs(1:(size(design.coillabellocs,1)),2);
        slotxpos = design.coillabellocs(1:(size(design.coillabellocs,1)),1);

        for i = 1:numel(slotypos)
            if exist('fpproc_interface_mex', 'file') == 3 && ~simoptions.usefemm
                % vector potential in slot on left hand side
                solution.clearblock();
                solution.selectblock(slotxpos(i,1), slotypos(i));
                design.intAdata.slotIntA(end+1,1,1) = solution.blockintegral(1);
            else
                % vector potential in slot on left hand side
                mo_clearblock();
                mo_selectblock(slotxpos(i,1), slotypos(i));
                design.intAdata.slotIntA(end+1,1,1) = mo_blockintegral(1);
            end
        end
        
    elseif design.CoilLayers == 2
        
        slotypos = [design.coillabellocs(1:2:(size(design.coillabellocs,1)),2), ...
                    design.coillabellocs((1:2:(size(design.coillabellocs,1)))+1,2)];
        slotxpos = [design.coillabellocs(1:2:(size(design.coillabellocs,1)),1), ...
                    design.coillabellocs((1:2:(size(design.coillabellocs,1)))+1,1)];
              
        for i = 1:size(slotypos,1)

            if exist('fpproc_interface_mex', 'file') == 3 && ~simoptions.usefemm
                % vector potential in slot on left hand side outer
                % (leftmost) layer
                solution.clearblock();
                solution.selectblock(slotxpos(i,1), slotypos(i,1));
                design.intAdata.slotIntA(end+1,1,1) = solution.blockintegral(1);
                
                % vector potential flux in slot on left hand side inner
                % (rightmost) layer
                solution.clearblock();
                solution.selectblock(slotxpos(i,2), slotypos(i,2));
                design.intAdata.slotIntA(end,2,1) = solution.blockintegral(1);
                
            else
                % vector potential in slot on left hand side outer
                % (leftmost) layer
                mo_clearblock();
                mo_selectblock(slotxpos(i,1), slotypos(i,1));
                design.intAdata.slotIntA(end+1,1,1) = mo_blockintegral(1);
                
                % vector potential in slot on left hand side inner
                % (rightmost) layer
                mo_clearblock();
                mo_selectblock(slotxpos(i,2), slotypos(i,2));
                design.intAdata.slotIntA(end,2,1) = mo_blockintegral(1);
                
            end
        end

    end
    
    % store the relative coil/slot positions. We use -ve slotypos as the
    % direction of sampling is the opposite of the direction of the fea
    % drawing, so choosing a slot in the +ve y direction is the same as
    % the magnets being in the opposite direction
    [thetaslotypos, ~] = cart2pol(slotxpos(:,1), slotypos(:,1));
    design.intAdata.slotPos = [design.intAdata.slotPos; 
                              (-thetaslotypos./design.thetap) + design.FirstSlotCenter + feapos];
    
    % sort the data in ascending position order
    [design.intAdata.slotPos, idx] = sort(design.intAdata.slotPos);
    design.intAdata.slotIntA = design.intAdata.slotIntA(idx,:,:);
    
end



function design = corelosssetup(design, feapos, solution)
    
    % get the number of positions
    npos = numel(feapos);
    
    if isa(solution, 'fpproc')
        % it's an xfemm fpproc object

        % get the volume of each element under consideration
        design.CoreLoss(1).dV = solution.getgroupareas (2) .* design.ls;
        
        % get the location we will use to estimate the flux in the element
        temp = solution.getgroupcentroids (2);
        design.CoreLoss(1).meshx = temp(:,1);
        design.CoreLoss(1).meshy = temp(:,2);
        
        % The values along the first dimension (down the columns) of the
        % arrays will contain values sampled in the x-direction in the fea
        % simulation. In this case values along the dimension hbi, or ht in
        % the case of the teeth. The values along the second dimension are
        % the values sampled at each value xRVWp. The values along the
        % third dimension of the arrays will contain values which are
        % sampled from the y direction of the fea simulation, in this case
        % along the dimension taupm, or tsb and tc in the case of the
        % teeth.
        design.CoreLoss(1).Bx = zeros([size(design.CoreLoss(1).meshx,1), npos, size(design.CoreLoss(1).meshx,2)]);
        design.CoreLoss(1).By = zeros([size(design.CoreLoss(1).meshx,1), npos, size(design.CoreLoss(1).meshx,2)]);
        design.CoreLoss(1).Bz = zeros([size(design.CoreLoss(1).meshx,1), npos, size(design.CoreLoss(1).meshx,2)]);

        % add xstep which is the size of the steps in xRVWp denormalised. This
        % is later used to find the value of dB/dx at each position
        design.CoreLoss(1).xstep = (feapos(2) - feapos(1)) * design.thetap;
        
    else
        error ('FEMM not currently supported for this, use xfemm.')
    end

end
