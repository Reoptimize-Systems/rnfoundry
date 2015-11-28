function [FemmProblem, rotorinfo, statorinfo] = slotlessfemmprob_radial(design, varargin)
% creates a FemmProblem structure for a slotless torus axial flux permanent
% magnet machine
               
    Inputs.DrawingType = 'MagnetRotation';
    Inputs.NBoundaryPositions = 10;
    Inputs.BoundaryShift = 0;
    Inputs.ArmatureType = 'external';
    Inputs.NWindingLayers = nan;
    Inputs.CoilCurrent = zeros (1,design.Phases);
    Inputs.MagArrangement = 'NN';
    Inputs.PolarisationType = 'constant';
    if isfield (design, 'MagnetPolarisation') && ischar (design.MagnetPolarisation)
        Inputs.PolarisationType = design.MagnetPolarisation;
    end
    Inputs.FemmProblem = newproblem_mfemm ('planar', 'Depth', design.ls, 'MinAngle', 15);
    Inputs.Position = 0;
    Inputs.FractionalPolePosition = [];
    Inputs.RotorAnglePosition = [];
    Inputs.MagnetGroup = [];
    Inputs.MagnetSpaceGroup = [];
    Inputs.RotorBackIronGroup = [];
    Inputs.RotorOuterRegionGroup = [];
    Inputs.CoilGroup = 0;
    Inputs.ArmatureBackIronGroup = [];
    Inputs.MagnetRegionMeshSize = choosemesharea_mfemm (design.tm, (design.Rmm*design.thetam), 1/10);
    Inputs.BackIronRegionMeshSize = choosemesharea_mfemm (min(design.tbi), 2*(design.Rbm*design.thetap), 1/10);
    Inputs.RotorOuterRegionsMeshSize = [choosemesharea_mfemm(design.tm, (design.Rbo*design.thetap), 1/5), -1];
    Inputs.StatorOuterRegionsMeshSize = [];
    Inputs.StatorOuterRegionMaterials = {};
    Inputs.AirGapMeshSize = choosemesharea_mfemm (design.g, (design.Rmm*design.thetap), 1/10);
    Inputs.DrawOuterRegions = true;
    Inputs.StatorOuterRegionSize = [];
    Inputs.DrawCoilInsulation = false;
    Inputs.CoilInsRegionMeshSize = -1;
    
    if design.tsg > 1e-5
        if design.tsb > 1e-5
            Inputs.ShoeGapRegionMeshSize = ...
                choosemesharea_mfemm (max(design.tsg, design.tsb), (design.Rmo*design.thetasg), 1/20);
        else
            Inputs.ShoeGapRegionMeshSize = ...
                choosemesharea_mfemm (design.tsb, (design.Rmo*design.thetasg), 1/20);
        end
    else
        if design.tsb > 1e-5
            Inputs.ShoeGapRegionMeshSize = ...
                choosemesharea_mfemm (design.tsb, (design.Rmo*design.thetasg), 1/20);
        else
            Inputs.ShoeGapRegionMeshSize = -1;
        end
    end
    Inputs.YokeRegionMeshSize = mean( [choosemesharea_mfemm(design.ty, 2*(design.Rym*design.thetap), 1/10), ...
                                       choosemesharea_mfemm(design.tc(1), (design.Rcm*(design.thetas-mean(design.thetac))), 1/10)] );
    Inputs.CoilRegionMeshSize = choosemesharea_mfemm (design.tc(1), (design.Rcm*mean(design.thetac)));
    Inputs.Tol = 1e-5;
    Inputs.SimType = 'Magnetics';
    Inputs.MaterialsLibrary = '';
    Inputs.NPolePairs = 1;

    Inputs = parse_pv_pairs (Inputs, varargin);
    
    % we'll draw a radial flux stator with the stator iron replaced with
    % air
    
    switch Inputs.ArmatureType
        
        case 'external'
            % single inner facing stator
            Rs = design.Rmo + design.g + design.tc(1) + design.tsb + design.ty/2;
        case 'internal'
            % single outer facing stator
            Rs = design.Rmi - design.g - design.tc(1) - design.tsb - design.ty/2;
        case 'di'
            % double internal stator (mags on outside)
%             Rs = design.Rmo(1) + design.g + design.tc(1) + design.ty/2;
            error('not yet supported');
        case 'do'
            % double outer/external stator (mags on inside)
            error('not yet supported');
            
        otherwise
            error('Unrecognised ArmatureType option.')
                
    end
    
    % the sizes of the stator outer air regions
    if isempty (Inputs.StatorOuterRegionSize)
        Inputs.StatorOuterRegionSize = [design.ty, 2*design.tm, 10*design.tm];
    else
        Inputs.StatorOuterRegionSize = [design.ty, Inputs.StatorOuterRegionSize];
    end
    if isempty (Inputs.StatorOuterRegionsMeshSize)
        Inputs.StatorOuterRegionsMeshSize = [ Inputs.YokeRegionMeshSize, ...
                                              choosemesharea_mfemm(design.tm, (Rs*design.thetap), 1/5), ...
                                              repmat(-1,1,numel(Inputs.StatorOuterRegionSize)-2) ];
    else
        Inputs.StatorOuterRegionsMeshSize = [ Inputs.YokeRegionMeshSize, ...
                                              Inputs.StatorOuterRegionsMeshSize ];
    end
    assert (samesize (Inputs.StatorOuterRegionSize, Inputs.StatorOuterRegionsMeshSize), ...
        'RENEWNET:slotlessfemmproble_radial:nstatormeshsizes', ...
        'Number of supplied stator outer region mesh sizes does not match number of outer region sizes.');

    
    % handle the outer region materials
    % get what the armature yoke material should be and store it
    yokemagsimmat = design.MagFEASimMaterials.ArmatureYoke;
    % swap the yoke material to be the same as the air gap material
    design.MagFEASimMaterials.ArmatureYoke = design.MagFEASimMaterials.AirGap;
    
    if isempty (Inputs.StatorOuterRegionMaterials)
        % FemmProblem.Materials(GapMatInd).Name
        Inputs.StatorOuterRegionMaterials = [ yokemagsimmat, ...
            repmat({design.MagFEASimMaterials.AirGap}, 1, numel(Inputs.StatorOuterRegionSize)-1) ];
    else
        Inputs.StatorOuterRegionMaterials =  [{ yokemagsimmat }, Inputs.StatorOuterRegionMaterials ];
    end
    assert (samesize (Inputs.StatorOuterRegionSize, Inputs.StatorOuterRegionMaterials), ...
        'RENEWNET:slotlessfemmproble_radial:nstatormaterials', ...
        'Number of supplied stator outer region material names does not match number of outer region sizes.');
        
    
    % make the yoke drawing thickness as small as possible, we will be
    % using the first outer region as the yoke
    design.ty = 3 * Inputs.Tol;
    
    [FemmProblem, rotorinfo, statorinfo] = slottedfemmprob_radial( design, Inputs );
    
end

