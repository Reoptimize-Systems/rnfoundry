function varargout = evaluatefestructure_PMSM(design, options, last_n_beams, BeamInfo)
% evaluatefestructure_PMSM: evaluates whether a given beam can withstand
% the forces in a linear permanent magnet synchronous machine design using
% the FAESOR finite element toolbox
%
%

    if nargin < 4
        BeamInfo = [];
    end

    if ~(isfield(design, 'fens') && isfield(design, 'gcells')) || (design.OuterStructureNBeams ~= last_n_beams)

        % create the finite element mesh for feAirGapClosure
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

        [design.fens, design.gcells, BeamInfo, design.OuterSupportCells] = createmeshfebeamdef2_FM(x, y, z, n, znodes, ...
            zguidecon, zsupp, zframe, zextreme, guideyshift, tolerance, BeamInfo);

    end
    
    if isempty(BeamInfo)

        [new_g, closingForces, BeamInfo] = feairgapclosure_PMSM(design, options);
        
    else

        [new_g, closingForces, BeamInfo] = feairgapclosure_PMSM(design, options, BeamInfo);

    end

    % Check if beam is successful or a failure
    if min(new_g(:)) >= options.gfactor * design.g
        varargout{1} = 1;
        varargout{2} = design;
        varargout{3} = BeamInfo;
    else
        varargout{1} = 0;
        varargout{2} = design;
        varargout{3} = BeamInfo;
    end

end