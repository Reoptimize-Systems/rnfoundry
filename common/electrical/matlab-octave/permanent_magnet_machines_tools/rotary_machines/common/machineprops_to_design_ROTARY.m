function [design, simoptions] = machineprops_to_design_ROTARY (Rphase, Lphase, Poles, Qc, varargin)
% creates a rotary PM machine design structure from basic properties
%
% Syntax
%
% [design, simoptions] = machineprops_to_design_ROTARY (Rphase, Lphase, Poles, Qc)
% [design, simoptions] = machineprops_to_design_ROTARY (..., 'Parameter', Value)
%
% Description
%
% machineprops_to_design_ROTARY creates a machine design structure and
% simoptions structure from som basic properties of an electrical machine
% that is suitible for use with the rotary machine simulation functions. It
% is intended to create a simulation of a machine when detailed information
% or the real design of the machine is not available.
%
% Input
%
%  Rphase - Phase resistance of the machine
%
%  Lphase - Phase inductance of the machine
%
%  Poles - number of magnetic poles on the field
%
%  Qc - total number of coils in machine 
%
% Addtional arguments may be supplied as parameter-value pairs. As well as
% the required arguments You must supply either the option 'FluxPhase'
% or 'FluxCoil'. The available options are:
%
%  'FluxPhase' - The per-phase flux linkage. This can be either a scalar
%    value, or a vector of values. If a scalar this is taken as the peak
%    value of a pure sinusoidal flux linkage waveform. If a vector, this is
%    assumed to be the flux linkage seen by a phase at multiple evenly
%    spaced displacements of the translator over two poles of the machine
%    (i.e. one traversal of a pole-pair).
%
%  'FluxCoil' - The per-coil flux linkage. This can be either a scalar
%    value, or a vector of values. If a scalar this is taken as the peak
%    value of a pure sinusoidal flux linkage waveform. If a vector, this is
%    assumed to be the flux linkage seen by a coil at multiple evenly
%    spaced displacements of the translator over two poles of the machine
%    (i.e. one traversal of a pole-pair). The per-phase flux linkage is
%    determined by multiplying this value (or all values) by the number of
%    coils per phase.
%
%  'Mphase' - Mutual inductance between phases. Default is zero if not
%    supplied.
%
%  'CoilTurns' - Number of turns per coil. Default is 1 if not supplied.
%    Since the phase resistance and inductance is supplied directly it is
%    generally safe to simply leave this value as 1, regardless of what the
%    real number of turns is. This, in combination with the wire diameter
%
%  'CoilFillFactor' - The fraction of the cross-sectional area of the coil
%    which is made up of copper. Default is 0.55. This will be used to
%    calculate the wire diameter if a wire diameter is not supplied (using
%    the 'Dc' option below).
%
%  'Dc' - Diameter of wire in the coils (or equivalent diameter if wire is
%    not round or made up of multiple strands). If not supplied, this is
%    calculated from the number of turns and coil fill factor. It is
%    required to calculate the correct current density in the wire.
%
%  'NStrands' - Number of parallel strands of wire making up each turn in
%    the coil. Default is 1 if not supplied.
%
%  'Branches' - Number of parallel branches of of series coils in the
%    winding. Default is 1, meaning all coils are wired in series.
%
%  'Phases' - The number of phases in the machine. Default is 3.
%
%  'yd' - coil pitch in terms of number of slots
%
%  'CoilLayers' - Number of layers in the coil.
%
% Output
%
%  design - A complete machine design sturcture
%
%  simoptions - A simoptions structure with information on the options used
%   to set up the machine design. You should append your own simulation
%   options to this structure and pass it to the simulation functions along
%   with the design structure.
%
%


    options.Mphase = 0;
    options.CoilTurns = 1;
    options.NStrands = 1;
    options.Branches = 1;
    options.Dc = [];
    options.CoilFillFactor = 0.55;
    options.Phases = 3;
    options.yd = [];
    options.CoilLayers = 2;
    options.CoilTurns = 1;
    options.FluxPhase = [];
    options.FluxCoil = [];
    
    options = parse_pv_pairs (options, varargin);
    
    if isempty (options.FluxPhase) && isempty (options.FluxCoil)
        error ('You must supply either FluxPhase or FluxCoil');
    end
    
    if ~isempty (options.FluxPhase) && ~isempty (options.FluxCoil)
        error ('You have supplied values for both FluxPhase and FluxCoil, you must use only one of these options');
    end

    design.Qc = Qc;
    design.CoilTurns = options.CoilTurns;
    design.Poles = Poles;
    design.PhaseResistance = Rphase;
    design.PhaseInductance = [ Lphase, options.Mphase ];
    design.Phases = options.Phases;
    design.Branches = options.Branches;
    
    if ~isempty (options.yd)
        design.yd = options.yd;
    end
    
    if ~isempty (options.Dc)
        design.Dc = options.Dc;
    end
    
    design.CoilFillFactor = options.CoilFillFactor;
    
    simoptions = struct ();
    
    design = completedesign_ROTARY (design, simoptions);
    
    simoptions.LoadModel = 'none';
    design = circuitprops_AM (design, simoptions);
    simoptions = rmfield (simoptions, 'LoadModel');
    
    if isempty (options.FluxCoil)
        options.FluxCoil = options.FluxPhase ./ design.CoilsPerBranch;
    end

    if isscalar (options.FluxCoil)
        
        design.psilookup = linspace(0, 2, 30);
        design.psilookup(2,:) = options.FluxCoil .* sin (design.psilookup(1,:) .* pi);
                             
    elseif isvector (options.FluxCoil)
        
        design.psilookup = linspace(0, 2, numel(options.FluxCoil));
        design.psilookup(2,:) = options.FluxCoil;
        
    elseif ( size (options.FluxCoil, 1) == 2 ) ...
            && ( max (options.FluxCoil(1,:) == 2) ) ...
            && ( min (options.FluxCoil(1,:) == 0) )
        
        design.psilookup = options.FluxCoil;
        
    else
        error ('FluxPhase must be a scalar value, or a vector of flux values over two poles, or an matrix with 2 rows')
    end
    
    simoptions.SkipCheckCoilProps = true;
    simoptions.SkipFEA = true;
    
    [design, simoptions] = simfun_ROTARY (design, simoptions);
    
end