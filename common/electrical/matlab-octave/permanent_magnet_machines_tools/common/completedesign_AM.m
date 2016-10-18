function design = completedesign_AM (design, simoptions)
% completes common design aspects to create a minimal design ready for
% simulation
%
% Syntax
%
% design = completedesign_AM (design)
% design = completedesign_AM (design, simoptions)
%
% Description
%
% completedesign_AM performs some design completion tasks common to all
% machines. Most notably it completes the winding specification based on a
% minimum spec provided in the design structure. The winding is described
% using the following variables:
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
%  CoilLayers - the number of layers in the coil slot
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
% a combination of some or all of the following variables:
%
%  Phases - The number of phases in the machine
%  Poles - The number of magnetic poles in the machine
%  NBasicWindings - the number of basic winding segments in the machine
%  qc - number of coils per pole and phase (as a fraction object)
%  Qc - total number of coils (in all phases) in the machine
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
% further information. If CoilLayers is not supplied, it defaults to 1.
%
% completedesign_AM adds various default design options common to all pm
% machines if they are not already present in the design structure. The
% following fields will be added to 'design' if they are not already present:
%
%  CoilLayers - default is 1
%
%  MagnetSkew - determines the amount of magnet skewwing as a ratio of a pole
%   width (i.e. it is expected to be between 0 and 1). Defaults to zero if not
%   supplied.
%
%  NStrands - number of strands making up the wire in the coils. Defaults to 1
%   if not supplied.
%
%  NStages - number of stages making up the machine. Defaults to 1 if not
%   supplied
%
% Not all machine types may make use of all the default options set here.
%

    % set the magnet skew to zero by default
    design = setfieldifabsent (design, 'MagnetSkew', 0);

    % use single layered coils by default
    design = setfieldifabsent (design, 'CoilLayers', 1);
    
    % no coil insulation by default
    design = setfieldifabsent (design, 'CoilInsulationThickness', 0);
    
    % use single-stranded conductors by default
    design = setfieldifabsent (design, 'NStrands', 1);
    
    % check the number of stages in the design
    design = setfieldifabsent (design, 'NStages', 1);
    
    % check the minimum set of winding specification
    if all (isfield (design, {'Qc', 'Phases', 'Poles'})) && ~isfield (design, 'qc')
        % create qc fraction
        design.qc = fr (design.Qc, design.Poles * design.Phases);
    end

    % the number of Poles should be supplied and the ratio of coils to
    % Poles and Phases
    if ~( xor( all (isfield (design, {'qc', 'Phases', 'NBasicWindings'})), all (isfield (design, {'qc', 'Phases', 'Poles'})) )  ...
          || all (isfield (design, {'qc', 'Phases', 'Poles', 'NBasicWindings'})))
    
        error ( ['You must specify either the fields ''qc'', ''Phases'' and ''NBasicWindings''', ...
                 ' or ''qc'', ''Phases'' and ''Poles'', or all of these for the winding design specification'] )
            
    end
    
    if ~isa (design.qc, 'fr')
        error ('LINEAR:qcnotfr', ...
               'The qc winding spec variable must be an object of the ''fr'' (fraction) class.');
    end
        
    % get the "basic winding" (the smallest repetitive segment of the
    % machine winding in terms of symmetry) based on the qc ratio (the
    % ratio of coils per pole and phase). The smallest repetitive segment
    % is the smallest part that can be modelled using symmetric boundaries
    [design.Qcb,design.pb] = rat (design.qc * design.Phases);
    
    if ~isfield (design, 'Poles') && isfield (design, 'NBasicWindings')
        % multiply the basic winding by the desired number to get the full
        % winding, this therefore sizes ths machine, also determining the
        % number of Poles, this option is useful for optimisation routines
        
        % ensure a magnetically neutral rotor by having an even number of 
        % poles
        if design.pb == 1 && ~iseven (design.NBasicWindings)
            design.NBasicWindings = design.NBasicWindings + 1;
        end
    
        design.Qc = design.Qcb * design.NBasicWindings;
        design.Poles = design.pb * design.NBasicWindings;
        
    else        
         % determine the number of coils 
        [design.Qc,~] = rat (design.qc * design.Phases * design.Poles);
        
        % calculate the number of basic windings in the design
        design.NBasicWindings = design.Poles / design.pb;
        
    end
    
    % determine the total number of slots in the machine
    if design.CoilLayers == 2
        design.Qs = design.Qc;
        design.Qsb = design.Qcb;
        [coillayout, phaselayout] = windinglayout (design.Phases, design.Qs, design.Poles, 0);
    elseif design.CoilLayers == 1
        design.Qs = 2 * design.Qc;
        design.Qsb = 2 * design.Qcb;
        [coillayout, phaselayout] = windinglayout (design.Phases, design.Qs, design.Poles, 1);
    else
        error ('Only coils with one or two layers are implemented.')
    end
    
    design.WindingLayout = struct ('Coils', coillayout, 'Phases', phaselayout);
        
    % get the numerator and denominator of qc
    [design.qcn,design.qcd] = rat (design.qc);
    
    % Average coil pitch as defined by (Qs/Poles)
    design.yp = fr(design.Qs, design.Poles);
    
    % get the numerator and denominator of the coil pitch in slots
    [design.ypn,design.ypd] = rat (design.yp);
    
    % calculate the actual coil pitch in slots if not supplied
    if ~isfield(design, 'yd')
        if design.ypd == 1
            % the coil pitch in slots will be the same as the numerator of
            % yp, being an integral slot winding
            design.yd = design.ypn;
        else
            error('You must specify the coil pitch in design.yd for fractional slot windings.')
        end
    end
    
    design.NCoilsPerPhase = design.Qc / design.Phases;
    
    if design.pb > design.Poles
        error ('LINEAR:badwinding', 'Impossible winding design, basic winding poles greater than total number of poles.')
    end
    
end