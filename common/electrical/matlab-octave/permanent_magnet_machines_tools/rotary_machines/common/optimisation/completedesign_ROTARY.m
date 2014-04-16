function design = completedesign_ROTARY(design, simoptions)
% performs design creation operations common to all rotary type machines
%
% Syntax
%
% design = completedesign_ROTARY(design, simoptions)
%
% 
% Description
%
% 

% Stator slot/coil/winding terminology is based on that presented in:
%
% J. J. Germishuizen and M. J. Kamper, "Classification of symmetrical
% non-overlapping three-phase windings," in The XIX International
% Conference on Electrical Machines - ICEM 2010, 2010, pp. 1-6.
% Single layer
%   q = 2qc
%   Qs = 2Qc
% Double layer
%   q = qc
%   Qs = Qc
%
% yp - Average coil pitch as defined by (Qs/Poles)
% yd - Actual coil pitch as defined by round(yp) +/- k
% Qs  -  total number of stator slots in machine
% Qc  -  total number of winding coils in machine 
% q  -  number of slots per pole and phase, (1)
% qn  -  numerator of q
% qd  -  denominator of q
% qc - number of coils per pole and phase, (2)
% qcn  -  numerator of qc
% qcd  -  denominator of qc

    % check the minimum set of winding specification
    if all(isfield(design, {'Qc', 'Phases', 'Poles'})) && ~isfield (design, 'qc')
        % create qc fraction
        design.qc = fr(design.Qc, design.Poles * design.Phases);
    end

    % the number of Poles should be supplied and the ratio of coils to
    % Poles and Phases
    if ~( xor( all(isfield(design, {'qc', 'Phases', 'NBasicWindings'})), all(isfield(design, {'qc', 'Phases', 'Poles'})) )  ...
          || all(isfield(design, {'qc', 'Phases', 'Poles', 'NBasicWindings'})))
    
        error( ['You must specify either the fields ''qc'', ''Phases'' and ''NBasicWindings''', ...
                ' or ''qc'', ''Phases'' and ''Poles'', or all of these for the winding design specification'] )
            
    end
        
    % get the "basic winding" (the smallest repetitive segment of the
    % machine winding in terms of symmetry) based on the qc ratio (the
    % ratio of coils per pole and phase). The smallest repetitive segment
    % is the smallest part that can be modelled using symmetric boundaries
    [design.Qcb,design.pb] = rat(design.qc * design.Phases);
    
    % ensure a magnetically neutral rotor by having at least 2 poles in the
    % basic winding segment
    if design.pb == 1
        design.pb = 2;
        design.Qcb = 2 * design.Qcb;
    end
    
    if ~isfield(design, 'Poles') && isfield(design, 'NBasicWindings')
        % multiply the basic winding by the desired number to get the full
        % winding, this therefore sizes ths machine, also determining the
        % number of Poles, this options is useful for optimisation routines
        design.Qc = design.Qcb *  design.NBasicWindings;
        design.Poles = design.pb *  design.NBasicWindings;
        
        % determine the total number of slots in the machine
        if design.CoilLayers == 2
            design.Qs = design.Qc;
        elseif design.CoilLayers == 1
            design.Qs = 2 * design.Qc;
        else
            error('Only coils with one or two layers are implemented.')
        end
        
    else        
         % determine the number of slots 
        [design.Qs,~] = rat(design.qc * design.Phases * design.Poles);
        
        % calculate the number of basic windings in the design
        design.NBasicWindings = design.Poles / design.pb;
        
        % determine the total number of coils in the machine
        if design.CoilLayers == 2
            design.Qc = design.Qs;
        elseif design.CoilLayers == 1
            design.Qc = design.Qs / 2;
        else
            error('Only coils with one or two layers are implemented.')
        end
        
    end
    
    % get the numerator and denominator of qc
    [design.qcn,design.qcd] = rat(design.qc);
    
    % Average coil pitch as defined by (Qs/Poles)
    design.yp = fr(design.Qs, design.Poles);
    
    % get the numerator and denominator of the coil pitch in slots
    [design.ypn,design.ypd] = rat(design.yp);
    
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
    
    design.thetap = 2*pi / design.Poles;
    
    design.NCoilsPerPhase = design.Qc / design.Phases;
    
end


