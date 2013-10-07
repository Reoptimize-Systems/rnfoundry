function design = completedesign_ROTARY(design, simoptions)
% performs design creation operations common to all rotary type machines
%
% Syntax
%
% design = completedesign_ROTARY(design, simoptions)
%
% 
    

    % the number of Poles should be supplied and the ratio of coils to
    % Poles and Phases
    if ~( xor( all(isfield(design, {'qc', 'Phases', 'NBasicWindings'})), all(isfield(design, {'qc', 'Phases', 'Poles'})) )  ...
        || ~all(isfield(design, {'qc', 'Phases', 'Poles', 'NBasicWindings'})))
    
        error( ['You must specify either the fields ''qc'', ''Phases'' and ''NBasicWindings''', ...
                ' or ''qc'', ''Phases'' and ''Poles'', or all of these for the winding design specification'] )
            
    end
        
    % get the "basic winding" (the smallest repetitive segment of the
    % machine winding in terms of symmetry) based on the qc ratio (the
    % ratio of coils per pole and phase). The smallest repetitive segment
    % is the smallest part that can be modelled using symmetric boudaries
    [design.Qcb,design.pb] = rat(design.qc * design.Phases);
    
    % ensure a magnetically neutral rotor
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
            design.Qc = design.Qc / 2;
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
    
end


