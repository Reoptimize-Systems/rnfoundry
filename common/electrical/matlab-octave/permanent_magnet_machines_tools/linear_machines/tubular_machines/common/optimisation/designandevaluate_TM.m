function [T, Y, results, design, simoptions] = designandevaluate_TM(design, simoptions)

    [T, Y, results, design, simoptions] = evalsim_linear(design, simoptions);
    
    % if simoptions.evaloptions.mlength is not supplied we simulate a single pole of the
    % machine
    if isempty(simoptions.evaloptions.mlength)
        simoptions.evaloptions.mlength = [design.PoleWidth design.PoleWidth];
    end
    
    if ~isfield(design, 'nBpoints')
        design.nBpoints = 1;
    end
    
    if isfield(design, 'minLongMemberPoles')
        if isfield(simoptions, 'StatorPoles')
            design.poles(simoptions.StatorPoles) = design.minLongMemberPoles;
        else
            % by default the armature is the stator which is stored in
            % design.poles(2) for the ACTIAM
            design.poles(2) = design.minLongMemberPoles;
        end
    end
    
    if isfield(design, 'bearingwidth')
        design.extraBpoles = ceil(design.nBpoints * design.bearingwidth / design.PoleWidth);
    else
        design.extraBpoles = 1 * design.nBpoints;
    end
    
    if design.poles(1) == 1 && design.poles(2) == 1                                         
        if simoptions.evaloptions.targetPower == 0;
            design.poles(1) = ceil(simoptions.evaloptions.mlength(1)/design.PoleWidth);
            design.PowerLoadMean = design.PowerLoadMean * design.poles(1);
            design.poles(2) = ceil(simoptions.evaloptions.mlength(2)/design.PoleWidth);
            
            if isfield(design, 'bearingwidth')
                design.extraBpoles = ceil(design.nBpoints * design.bearingwidth / design.PoleWidth);
            else
                design.extraBpoles = 1 * design.nBpoints; 
            end
        else
            design.poles(1) = ceil(simoptions.evaloptions.targetPower / design.PowerLoadMean);
            design.PowerLoadMean = design.PowerLoadMean * design.poles(1);
            % if target power specified, simoptions.evaloptions.mlength is a scalar containing the
            % overlap required (in m) between the field and armature, i.e.
            % how much longer the armature is than the field
            design.poles(2) = ceil(((design.poles(1)*design.PoleWidth)+simoptions.evaloptions.mlength)/design.PoleWidth);
          
            % Now add in extra poles required to accomodate a bearing while
            % maintaining the same power
            if isfield(design, 'bearingwidth')
                design.extraBpoles = ceil(design.nBpoints * design.bearingwidth / design.PoleWidth);
            else
                design.extraBpoles = 1 * design.nBpoints; 
            end
            
            design.poles = design.poles + design.extraBpoles;
            
        end
    else
        if numel(design.poles) == 1 && design.poles(1) == design.PowerPoles && numel(simoptions.evaloptions.mlength) == 1
            % if only the number of field poles, which are also the number
            % of power poles is supplied, and an overlap length in
            % simoptions.evaloptions.mlength is suppplied, the overlap required (in m)
            % between the field and armature, is calculated in poles and
            % assigned to design.poles(2)
            design.poles(2) = ceil(((design.poles(1)*design.PoleWidth)+simoptions.evaloptions.mlength)/design.PoleWidth);
        end
        
        if isfield(design, 'bearingwidth')
            design.extraBpoles = ceil(design.nBpoints * design.bearingwidth / design.PoleWidth);
        else
            design.extraBpoles = 1 * design.nBpoints;
        end
    end
    
    % Now design the machine structure
    design.fLength = design.poles(1) * design.PoleWidth;
    design.aLength = design.poles(2) * design.PoleWidth;
    
    design.totalLength = [(design.poles(2) - design.extraBpoles) * design.PoleWidth, ...
                          (design.poles(2) - design.extraBpoles) * design.PoleWidth] ./ (design.nBpoints+1);
    
    % determine the
    if  design.PoleWidth * design.poles(1) < design.totalLength(1)
        % the length of the field magnts and discs is less than the spacing
        % between bearing points
        fsupp = round2((design.totalLength(1) - (design.PoleWidth * design.poles(1))) / 2, design.totalLength(1)/10000);
        asupp = fsupp;
        
        design.supportLengths = [fsupp, fsupp;
                                 asupp, asupp];
    else
        % the length of the field magnets and discs is greater than the
        % spacing between bearing points
        design.supportLengths = [0, 0; 
                                 0, 0];
    end
    
end