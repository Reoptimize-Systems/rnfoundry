function [T, Y, results, design, simoptions] = designandevaluate_TM(design, simoptions)

    [T, Y, results, design, simoptions] = evalsim_linear(design, simoptions);
    
    % if simoptions.Evaluation.mlength is not supplied we simulate a single pole of the
    % machine
    if isempty(simoptions.Evaluation.mlength)
        simoptions.Evaluation.mlength = [design.PoleWidth design.PoleWidth];
    end
    
    if ~isfield(design, 'nBpoints')
        design.nBpoints = 1;
    end
    
    if isfield(design, 'minLongMemberPoles')
        if isfield(simoptions, 'StatorPoles')
            design.Poles(simoptions.StatorPoles) = design.minLongMemberPoles;
        else
            % by default the armature is the stator which is stored in
            % design.Poles(2) for the ACTIAM
            design.Poles(2) = design.minLongMemberPoles;
        end
    end
    
    if isfield(design, 'bearingwidth')
        design.extraBpoles = ceil(design.nBpoints * design.bearingwidth / design.PoleWidth);
    else
        design.extraBpoles = 1 * design.nBpoints;
    end
    
    if design.Poles(1) == 1 && design.Poles(2) == 1                                         
        if simoptions.Evaluation.targetPower == 0
            design.Poles(1) = ceil(simoptions.Evaluation.mlength(1)/design.PoleWidth);
            design.PowerLoadMean = design.PowerLoadMean * design.Poles(1);
            design.Poles(2) = ceil(simoptions.Evaluation.mlength(2)/design.PoleWidth);
            
            if isfield(design, 'bearingwidth')
                design.extraBpoles = ceil(design.nBpoints * design.bearingwidth / design.PoleWidth);
            else
                design.extraBpoles = 1 * design.nBpoints; 
            end
        else
            design.Poles(1) = ceil(simoptions.Evaluation.targetPower / design.PowerLoadMean);
            design.PowerLoadMean = design.PowerLoadMean * design.Poles(1);
            % if target power specified, simoptions.Evaluation.mlength is a scalar containing the
            % overlap required (in m) between the field and armature, i.e.
            % how much longer the armature is than the field
            design.Poles(2) = ceil(((design.Poles(1)*design.PoleWidth)+simoptions.Evaluation.mlength)/design.PoleWidth);
          
            % Now add in extra Poles required to accomodate a bearing while
            % maintaining the same power
            if isfield(design, 'bearingwidth')
                design.extraBpoles = ceil(design.nBpoints * design.bearingwidth / design.PoleWidth);
            else
                design.extraBpoles = 1 * design.nBpoints; 
            end
            
            design.Poles = design.Poles + design.extraBpoles;
            
        end
    else
        if numel(design.Poles) == 1 && design.Poles(1) == design.PowerPoles && numel(simoptions.Evaluation.mlength) == 1
            % if only the number of field Poles, which are also the number
            % of power Poles is supplied, and an overlap length in
            % simoptions.Evaluation.mlength is suppplied, the overlap required (in m)
            % between the field and armature, is calculated in Poles and
            % assigned to design.Poles(2)
            design.Poles(2) = ceil(((design.Poles(1)*design.PoleWidth)+simoptions.Evaluation.mlength)/design.PoleWidth);
        end
        
        if isfield(design, 'bearingwidth')
            design.extraBpoles = ceil(design.nBpoints * design.bearingwidth / design.PoleWidth);
        else
            design.extraBpoles = 1 * design.nBpoints;
        end
    end
    
    % Now design the machine structure
    design.fLength = design.Poles(1) * design.PoleWidth;
    design.aLength = design.Poles(2) * design.PoleWidth;
    
    design.totalLength = [(design.Poles(2) - design.extraBpoles) * design.PoleWidth, ...
                          (design.Poles(2) - design.extraBpoles) * design.PoleWidth] ./ (design.nBpoints+1);
    
    % determine the
    if  design.PoleWidth * design.Poles(1) < design.totalLength(1)
        % the length of the field magnts and discs is less than the spacing
        % between bearing points
        fsupp = round2((design.totalLength(1) - (design.PoleWidth * design.Poles(1))) / 2, design.totalLength(1)/10000);
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