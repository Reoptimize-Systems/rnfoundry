function evaloptions = designandevaloptions_FM(evaloptions)

    if nargin == 0
        evaloptions = [];
    end
    
    evaloptions = designandevaloptions_linear(evaloptions);
    
    if isempty(evaloptions)

        evaloptions.gfactor = 0.9;
        evaloptions.sections = 10;
        evaloptions.alphab = 1;
        evaloptions.IMethod = '1.3';
        evaloptions.InitDef = [0;0;0];
        evaloptions.maxPoleSupportBeams = inf;
        
    else

        % If the evaloptions structure is supplied, we check each possible
        % optional values and set any not supplied to the default value
                
        if ~isfield(evaloptions, 'gfactor')
            evaloptions.gfactor = 0.9;
        elseif evaloptions.gfactor >= 1.0
            error('gfactor must be less than 1.0 as structure will always deflect')
        end
        
        if ~isfield(evaloptions, 'sections')
            evaloptions.sections = 10;
        end
        
        if ~isfield(evaloptions, 'alphab')
            evaloptions.alphab = 1;
        elseif evaloptions.alphab < 1.0
            error('alphab (extra bearing length factor) must be greater than or equal to one')
        end
        
        % Beam type used to stiffen structure, default to hollow
        % rectangular section, see MomentOfInertiaY1 for other evaloptions
        if ~isfield(evaloptions, 'IMethod')
            evaloptions.IMethod = '1.3';
        end
        
        if ~isfield(evaloptions, 'InitDef')
            evaloptions.InitDef = [0;0;0];
        end
        
        evaloptions = setfieldifabsent(evaloptions, 'maxPoleSupportBeams', inf);
        
    end

end