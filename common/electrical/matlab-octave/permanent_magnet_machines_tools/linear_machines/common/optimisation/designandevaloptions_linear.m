function evaloptions = designandevaloptions_linear(evaloptions)

    if nargin == 0
        evaloptions = [];
    end
    
    evaloptions = designandevaloptions_common(evaloptions);
    
    if isempty(evaloptions)

        evaloptions.targetPower = 500e3;
        evaloptions.mlength = 6;
        evaloptions.presimodeevfun = 'prescribedmotode_linear';
        evaloptions.presimresfun = 'prescribedmotresfun_linear';
        evaloptions.presimforcefun = 'forcefcn_linear_pscbmot';

    else

        % If the evaloptions structure is supplied, we check each possible
        % optional values and set any not supplied to the default value

        if ~isfield(evaloptions, 'targetPower')
            % if a target power is not supplied, we set this to zero and check
            % for a set machine length vector
            if ~isfield(evaloptions,'mlength')
                % if there's no machine length supplied either use a single
                % pole for the machine
                evaloptions.mlength = [];
            else
%                 % if a vector for the machine length is supplied but is the
%                 % wrong size raise an error
%                 if numel(evaloptions.mlength) ~= 2 || ~(numel(evaloptions.mlength) == 1 && design.Poles) 
%                     error('If a target power is not supplied in the evaloptions structure, evaloptions.mlength must be a vector of 2 values,\n the length of the field and armature respectively')
%                 end

            end
            evaloptions.targetPower = 0;
        else
            if isfield(evaloptions, 'mlength') && ~isempty(evaloptions.mlength)
                
                if evaloptions.targetPower == 0 && size(evaloptions.mlength,2) ~= 2
                    error('If a target power is supplied in the evaloptions structure but is zero, and mlength is also supplied,mlength must be a vector of 2 values, the length of the field and armature respectively')
                end

                if evaloptions.targetPower > 0 && size(evaloptions.mlength,2) ~= 1
                    error('If a target power greater than zero is supplied in the evaloptions structure, and mlength is also supplied,mlength must be a scalar value of the overlap of the field and armature')
                end
                
            else
                
                if evaloptions.targetPower > 0
                    % if there's no machine length supplied use a single
                    % pole for the machine
                    evaloptions.mlength = [];
                else
                    % Otherwise set mlength to zero to give no extra
                    % overlap between field and armature
                    evaloptions.mlength = 0;
                end
                
            end
        end
        
        evaloptions = setfieldifabsent(evaloptions, 'presimodeevfun' ,'prescribedmotode_linear');
        evaloptions = setfieldifabsent(evaloptions, 'presimresfun', 'prescribedmotresfun_linear');
        evaloptions = setfieldifabsent(evaloptions, 'presimforcefun', 'forcefcn_linear_pscbmot');

    end


end