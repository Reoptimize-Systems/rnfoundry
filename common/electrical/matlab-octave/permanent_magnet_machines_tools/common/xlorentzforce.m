function Fx = xlorentzforce(Jz, p_intBy, indepvar, wirelen)

    if iscell(p_intBy)

        Fx = zeros(size(p_intBy));

        if isscalar(Jz)

            Jz = repmat(Jz, size(p_intBy));
            
        end

        if all(size(Jz)==size(p_intBy))

            for k = 1:size(p_intBy,3)
                for j = 1:size(p_intBy,2)
                    for i = 1:size(p_intBy,1)

                        Fx(i,j,k) = ylorentzforce(Jz(i,j,k), p_intBy{i,j,k}, indepvar{i,j,k}, wirelen);

                    end
                end
            end

        else
            error('Jz must be a scalar or matrix of same dimensions as cell array of polynomials')
        end

    elseif isstruct(p_intBy)

        % Multiply integral of flux density by current density and mean
        % turn length to get forces
        Fx = polyvaln(p_intBy, indepvar) .* Jz .* wirelen;

    elseif ismatnotvec(p_intBy)

        error('matpoly form not yet supported')

    end


end

