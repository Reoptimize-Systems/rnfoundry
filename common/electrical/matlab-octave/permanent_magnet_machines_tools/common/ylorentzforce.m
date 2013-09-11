function Fy = ylorentzforce(Jz, intBxvar, indepvar, wirelen)
% calculates lorentz forces current density in a magnetic field, typically
% to find the reaction, or shear, force on the winding
%
% Syntax
%
% Fy = ylorentzforce(Jz, intBxvar, indepvar, wirelen)
%
% Input
%
%   Jz, 
%
%   intBxvar, 
%
%   indepvar, 
%
%   wirelen
%
% Output
%
%   

    if iscell(intBxvar)
        
        Fy = zeros(size(intBxvar));
        
        if isscalar(Jz)
            
            Jz = repmat(Jz, size(intBxvar));
            
        end
            
        if all(size(Jz)==size(intBxvar))
        
        for k = 1:size(intBxvar,3)
            for j = 1:size(intBxvar,2)
                for i = 1:size(intBxvar,1)
                    
                    Fy(i,j,k) = ylorentzforce(Jz(i,j,k), intBxvar{i,j,k}, indepvar{i,j,k}, wirelen);
                    
                end
            end
        end
        
        else
            error('Jz must be a scalar or matrix of same dimensions as cell array of polynomials')
        end
        
    elseif isstruct(intBxvar)
        
        if isslm(intBxvar)
            
            Fy = evalhpslmodd(intBxvar, 0, indepvar, 1) .* Jz .* wirelen;

        else

            % Multiply integral of flux density by current density and mean
            % turn length to get forces
            Fy = polyvaln(intBxvar, indepvar) .* Jz .* wirelen;

        end
        
        
    elseif ismatnotvec(intBxvar)
        
        error('matpoly form not yet supported')
        
    end


end

