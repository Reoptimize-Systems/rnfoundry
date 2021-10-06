function varargout = coilclosingforce_TM(I, avAxBvals, RoVRm, N, Rm, g, div, avFinfo)
% Calculates the air-gap closing force in a tubular machine coil.
%
% Input:
%
%   I - (n x 1) column vector of current values
%
%   avAxBvals - (n x k) matrix of flux density values. These are the
%               average flux densities in n coil sections at k positions.  
%
%   RoVRm - Ro / Rm ratio of machine
%
%   N - Number of turns in coil
%
%   Rm - The translator radius, can be in any metric units provided the
%        value of the air-gap is given in the same units.
%
%   g - The length of the air-gap, must be in the same units as Rm.
%
%   div - This is the number of sections into which the coil is split,
%         usually used for calculating the net force when the coil is
%         displaced from centre
%
%   avFinfo - optional (1 x 3) vector giving info required for the
%             calculation of the per unit volume forces at the inside and
%             outside layers of the coil for calculation of stresses. 
%        
%           avFinfo(1) - Sz, the number of sections per row
%           avFinfo(2) - Sr, the number of sections per column
%           avFinfo(3) - Wp, the pole width
%
% Output:
%
%   varargout - if avFinfo is not present, this is a row vector of the
%   total forces due to the B field in each column of avAxBvals and the
%   current in the coils at the corresponding postion.
%
%   If avFinfo is present, there are three outputs, the first is as described
%   previously, the second is a vector containing the per unit volume
%   forces at the outer coil sections at each position considered and the
%   third is the per unit volume forces at the inner coil sections
%   
    Stot = size(avAxBvals, 1);

    length = wirelength_TM(RoVRm, N, Rm, g);
    
    % Get length of one turn
    length = length / N; 
    
    % get length of wire to be considered for force calculation
    delLength = length / div; 

    if nargin == 8
        
       Ro = RoVRm * Rm;
       % find the section height by dividing the coil height by the number
       % of radial sections, passed in in avFinfo(2)
       sHeight = (Ro - Rm - g) / avFinfo(2);
       % Find the volume of a row of coil sections, (same for all rows)
       % using the cross-sectional area of the coil section and the pole
       % width stored in avFinfo(3)
       sRowVolume = (avFinfo(3)/3) * pi * (Ro^2 - (Ro-sHeight)^2);
       
    end
    
    % Preallocate vectors with zeros
    Force = zeros(size(avAxBvals,2),1);
    avFro = Force;
    avFri = Force;
    
    for n = 1:size(avAxBvals,2)
        
        % Get the force in each section for position 'n'
        temp = avAxBvals(:,n) .* I .* delLength .* (N/Stot);
        % Total the forces
        Force(n) = sum(temp);
        
        if nargin == 8

            % Reshape the temp force vector to easily extract the desired
            % section forces, the columns of this matrix will each be a
            % layer of sections of the coil working from the outside in as
            % we move across the matrix columns 
            temp = reshape(temp(:,n), avFinfo(1), []);
            
            % calculate the total force along the outside of the coil, we
            % will consider the total force, so scale back up by div
            avFro(n) = sum(temp(:,1).*div);
            
            % calculate the total force per uit volume along the inside
            % of the coil, again scale back up by div
            avFri(n) = sum(temp(:,end).*div);

        end
        
    end

    varargout{1} = Force;
    
    if nargin == 8
        
       varargout{2} = avFro ./ sRowVolume;
       
       varargout{3} = avFri ./ sRowVolume;
       
    end
    
end