function varargout = coilgapclosingforce_TM(varargin)
% Calculate the air-gap closing force in the tubular permanent magnet
% machine.
%
% Syntax
%
% [Force, avFro, avFri] = coilgapclosingforce_TM(I, avAxBvals, RoVRm, N, Rm, g, div, avFinfo)
% 
% [Force, avFro, avFri] = coilgapclosingforce_TM(design, Jz)
% 
% [Force, avFro, avFri] = coilgapclosingforce_TM(...,y)
% 
% [Force, avFro, avFri] = coilgapclosingforce_TM(...,xshift)
%     
% [Force, avFro, avFri] = coilgapclosingforce_TM(...,div)
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

    if nargin >= 7

        I = varargin{1};
        avAxBvals = varargin{2};
        RoVRm = varargin{3};
        N = varargin{4};
        Rm = varargin{5};
        g = varargin{6};
        div = varargin{7};        

        Stot = size(avAxBvals, 1);

        length = wirelength_TM(RoVRm, N, Rm, g);

        % Get length of one turn
        length = length / N;

        % get length of wire to be considered for force calculation
        delLength = length / div;

        if nargin == 8
        
            avFinfo = varargin{8};
            
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

    elseif nargin <= 5

        design = varargin{1};
        Jz = varargin{2};

        if nargin < 3
            y = 0;
        else
            y = varargin{3};
        end

        if nargin < 4
            xshift = 0;
        else
            xshift = varargin{4};
        end

        if nargin < 5
            div = 1;
        else
            div = varargin{5};
        end

        % get length of wire to be considered for force calculation
        delLength = design.MTL / div;

%         ymin = y - (design.WcVWp / 2);
%         ymax = y + (design.WcVWp / 2);
% 
%         xmin = design.Rm + design.g + xshift;
%         xmax = xmin + design.Hc;
%         
%         % if we are at the limit of the data in the x direction set xmax to
%         % this limit
%         if xmax > max(design.X(:))
%             xmax = max(design.X(:));
%         end
%         
%         if xmin > max(design.X(:))
%             xmin = max(design.X(:));
%         end
%         
%         if xmin == xmax
%             xmin = 0.99999 * xmax;
%         end
% 
%         tol = 1e-8;
        
        intBy = slmeval(xshift, design.slm_rposintBy, 0, false);
%         intBy = integratehalfperiod2ddata(design.X, design.Y, design.By, xmin, ymin, xmax, ymax, tol, 1) * design.Wp;

        Force =  intBy .* Jz .* delLength;

        varargout{1} = Force;

        if nargout > 1

            avFro = quad(@(z) lineBy(design.Ro, z, design.p_By), ymin, ymax) * design.Wp .* Jz  .* delLength ./ (pi * 2 * design.Ro * design.Wc );

            avFri = quad(@(z) lineBy(design.Ri, z, design.p_By), ymin, ymax) * design.Wp .* Jz .* delLength ./ (pi * 2 * design.Ri * design.Wc );

            varargout{2} = avFro;

            varargout{3} = avFri;
        end

    else
        error('Incorrect number of arguments.')

    end

    

end

function By = lineBy(xtemp, ytemp, p_By)

    xtemp = repmat(xtemp, size(ytemp));

    By = polyvaln(p_By, [xtemp', ytemp'])';

end