function [BCoordinates, r] = GetCoilCoordinates(RoVRm, cwWp, Sr, Sz, pos, Wp, g)
% GetCoilCoordinates: a function to calculate the positions of the
% coordinates above the translator which will be required to determine the
% average flux density in a region of coil.
%
% Arguments: (input)
%
%   RoVRm - scalar value of Rm/Ro Ratio for machine to be evaluated, in
%   order to define the coil height
%
%   chWp - scalar value of Ratio of coil width to pole width (should 
%   normally be 1/3 but in the interests of generalisation will make 
%   other values up to ch/Wp = 1 possible) for machine to be evaluated
%        
%   Sr - number of sections in r direction into which coil is to be 
%   split (must be even)
%        
%   Sz - number of sections in z direction into which coil is to be 
%   split (must be even)
%
%   pos - array of positions of centre of coil relative to centre of a
%   steel piece (radial North position), defined as actual pos / Wp
%
%   g - Scalar normalised value of distance of coils above the translator
%   surface, i.e. actual g/Rm.
%
% Arguments: (output)
%
%   BCoordinates - (n x 3) array of r,z coordinates for calculating
%   the average flux at several coil positions. The rows of the first
%   column of the array are the r coordinates and the rows of the
%   second column are the corresponding z coordinates. Each matrix in
%   the third dimension corresponds to each position in pos.
%
% Copyright 2007 Richard Crozier and The Institute For Energy Systems at
% The University of Edinburgh

    % The coil is to be split up into a number of rectangles, or a grid,
    % the points of which will be described in relation to the translator.
    % This is to be done for each of the positions in the vector pos, which
    % defines the position of the centre of the coil in relation to the
    % centre of a North facing steel piece.
    
    % First check the coil is not wider than a pole, this is not currently
    % allowed.
    if cwWp > 1
       error(['The coil is wider than a pole width, this is not currently '...
           'supported by the program, please restate coil width and start again.']) 
    end
    
    % The position of the top of the coil is simply the Rm/Ro ratio. These
    % positions can be assigned to a vector of r position variables.
    
    r(1) = 1/RoVRm;
    r(Sr + 1) = g;
    diff = ((1/RoVRm) - g) / Sr;
    
    for n = 2:Sr
        
        r(n) = r(n-1) - diff;
                 
    end
    
    % Next work out Wp positions. These are first calculated for a coil
    % centred over a steel piece and will be modified  using the supplied
    % position data and normalised. 

    z = zeros(Sz + 1); % Preallocate dimensions for speed
    
    cw = cwWp * Wp; 
    
    z(1) = cw / 2;
    
    for n = 2:(Sz + 1)
        
       z(n) = z(n - 1) - (cw / Sz);
       
    end
    
    % Next we combine the r and z vectors into an r,z matrix of coordinates
    % which can then be modified into a 3 dimensional array of coordinates
    % for all the positions supplied in pos. The coordinates in the matrix
    % will be compiled from the top left of the coil moving across and then
    % downward for each line of points as shown below (example shows Sr = 3, 
    % Sz = 5).
    %          1   2   3   4   5   6
    %          .___.___.___.___.___.
    %         |7   8   9  10  11  12|
    %         |.   .   .   .   .   .|
    %         |13 14  15  16  17 -->|
    %         |.   .   .   .   .   .|
    %         |.___.___.___.___.___.|
    %
    %_____________________________________________
    
    n = 1;
    
    temp = zeros(2, ((Sz + 1) * (Sr + 1)) ); % Preallocate dimensions for speed
    
    for i = 1:(Sr + 1)

        for j = 1:(Sz +1)
            temp(1, n) = r(i);
            temp(2, n) = z(j);
            n = n + 1;
        end
        
    end
    
    % Now we can modify the z coordinates to create an array of temp
    % matrices for each of the coil positions given in pos. To do this,
    % generate a matrix of z values for each coil position and check that
    % coil position does not produce z coordinates outside the normalised
    % design space. 
    
    % Next recombine the coordinates to give the appropriate output
    
    BCoordinates = zeros(2, ((Sz + 1) * (Sr + 1)), size(pos, 2)); % Preallocate dimensions for speed
    
    for i = 1:size(pos, 2)
        
        BCoordinates(1, :, i) = temp(1, :);
        
        BCoordinates(2, :, i) = relativepoleposition_TM(temp(2,:) + pos(i), Wp);
        
        %BCoordinates(2, :, i) = relativepoleposition_TM(BCoordinates(2, :, i),1);
        
    end
    
    % Lastly, fix r to give the positions of the centre of sections for
    % GetFluxLinkage
    temp = [];
    for i = 1:(size(r, 2)-1)
        
       temp(i) = ((r(i) - r(i+1)) / 2 ) + r(i+1);
        
    end
    
    r = temp;

end