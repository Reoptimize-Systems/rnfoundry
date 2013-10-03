function avBvals = meanBincoilsections_TM(Bvals, Sr, Sz)
% meanBincoilsections_TM: A function to find the average B in each coil
% section for all positions using the B values at each point as reported by
% the GetFlux function.
%
% Arguments: (input)
%       
%       Bvals - (j x i) Matrix of Values of B as given by GetFlux where the columns
%       denote 'i' different coil positions and the values in those columns are
%       the 'j' B values at the points of the coil at the corners of all
%       sections. 
%
% Arguments: (output)
%
%       avBvals - (n x i) Matrix of average values of B for the sections of
%       coil at each of i positions. Sections are numbered from the top
%       left hand corner moving across the coil rowwise.

    avBvals = zeros((Sr*Sz), size(Bvals, 2)); % Preallocate for speed

    for i = 1:size(Bvals, 2)

        k = 1;
        for j = 1:(size(Bvals,1)-(Sz + 1))

            if mod(j, (Sz+1)) == 0
                % Skip if j counter is equal to the number of sections in order
                % to move to the next row of sections in coil.
            else
                % Otherwise store the average of the appropriate four B values
                % located at the four points of the section in the avBvals
                % matrix.
                avBvals(k,i) = (Bvals(j,i) + Bvals(j + 1,i) + Bvals(Sz + j + 1,i) + Bvals(Sz + 2 + j,i)) / 4;
                k = k + 1;
                
            end
            
        end

    end

end