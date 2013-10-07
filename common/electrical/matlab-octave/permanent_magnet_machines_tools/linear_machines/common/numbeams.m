function n_beams = numbeams(b, PoleWidth, fieldLen, bpfactor)
% numbeams: determines the number of beams that are appropriate for a
% design
%
% Input:
%
%   b - beam width
%
%   PoleWidth - width of one pole
%
%   armLen - total length of all the Poles
%
%   bpfactor - beam spread factor, defines the desired spacing between 
%              beams as a factor of the pole width. This is the fraction of
%              one pole width which lies between adjacent beams.
%
% Ouput:
%
%   n_beams - the number of beams appropriate for the desired design. A
%             minimum of one beam per pole is assumed. 


    % first determine the fraction of a pole width taken up by beam
    bpfrac = b / PoleWidth;
    
    % fraction of pole width which is empty
    spacefrac = 1 - bpfrac;
    
    % determine the desired space between beams
    desSpace = (bpfactor * PoleWidth) - b;
    
    if bpfrac >=0 && bpfrac < 1.0 && desSpace > 0

        % determine the available space between beams
        avBeamSpace = spacefrac * PoleWidth;
        
        if avBeamSpace < desSpace
            % in this case, space the beams by the available space, i.e.
            % just have one beam per pole
            n_beams = round(fieldLen / PoleWidth);
        else
            % in this case space out the beams by the desired amount
            totalb = b + desSpace;
            n_beams = ceil(fieldLen / totalb);
        end

    elseif bpfrac >= 1.0 || desSpace <= 0
        % the beam width is more than one pole width, or the desired space
        % between beams is less than or equal to zero so we just lay out
        % the beams flush agains one another with no regard for division by
        % electrical pole, this is the closest we will get to the minimum
        % spec of one beam per pole in this case
        n_beams = ceil(fieldLen / b);
    else
        % the beam width or pole must have been negative
        error('Beam width or pole must have been negative, check inputs')
    end

end