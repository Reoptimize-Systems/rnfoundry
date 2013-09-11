function sfdslm = makesfdeddyslm(rho, length, dc, nturns, pos, meanB1, meanB2, nstrands)
% creates an slm object which is fitted to precalculated values which can
% be used to calculate winding losses due to edddy currents in a coil due
% to the changing externally applied fields using the SFD method.
%
% Syntax
%
% sfdslm = makesfdeddyslm(rho, length, dc, pos, nturns, meanB)
% sfdslm = makesfdeddyslm(rho, length, dc, pos, nturns, meanB1, meanB2)
% sfdslm = makesfdeddyslm(rho, length, dc, pos, nturns, meanB1, meanB2, nstrands)
%
% Description
%
% makesvdeddyslm produces an slm object fitted to the part calculation of
% winding eddy current losses in a coil due to an externally applied
% time-varying field for use by the function lossforces_AM.m. The resulting
% calculation of losses is based on the Squared Field Derivative (SFD)
% method presented in [1]. The SFD method uses the formula:
%                              
%         pi * l * N * d_c^4   
%     P = ------------------ . < ( dB / dt )^2 > 
%              64 rho_c         
%                              
% Where l is the length of a strand in a winding in the changing field, N
% is the product of the number of strands and number of turns, d_c is the
% diameter of each strand, rho_c is the resistivity of the wire material,
% and < ( dB / dt )^2 >  is the spacial average of the value of (dB/dt)^2
% over the winding cross-section, to calculate the instantaneous power
% losses in the coil. lossforces_AM exploits the fact that for electrical
% machines, the value of dB / dt at any time is closely related to the
% motion of the machine such that it is equivalent to (dB / dx)*(dx/dt) =
% (dB / dx)*(v). Therefore, for a machine winding the losses can be
% calculated using the formula:
%
%               pi * l * N * d_c^4   
%     P = v^2 * ------------------ . < ( dB / dx )^2 > 
%                   64 rho_c         
%                        
% makesvdeddyslm produces a periodic SLM object fitted to the value of the
% part of this expression to the left of v^2 at positions normalised to the
% pole width (i.e. fitted to values of the expression between normalised x
% values between 0 and 2), such that multiplying the value returned by
% evaluating the slm by v^2 yields the instantaneous losses in a winding
% coil due to eddy currents in the conductors.
%
% Input
%
%   rho - the resistivity of the wire material
% 
%   length - the length of the piece of wire interacting with the magnetic
%     field
% 
%   dc - the diamter of the wire (if the wire is stranded this should be
%    the diameter of the bundle of wire, not the individual strands, the
%    diamter of each strand is calculated from the area of each strand,
%    where this area is calculated as the wire area / nstrands).
% 
%   pos - a vector of magnet positions at which values of the applied
%     magnetic field int coil parts have been sampled. For compatibility
%     with lossforces_AM.m, these positions should cover fully two poles of
%     displacement.
% 
%   nturns - the number of turns in the coil
%
%   meanB1 - either a vector or a matrix of two column vectors containing
%     values of the mean value of the flux density in the cross-section of
%     a coil. If two columns, each colum contains the values of orthogonal
%     components of the applied magnetic field, e.g. the Bx and By. 
% 
%   meanB2 - a second optional vector or matrix of two column vectors
%     containing values of the mean value of the flux density in the
%     cross-section of a second part of a coil. This is expected to be used
%     for coils which have two parts of equal size interacting with a
%     field, displaced from each other by the coil pitch.
% 
%   nstrands - optional, for stranded wire, the number strands per turn,
%     defaults to one, a single strand per turn.
%
% Output
%
%   svdslm - a periodic slm object fitted to the part calculation of losses
%     in a coil
%
%
% [1] C. R. Sullivan, "Computationally efficient winding loss calculation
% with multiple windings, arbitrary waveforms, and two-dimensional or
% three-dimensional field geometry," IEEE Transactions on Power
% Electronics, vol. 16, no. 1, pp. 142–150, 2001.
%
%
% See also: lossforces_AM
%

% Copyright Richard Crozier 2012


    if nargin < 6
        error('Insufficient input arguments.');
    end
    
    if nargin < 8
        nstrands = 1;
    end
    
    if isvector(meanB1)
        % if meanB1 is a vector, convert it to a column vector
        meanB1 = meanB1(:);
    end
    
    if size(meanB1,2) > 2
        error('meanB1 is wrong size, must be a vector or a matrix of one or two column vectors');
    end
    
    % get equivalent diameter for strands of stranded wire
    wirearea = pi * (dc/2)^2;
    strandarea = wirearea / nstrands;
    stranddiameter = 2 * area2radius(strandarea);
    
    % fit slms to the mean flux density data 
    
    % first component
    meanBx1slm = slmengine(pos, meanB1(:,1), ...
            'knots', round(numel(meanB1(:,1)) / 2), ...
            'EndCon', 'periodic', ...
            'Plot', 'off');
    
    if size(meanB1, 2) == 2
        % other component, if supplied
        meanBy1slm = slmengine(pos, meanB1(:,2), ...
                'knots', round(numel(meanB1(:,2)) / 2), ...
                'EndCon', 'periodic', ...
                'Plot', 'off');
    else
        % fit constant slm which evaluates to zero at all positions
        meanBy1slm = slmengine([pos(1), pos(end)], [0,0], ...
                'knots', 2, ...
                'Degree', 0, ...
                'EndCon', 'periodic', ...
                'Plot', 'off');
    end
        
    % make new positions at which the derivative of the slms fitted to the
    % meanB data will be evaluated
    nderivpos = 100;
    derivpos = linspace(pos(1), pos(end), nderivpos);
    
    % evaluate the field derivative for the first coil part
    dBx1dx = periodicslmeval(derivpos, meanBx1slm, 1);
    dBy1dx = periodicslmeval(derivpos, meanBy1slm, 1);  
    
    % if supplied also create slms for the 
    if nargin > 7 && ~isempty(meanB2)
        
        if isvector(meanB2)
            % if meanB1 is a vector, convert it to a column vector
            meanB2 = meanB2(:);
        end

        if size(meanB2,2) > 2
            error('meanB2 is wrong size, must be a vector or a matrix of one or two column vectors');
        end
    
        % first component
        meanBx2slm = slmengine(pos, meanB2(:,1), ...
                'knots', max(2,round(numel(meanB2(:,1)) / 2)), ...
                'EndCon', 'periodic', ...
                'Plot', 'off');

        if size(meanB2, 2) == 2
            % other component if supplied
            meanBy2slm = slmengine(pos, meanB2(:,2), ...
                    'knots', max(2,round(numel(meanB2(:,2)) / 2)), ...
                    'EndCon', 'periodic', ...
                    'Plot', 'off');
        else
            % fit constant slm which evaluates to zero at all positions
            meanBy2slm = slmengine([pos(1), pos(end)], [0,0], ...
                    'knots', 2, ...
                    'Degree', 0, ...
                    'EndCon', 'periodic', ...
                    'Plot', 'off');
        end

        % get the derivative at the desired points
        dBx2dx = periodicslmeval(derivpos, meanBx2slm, 1);
        dBy2dx = periodicslmeval(derivpos, meanBy2slm, 1);
        
    else
        dBx2dx = zeros(size(dBx1dx));
        dBy2dx = zeros(size(dBy1dx));
    end
    
    % calculate the information necessary to reproduce the SVD funtion
    % later (the result of the fit to this data must be multiplied by the
    % square of the velocity to calculate the instantaneous eddy current
    % loss)
    SFDpart = sum([dBx1dx.^2; dBx2dx.^2; dBy1dx.^2; dBy2dx.^2], 1) ...
                .* (pi * stranddiameter.^4 * length ./ (64 * rho)) * nturns * nstrands;
    
    % fit an slm to this part calculation for ease of reconstruction later,
    % the positions in this slm are normalised to the interval [0,2] for
    % compatibility with lossforces_AM.m
    sfdslm = slmengine(linspace(0,2,nderivpos), SFDpart, ...
                       'knots', round(nderivpos/4), ...
                       'EndCon', 'periodic', ...
                       'minvalue', 0, ...
                       'Plot', 'off');
                   
end